program main

    use mpi
    use constants, only: file_list, record_unit, nx_max, nz_max, substrate_list, hc, e0
    use variable, only: max_files, snapshot_index, species, augmentation_index, atom_typ1, atom_typ2, placed, &
                        nx, ny, nz, ht, data_file, found, n_ranks, output_data, substrate_dir, &
                        lambda
    use utils_io
    use nano_process
    use random_utils
    use slc, only: run_slc
    use msa, only: run_msa
    use wavimg, only: run_wavimg
    use descriptor

    implicit none

    integer             ::  ierr, rank, first_file, last_file, accepted_local, accepted_total
    integer             ::  n_images_local, n_accepted
    character(len=255)  ::  rank_suffix, rank_records
    character(len=255)  ::  xyz_dir, env_var
    logical             ::  accepted

    call mpi_init(ierr)
    call mpi_comm_rank(MPI_COMM_WORLD, rank, ierr)
    call mpi_comm_size(MPI_COMM_WORLD, n_ranks, ierr)

    call get_environment_variable('xyz_dir', xyz_dir)
    call get_environment_variable('data_file', data_file)
    call get_environment_variable('substrate_dir', substrate_dir)
    call get_environment_variable('output_data', output_data)
    call get_environment_variable('max_files', env_var)
    read(env_var, *) max_files
    call get_environment_variable('atom_typ1', atom_typ1)
    call get_environment_variable('atom_typ2', atom_typ2)
    call get_environment_variable('n_px', env_var)
    read(env_var, *) nx
    ny = nx ! grille carree uniquement
    call get_environment_variable('nz', env_var)
    read(env_var, *) nz
    call get_environment_variable('ht', env_var)
    read(env_var, *) ht

    if (nx .lt. 1 .or. nx .gt. nx_max) stop 'error: n_px out of range [1, 256]'
    if (nz .lt. 1 .or. nz .gt. nz_max) stop 'error: nz out of range [1, 20]'
    if (len_trim(output_data) .eq. 0) stop 'error: output_data not set'


    if(rank.eq.0) call execute_command_line('find ' // xyz_dir // & 
            ' -maxdepth 1 -name "*.xyz" | sort -V > ' // file_list)
    if(rank.eq.0) call execute_command_line('find ' // trim(substrate_dir) // &
        ' -maxdepth 1 -name "*.xyz" -printf "%f\n" | sort -V > ' // substrate_list)

    if (rank .eq. 0) call execute_command_line('mkdir -p xyz_substrate')
    call mpi_barrier(MPI_COMM_WORLD, ierr)

    call read_file_list
    call read_substrate_list
    call load_data
    first_file = rank * max_files / n_ranks + 1
    last_file  = (rank + 1) * max_files / n_ranks
    write(rank_suffix, *) rank
    rank_records = 'records_rank_' // trim(adjustl(rank_suffix)) // '.tmp'

    open(record_unit, file=rank_records, action='write', status='replace', access='stream', form='unformatted')
    n_images_local = (last_file - first_file + 1) * 10
    accepted_local = 0
    lambda = hc / sqrt(ht * (2.0d0 * e0 + ht))
    do snapshot_index = first_file, last_file
        n_accepted = 0
        call read_xyz
        call read_data        
        if (.not. found) cycle
        call extract_largest_cluster
        if (count(species(1:n_atoms) .eq. atom_typ1) .eq. 0 .or. count(species(1:n_atoms) .ne. atom_typ1) .eq. 0) then
            write(*,*) "Rank ", rank, ": snapshot ", snapshot_index, " does not contain both atom types. Skipping."
            cycle
        endif
        call compute_descriptors
        do augmentation_index = 1, 10
            call stable_seed
            accepted = .false.
            call seed_random            
            call compute_rotation
            call place_in_box
            if (.not. placed) cycle
            call choose_substrate
            call read_substrate_file
            call add_substrate_to_cluster
            if (.not. placed) cycle
            call read_input
            call prepare_hrtem_particle
            
            call run_slc
            call run_msa
            call run_wavimg

            call save_data_row(accepted)
            if (accepted) n_accepted = n_accepted + 1
        enddo
        accepted_local = accepted_local + n_accepted
    enddo
    close(record_unit)

    call mpi_reduce(accepted_local, accepted_total, 1, mpi_integer, mpi_sum, 0, mpi_comm_world, ierr)
    call mpi_barrier(mpi_comm_world, ierr)

    if (rank .eq. 0) then
        call save_data
        write(*, '(A,1X,I0,1X,A,1X,I0,1X,A)') 'prep: ', accepted_total, 'accepted images from', max_files, 'snapshots'
        write(*, '(A,1X,A)') 'prep: wrote descriptors and images to', trim(output_data)
    endif

    call mpi_barrier(mpi_comm_world, ierr)
    call execute_command_line('rm -f ' // trim(rank_records), exitstat=ierr)
    if (rank .eq. 0) call execute_command_line('rm -f ' // file_list, exitstat=ierr)
    if (rank .eq. 0) call execute_command_line('rm -f ' // substrate_list, exitstat=ierr)


    call mpi_finalize(ierr)

endprogram