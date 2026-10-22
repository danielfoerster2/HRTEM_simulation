module utils_io

    implicit none
    
    integer             :: id_sim_tab(48000)
    integer             :: n_atoms_tab(48000), n_steps_tab(48000)
    double precision    :: composition_tab(48000), initial_temperature_tab(48000)

    ! Champs d'un enregistrement du fichier de sortie, dans l'ordre d'ecriture de save_data_row (voir README.md):
    ! les descripteurs entiers (int32), les descripteurs reels (float64), puis l'image (uint8, ny x nx)
    character(len=8), parameter     :: int_fields(9) = [character(len=8) :: 'id_sim', 'n_atoms', 'n_steps', &
                                        'nat1', 'nat2', 'nat1_out', 'nat2_out', 'nat1_in', 'nat2_in']
    character(len=19), parameter    :: real_fields(28) = [character(len=19) :: 'initial_temperature', 'epot_total', &
                                        'composition', 'gyration_radius', 'd_com', 'coreshell_index', 'phi', 'theta', 'psi', &
                                        'position_scale', 'defocus', 'A1_re', 'A1_im', 'B2_re', 'B2_im', 'A2_re', 'A2_im', &
                                        'C3', 'S3_re', 'S3_im', 'A3_re', 'A3_im', 'B4_re', 'B4_im', 'D4_re', 'D4_im', &
                                        'A4_re', 'A4_im']

    contains


    subroutine read_file_list

        use constants, only: file_list
        use variable, only: max_files, xyz_files

        integer :: ierr, i_file, n_files

        n_files = 0
        open(10, file=file_list, status='old', action='read')
        do i_file = 1, size(xyz_files)
            read(10, '(A)', iostat=ierr) xyz_files(i_file)
            if (ierr /= 0) exit
            n_files = n_files + 1
        end do
        close(10)
        if (n_files .eq. 0) stop 'error: no xyz files found'
        if (max_files .gt. 0) then
            max_files = min(max_files, n_files) 
        else
            max_files = n_files
        end if

    endsubroutine

    subroutine read_substrate_list
    
        use constants, only: substrate_list
        use variable, only: substrate_files, n_substrate_files

        integer             :: ierr, i_file, len_max
        character(len=255)  :: line

        ! Premier passage: nombre de fichiers et longueur maximale des noms, pour allouer la liste
        n_substrate_files = 0
        len_max = 1
        open(10, file=substrate_list, status='old', action='read')
        do
            read(10, '(A)', iostat=ierr) line
            if (ierr /= 0) exit
            n_substrate_files = n_substrate_files + 1
            len_max = max(len_max, len_trim(line))
        enddo
        if (n_substrate_files .eq. 0) stop 'error: no substrate files found'
        allocate(character(len=len_max) :: substrate_files(n_substrate_files))
        rewind(10)
        do i_file = 1, n_substrate_files
            read(10, '(A)') substrate_files(i_file)
        enddo
        close(10)

    endsubroutine


    subroutine load_data

        use variable, only: xyz_files, data_file

        integer :: i_file, ierr

        open(10, file=data_file, status='old', action='read')
        do i_file = 1, size(xyz_files)
            read(10, *, iostat=ierr) id_sim_tab(i_file), n_atoms_tab(i_file), composition_tab(i_file), n_steps_tab(i_file), &
                        initial_temperature_tab(i_file)
            if (ierr /= 0) exit
        enddo
        close(10)

    endsubroutine


    subroutine read_xyz

        use constants, only: n_atoms_max
        use variable, only: xyz_files, snapshot_index, pos, species, epot, box
        use descriptor, only: n_atoms

        integer             :: i_atom, index0, index1
        double precision    :: lattice(9)
        character(len=255)  :: line

        open(10, file=xyz_files(snapshot_index), status='old', action='read')
        read(10, *) n_atoms
        if (n_atoms .gt. n_atoms_max) then
            write(*,*) 'error: n_atoms = ', n_atoms, ' exceeds n_atoms_max = ', n_atoms_max, &
                ' in ', trim(xyz_files(snapshot_index))
            stop
        endif
        read(10, '(A)') line
        do i_atom = 1, n_atoms
            read(10, *) species(i_atom)(1:2), pos(1, i_atom), pos(2, i_atom), pos(3, i_atom), epot(i_atom)
        end do
        close(10)
        index0 = index(line, 'Lattice="') + len('Lattice="')
        index1 = index(line(index0:), '"')
        index1 = index1 + index0 - 2
        read(line(index0:index1), *) lattice
        box = (/ lattice(1), lattice(5), lattice(9) /)
        pos(1, :) = modulo(pos(1, :), box(1))
        pos(2, :) = modulo(pos(2, :), box(2))
        pos(3, :) = modulo(pos(3, :), box(3))
        
    endsubroutine


    subroutine read_data
        
        use variable, only: xyz_files, snapshot_index, found
        use descriptor, only: id_sim, n_atoms, composition, n_steps, initial_temperature

        integer             :: i_file, idx
        character(len=255)  :: filename, basename
        integer             :: id_from_file

        filename = xyz_files(snapshot_index)
        idx = index(filename, '/', back=.true.)
        if (idx .gt. 0) then
            basename = filename(idx+1:)
        else
            basename = filename
        end if
        read(basename(1:index(basename, '.xyz')-1), *) id_from_file
        found = .false.
        do i_file = 1, size(xyz_files)
            if (id_sim_tab(i_file) .eq. id_from_file) then
                if (n_atoms_tab(i_file) .ne. n_atoms) then
                    write(*,*) "error: n_atoms mismatch for snapshot ", snapshot_index, ": n_atoms in xyz file = ", n_atoms, &
                        ", n_atoms in data file = ", n_atoms_tab(i_file), " (id_sim = ", id_sim_tab(i_file), ")"
                else
                    found = .true.
                    id_sim = id_sim_tab(i_file)
                    n_atoms = n_atoms_tab(i_file)
                    composition = composition_tab(i_file)
                    n_steps = n_steps_tab(i_file)
                    initial_temperature = initial_temperature_tab(i_file)
                end if
                exit
            end if
        enddo

    endsubroutine


    subroutine read_input

        use constants, only: box_hrtem
        use variable, only: fs, edge, sc_mrad, vib1, vib2, vibdir, oapr, dose_e_per_a2, readout_noise_e, &
                            doptc, dopsc, dovib, aberr_re, aberr_im, lambda, pos_cluster
        use random_utils, only: random_uniform, sample_aberration
        use descriptor, only: n_atoms


        doptc = 1 ! focal spread (1 = on, 0 = off)
        fs = 0.7d0 ! ecart-type de focal spread (nm)
        dopsc = 1 ! Active enveloppe de cohérence spatiale partielle (envs) (1 = on, 0 = off)
        dovib = 2 ! Mode vibration (1 = isotrope, 2 = anisotrope)
        edge = 0.0d0 ! Largeur du bord de transition de l'ouverture objective

        sc_mrad = random_uniform(0.05d0, 0.25d0) ! Semi angle de convergence de la source (en mrad)
        vib1 = random_uniform(0.005d0, 0.04d0) ! Amplitude de vibration direction 1 (en nm)
        vib2 = random_uniform(0.005d0, 0.04d0) ! Amplitude de vibration direction 2 (en nm)
        vibdir = random_uniform(0.0d0, 180.0d0) ! Direction de vibration (en degrés)
        oapr = random_uniform(20.0d0, 30.0d0) ! Rayon de l'ouverture objective (en mrad)
        dose_e_per_a2 = random_uniform(3000.0d0, 15000.0d0) ! Dose électronique par unité de surface (en e-/A^2) (intensité bruit de poisson)
        readout_noise_e = random_uniform(0.0d0, 1.0d0) ! bruit de lecture du détecteur

        aberr_re = 0.0d0
        aberr_im = 0.0d0
        
        call sample_aberration(3, 0.0d0, 6.0d0) ! A1 2-fold astigmatism
        call sample_aberration(4, 0.0d0, 50.0d0) ! B2 Axial coma
        call sample_aberration(5, 0.0d0, 50.0d0) ! A2 3-fold astigmatism
        call sample_aberration(6, 5000.0d0, 20000.0d0) ! C3 Spherical aberration (Cs)
        
        call sample_aberration(7, 0.0d0, 700.0d0) ! S3 Star aberration
        call sample_aberration(8, 0.0d0, 700.0d0) ! A3 4-fold astigmatism
        call sample_aberration(9, 0.0d0, 1500.0d0) ! B4 5th-order term
        call sample_aberration(10, 0.0d0, 1500.0d0) ! D4 5th-order term
        call sample_aberration(11, 0.0d0, 1500.0d0) ! A4 5th-order term

        ! Defocus de Scherzer + U(-3, 3) nm, rapporte au plan median de la particule (z moyen, nul apres place_in_box)
        ! et non au plan de sortie de la boite: on retire la distance entre ce plan et la sortie de la boite
        aberr_re(2) = -1.2d0 * sqrt(abs(aberr_re(6)) * lambda) + random_uniform(-3.0d0, 3.0d0)
        aberr_re(2) = aberr_re(2) - (0.5d0*box_hrtem(3)*10.0d0 - sum(pos_cluster(3,1:n_atoms)) / dble(n_atoms)) * 0.1d0
        aberr_im(2) = 0.0d0

    endsubroutine


    subroutine save_data

        use, intrinsic :: iso_fortran_env, only: int64
        use variable, only: n_ranks, output_data, nx, ny

        integer                         :: i_file, ierr
        integer(int64)                  :: file_size, total_size, record_size
        character(len=10)               :: rank_suffix
        character(len=:), allocatable   :: records_tmp

        ! Les fichiers par rang contiennent deja les enregistrements finaux: en-tete NPY puis concatenation dans l'ordre des rangs
        record_size = int(4 * size(int_fields) + 8 * size(real_fields) + nx * ny, int64)
        records_tmp = ''
        total_size = 0
        do i_file = 0, n_ranks-1
            write(rank_suffix, '(I0)') i_file
            records_tmp = records_tmp // ' records_rank_' // trim(rank_suffix) // '.tmp'
            inquire(file='records_rank_' // trim(rank_suffix) // '.tmp', size=file_size)
            total_size = total_size + file_size
        enddo
        if (mod(total_size, record_size) .ne. 0) stop 'error: truncated records'
        call write_npy_header(trim(output_data), int(total_size / record_size))
        call execute_command_line('cat' // records_tmp // ' >> ' // trim(output_data), exitstat=ierr)
        if (ierr .ne. 0) stop 'error: could not write output_data'

    endsubroutine


    subroutine write_npy_header(filename, n_images)

        use variable, only: nx, ny

        integer, intent(in)             :: n_images
        character(len=*), intent(in)    :: filename
        integer                         :: unit_npy, len_header, i
        character(len=16)               :: s_nx, s_ny, s_n
        character(len=:), allocatable   :: descr, dict

        ! En-tete NPY 1.0 d'un tableau structure (int_fields, real_fields, image): numpy.load(filename, mmap_mode='r')
        write(s_nx, '(I0)') nx
        write(s_ny, '(I0)') ny
        write(s_n, '(I0)') n_images
        descr = '['
        do i = 1, size(int_fields)
            descr = descr // "('" // trim(int_fields(i)) // "', '<i4'), "
        enddo
        do i = 1, size(real_fields)
            descr = descr // "('" // trim(real_fields(i)) // "', '<f8'), "
        enddo
        descr = descr // "('image', '|u1', (" // trim(s_ny) // ", " // trim(s_nx) // "))]"
        dict = "{'descr': " // descr // ", 'fortran_order': False, 'shape': (" // trim(s_n) // ",), }"
        ! magic (6) + version (2) + longueur (2) + dictionnaire complete par des blancs et '\n': multiple de 64 octets
        len_header = 64 * ((10 + len(dict) + 1 + 63) / 64) - 10
        dict = dict // repeat(' ', len_header - len(dict) - 1) // achar(10)
        open(newunit=unit_npy, file=filename, access='stream', form='unformatted', status='replace')
        write(unit_npy) char(147) // 'NUMPY' // char(1) // char(0) // char(mod(len_header, 256)) // char(len_header / 256), dict
        close(unit_npy)

    endsubroutine


    subroutine save_data_row(written)

        use, intrinsic :: iso_fortran_env, only: int8, int32, real64
        use constants, only: record_unit
        use descriptor
        use variable, only: nx, ny, image, aberr_re, aberr_im, phi, theta, psi, position_scale

        logical, intent(out)    :: written
        integer             :: i_px, j_px, q
        integer(int8)       :: qimage(nx, ny)
        integer(int32)      :: int_values(size(int_fields))
        real(real64)        :: real_values(size(real_fields))
        double precision    :: scaled, image_min, image_max

        ! Un enregistrement binaire par image: descripteurs (int_fields, real_fields) puis l'image en octets non signes
        ! 0-255 (x le plus rapide), ecrits en une seule instruction (voir README.md)
        written = .false.
        image_min = minval(image(1:nx, 1:ny))
        image_max = maxval(image(1:nx, 1:ny))
        if (image_max .le. image_min) return
        do j_px = 1, ny
            do i_px = 1, nx
                scaled = 255.0d0 * (image(i_px, j_px) - image_min) / (image_max - image_min)
                q = max(0, min(255, nint(scaled)))
                qimage(i_px, j_px) = int(q - 256 * (q / 128), int8) ! octet non signe stocke dans un entier 8 bits
            enddo
        enddo
        ! Meme ordre que int_fields et real_fields
        int_values = [id_sim_bis, n_atoms, n_steps, nat1, nat2, nat1_out, nat2_out, nat1_in, nat2_in]
        real_values = [initial_temperature, epot_total, composition, gyration_radius, d_com, coreshell_index, &
            phi, theta, psi, position_scale, aberr_re(2), aberr_re(3), aberr_im(3), aberr_re(4), aberr_im(4), &
            aberr_re(5), aberr_im(5), aberr_re(6), aberr_re(7), aberr_im(7), aberr_re(8), aberr_im(8), &
            aberr_re(9), aberr_im(9), aberr_re(10), aberr_im(10), aberr_re(11), aberr_im(11)]
        write(record_unit) int_values, real_values, qimage
        written = .true.

    endsubroutine


endmodule