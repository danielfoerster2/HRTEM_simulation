module nano_process
    
    implicit none
    
    contains

    subroutine extract_largest_cluster
        
        use descriptor, only: n_clusters, n_atoms
        use constants, only: cluster_cutoff
        use variable, only: pos, box, species, epot

        integer             ::  i_atom, j_atom, stack_top, current_atom, largest_cluster, i_cluster, largest_size
        integer             ::  cluster_id(n_atoms), cluster_size(n_atoms), stack(n_atoms)
        double precision    ::  cutoff2, delta(3), pos_tmp(3, n_atoms), epot_tmp(n_atoms)
        logical             ::  visited(n_atoms)
        character(len=2)    ::  species_tmp(n_atoms)

        cutoff2 = cluster_cutoff * cluster_cutoff
        cluster_id = 0
        cluster_size = 0
        visited = .false.
        pos_tmp = 0.0d0
        n_clusters = 0

        do i_atom = 1, n_atoms
            if (visited(i_atom)) cycle
            n_clusters = n_clusters + 1
            stack_top = 1
            stack(1) = i_atom
            visited(i_atom) = .true.
            cluster_id(i_atom) = n_clusters
            pos_tmp(:, i_atom) = pos(:, i_atom)

            do while(stack_top .gt. 0)
                current_atom = stack(stack_top)
                stack_top = stack_top - 1
                cluster_size(n_clusters) = cluster_size(n_clusters) + 1
            
                do j_atom = 1, n_atoms
                    if (visited(j_atom)) cycle
                    delta = pos(:, j_atom) - pos(:, current_atom)
                    delta = delta - box * nint(delta / box)
                    if(dot_product(delta, delta) .gt. cutoff2) cycle
                    stack_top = stack_top + 1
                    stack(stack_top) = j_atom
                    visited(j_atom) = .true.
                    cluster_id(j_atom) = n_clusters
                    !pos_tmp(:, j_atom) = pos(:, j_atom)
                    pos_tmp(:, j_atom) = pos_tmp(:, current_atom) + delta
                enddo
            enddo
        enddo

        largest_cluster = 1
        largest_size = cluster_size(1)
        do i_cluster = 2, n_clusters
            if (cluster_size(i_cluster) .gt. largest_size) then
                largest_cluster = i_cluster
                largest_size = cluster_size(i_cluster)
            endif
        enddo
        epot_tmp = epot(1:n_atoms)
        species_tmp = species(1:n_atoms)
        i_cluster = 0
        pos = 0.0d0
        species = '  '
        epot = 0.0d0
        do i_atom = 1, n_atoms
            if (cluster_id(i_atom) .ne. largest_cluster) cycle
            i_cluster = i_cluster + 1
            species(i_cluster) = species_tmp(i_atom)
            pos(:, i_cluster) = pos_tmp(:, i_atom)
            epot(i_cluster) = epot_tmp(i_atom)
        enddo
        n_atoms = largest_size
        
        ! A verifier que ce calcul centre bien le cluster sur la boite de simulation
        !pos(1, 1:n_atoms) = pos(1, 1:n_atoms) - box(1) * nint((pos(1, 1:n_atoms) - pos(1, 1) )/ box(1))
        !pos(2, 1:n_atoms) = pos(2, 1:n_atoms) - box(2) * nint((pos(2, 1:n_atoms) - pos(2, 1) )/ box(2))
        !pos(3, 1:n_atoms) = pos(3, 1:n_atoms) - box(3) * nint((pos(3, 1:n_atoms) - pos(3, 1) )/ box(3))
        pos(1, 1:n_atoms) = pos(1, 1:n_atoms) - sum(pos(1, 1:n_atoms)) / n_atoms
        pos(2, 1:n_atoms) = pos(2, 1:n_atoms) - sum(pos(2, 1:n_atoms)) / n_atoms
        pos(3, 1:n_atoms) = pos(3, 1:n_atoms) - sum(pos(3, 1:n_atoms)) / n_atoms
        !pos(1, 1:n_atoms) = pos(1, 1:n_atoms) - box(1) * nint(pos(1, 1:n_atoms) / box(1))
        !pos(2, 1:n_atoms) = pos(2, 1:n_atoms) - box(2) * nint(pos(2, 1:n_atoms) / box(2))
        !pos(3, 1:n_atoms) = pos(3, 1:n_atoms) - box(3) * nint(pos(3, 1:n_atoms) / box(3))
        !pos(1, 1:n_atoms) = pos(1, 1:n_atoms) - sum(pos(1, 1:n_atoms)) / n_atoms
        !pos(2, 1:n_atoms) = pos(2, 1:n_atoms) - sum(pos(2, 1:n_atoms)) / n_atoms
        !pos(3, 1:n_atoms) = pos(3, 1:n_atoms) - sum(pos(3, 1:n_atoms)) / n_atoms
        !pos(1, 1:n_atoms) = pos(1, 1:n_atoms) - box(1) / 2.0d0
        !pos(2, 1:n_atoms) = pos(2, 1:n_atoms) - box(2) / 2.0d0
        !pos(3, 1:n_atoms) = pos(3, 1:n_atoms) - box(3) / 2.0d0

    endsubroutine
    

    subroutine compute_descriptors
        
        use, intrinsic :: ieee_arithmetic
        use descriptor
        use variable, only: species, pos, atom_typ1, epot
        
        integer             ::  i_atom
        double precision    ::  dist2, radius_limit, term1, term2, term3, radius_limit2
        
        nat1 = 0
        nat2 = 0
        nat1_out = 0
        nat1_in = 0
        nat2_out = 0
        nat2_in = 0

        mass_center = 0.0d0
        mass_center_nat1 = 0.0d0
        mass_center_nat2 = 0.0d0

        do i_atom = 1, n_atoms
            mass_center = mass_center + pos(:, i_atom)
            if (species(i_atom) .eq. atom_typ1) then
                nat1 = nat1 + 1
                mass_center_nat1 = mass_center_nat1 + pos(:, i_atom)
            else
                nat2 = nat2 + 1
                mass_center_nat2 = mass_center_nat2 + pos(:, i_atom)
            end if
        enddo
        mass_center = mass_center / dble(n_atoms)
        if (nat1 .gt. 0) mass_center_nat1 = mass_center_nat1 / dble(nat1)
        if (nat2 .gt. 0) mass_center_nat2 = mass_center_nat2 / dble(nat2)
        epot_total = sum(epot(1:n_atoms))

        gyration_radius = 0.0d0
        do i_atom = 1, n_atoms
            gyration_radius = gyration_radius + dot_product(pos(:, i_atom) - mass_center, &
                    pos(:, i_atom) - mass_center)
        enddo
        gyration_radius = sqrt(gyration_radius / dble(n_atoms))
        radius_limit = 0.8d0 * gyration_radius
        radius_limit2 = radius_limit * radius_limit
        do i_atom = 1, n_atoms
            dist2 = dot_product(pos(:, i_atom) - mass_center, pos(:, i_atom) - mass_center)
            if(dist2 .ge. radius_limit2) then
                if(species(i_atom) .eq. atom_typ1) then
                    nat1_out = nat1_out + 1
                else
                    nat2_out = nat2_out + 1
                endif
            else
                if(species(i_atom) .eq. atom_typ1) then
                    nat1_in = nat1_in + 1
                else
                    nat2_in = nat2_in + 1
                endif
            endif
        enddo
        if (nat1 .gt. 0 .and. nat2 .gt. 0) then
            d_com = sqrt(dot_product(mass_center_nat1 - mass_center_nat2, mass_center_nat1 - mass_center_nat2))
        else
            d_com = ieee_value(d_com, ieee_quiet_nan)
        endif
        if(gyration_radius .gt. 0.0d0 .and. nat1 .gt. 0 .and. nat2 .gt. 0 .and. &
                nat1_out + nat2_out .gt. 0 .and. nat1_in + nat2_in .gt. 0) then
            term1 = dble(nat1_out) / dble(nat1_out + nat2_out) - composition
            term2 = dble(nat1_in) / dble(nat1_in + nat2_in) - composition
            term3 = d_com / (2.0d0 * gyration_radius)
            coreshell_index = 2.0d0 * abs(term1) + 2.0d0 * abs(term2) - term3
        else
            coreshell_index = ieee_value(coreshell_index, ieee_quiet_nan)
        endif

    endsubroutine


    subroutine compute_rotation
    
        use descriptor, only: n_atoms, id_sim, id_sim_bis, mass_center
        use variable, only: pos, pos_cluster, augmentation_index, phi, theta, psi, position_scale
        use constants, only: pi

        integer             ::  i_atom
        double precision    ::  u, rot_matrix(3, 3), cp, sp, ct, st, cs, ss

        do i_atom = 1, n_atoms
            pos_cluster(:, i_atom) = pos(:, i_atom) - mass_center
        enddo
        call random_number(u)
        position_scale = 0.95d0 + 0.1d0 * u
        pos_cluster = pos_cluster * position_scale
        ! Rotation uniforme sur SO(3): R = Rz(phi) Ry(theta) Rz(psi), avec cos(theta) uniforme sur [-1, 1]
        call random_number(phi)
        call random_number(theta)
        call random_number(psi)
        phi = 2.0d0 * pi * phi
        theta = acos(1.0d0 - 2.0d0 * theta)
        psi = 2.0d0 * pi * psi

        cp = cos(phi)
        sp = sin(phi)
        ct = cos(theta)
        st = sin(theta)
        cs = cos(psi)
        ss = sin(psi)

        rot_matrix(1, 1) = cp * ct * cs - sp * ss
        rot_matrix(1, 2) = -cp * ct * ss - sp * cs
        rot_matrix(1, 3) = cp * st
        rot_matrix(2, 1) = sp * ct * cs + cp * ss
        rot_matrix(2, 2) = -sp * ct * ss + cp * cs
        rot_matrix(2, 3) = sp * st
        rot_matrix(3, 1) = -st * cs
        rot_matrix(3, 2) = st * ss
        rot_matrix(3, 3) = ct

        do i_atom = 1, n_atoms
            pos_cluster(:, i_atom) = matmul(rot_matrix, pos_cluster(:, i_atom))
        enddo
        id_sim_bis = (id_sim-1) * 10 + augmentation_index - 1
    endsubroutine


    subroutine place_in_box
    
        use constants, only: lateral_border_margin, box_hrtem
        use variable, only: placed, shift, pos_cluster
        use descriptor, only: n_atoms

        integer             ::  i_axis
        double precision    ::  min_coord(3), max_coord(3), span(3), free_lateral(2), u(2), box_size(3)

        placed = .false.
        shift = 0.0d0
        box_size = box_hrtem * 10.0d0
        do i_axis = 1, 3
            min_coord(i_axis) = minval(pos_cluster(i_axis, 1:n_atoms))
            max_coord(i_axis) = maxval(pos_cluster(i_axis, 1:n_atoms))
        enddo
        span = max_coord - min_coord
        if (any(span .ge. box_size - 1.0d-8)) return
        free_lateral = box_size(1:2) - span(1:2) - 2.0d0 * lateral_border_margin
        if (any(free_lateral .lt. 0.0d0)) return

        call random_number(u)
        shift(1:2) = -min_coord(1:2) + lateral_border_margin + u * free_lateral
        shift(3) = -min_coord(3) + 0.5d0 * (box_size(3) - span(3))
        pos_cluster(:, 1:n_atoms) = pos_cluster(:, 1:n_atoms) + spread(shift, dim=2, ncopies=n_atoms)
        pos_cluster(3, 1:n_atoms) = pos_cluster(3, 1:n_atoms) - sum(pos_cluster(3, 1:n_atoms)) / dble(n_atoms)
        placed = .true.

    endsubroutine


    subroutine choose_substrate

        use variable, only: substrate_file, substrate_files, n_substrate_files, substrate_dir

        integer :: k
        double precision :: u

        call random_number(u)
        k = min(int(u * n_substrate_files) + 1, n_substrate_files)
        substrate_file = trim(substrate_dir) // '/' // trim(substrate_files(k))

    endsubroutine


    subroutine read_substrate_file

        use constants, only: n_atoms_substrate_max, box_hrtem
        use variable, only: substrate_file, n_atoms_substrate, pos_substrate, species_substrate, dbf_c

        integer ::  i_atom, index0, index1
        double precision    ::  lattice(9)
        character(len=255)  ::  line

        open(10, file=substrate_file, status='old', action='read')
        read(10, *) n_atoms_substrate
        if (n_atoms_substrate .gt. n_atoms_substrate_max) then
            write(*,*) 'error: n_atoms_substrate = ', n_atoms_substrate, ' exceeds n_atoms_substrate_max = ', &
                n_atoms_substrate_max, ' in ', trim(substrate_file)
            stop
        endif
        read(10, '(A)') line
        do i_atom = 1, n_atoms_substrate
            read(10, *) species_substrate(i_atom)(1:2), pos_substrate(1, i_atom), pos_substrate(2, i_atom), &
                            pos_substrate(3, i_atom), dbf_c(i_atom)
        enddo
        close(10)
        ! Le decalage lateral aleatoire (add_substrate_to_cluster) suppose un substrat periodique de la taille de la boite
        index0 = index(line, 'Lattice="') + len('Lattice="')
        index1 = index(line(index0:), '"') + index0 - 2
        read(line(index0:index1), *) lattice
        if (abs(lattice(1) - box_hrtem(1)*10.0d0) .gt. 1.0d-3 .or. abs(lattice(5) - box_hrtem(2)*10.0d0) .gt. 1.0d-3) then
            write(*,*) 'error: substrate cell ', lattice(1), lattice(5), ' does not match box_hrtem in ', trim(substrate_file)
            stop
        endif
        dbf_c(1:n_atoms_substrate) = 0.01d0 *sum(dbf_c(1:n_atoms_substrate)) / dble(n_atoms_substrate)

    endsubroutine


    subroutine add_substrate_to_cluster
        
        use constants, only: box_hrtem
        use variable, only: n_atoms_tot, pos_cluster, pos_substrate, species_substrate, species, n_atoms_substrate, placed
        use descriptor, only: n_atoms, id_sim_bis

        integer ::  i_atom
        double precision :: z_shift, u(2), xy_shift(2), box_xy(2)

        ! Decalage lateral aleatoire du substrat: exact car le substrat est periodique en x, y avec la periode de la boite
        box_xy = box_hrtem(1:2) * 10.0d0
        call random_number(u)
        xy_shift = u * box_xy
        z_shift = minval(pos_cluster(3, 1:n_atoms)) - maxval(pos_substrate(3, 1:n_atoms_substrate)) - 2.0d0

        do i_atom = 1, n_atoms_substrate
            pos_cluster(1:2, n_atoms + i_atom) = modulo(pos_substrate(1:2, i_atom) + xy_shift, box_xy)
            pos_cluster(3, n_atoms + i_atom) = pos_substrate(3, i_atom) + z_shift
            species(n_atoms + i_atom) = species_substrate(i_atom)
        enddo
        n_atoms_tot = n_atoms + n_atoms_substrate

        if (maxval(pos_cluster(3, 1:n_atoms_tot)) .ge. 0.5d0 * box_hrtem(3) * 10.0d0 .or. minval(pos_cluster(3, 1:n_atoms_tot)) &
                    &.lt. -0.5d0 * box_hrtem(3) * 10.0d0) then
            placed = .false.
            call save_xyz
            call save_xyz2
            write(*,*) "Cluster with substrate does not fit in the box height. Skipping. (id_sim = ", id_sim_bis, ")"
        endif

    endsubroutine
    

    subroutine save_xyz

        use constants, only: box_hrtem
        use variable, only: pos_cluster, species, n_atoms_tot
        use descriptor, only: id_sim_bis

        integer             ::  i_atom, unit_xyz
        double precision    ::  lattice(3)
        character(len=255)  ::  filename, header
        character(len=255)  ::  id_str

        lattice = box_hrtem * 10.0d0
        write(id_str, '(I0)') id_sim_bis
        filename = 'xyz_substrate/' // trim(id_str) // '.xyz'

        write(header, '(A,F0.3,A,F0.3,A,F0.3,A,F0.3,A)') 'Lattice="', lattice(1), ' 0.0 0.0 0.0 ', lattice(2), &
            ' 0.0 0.0 0.0 ', lattice(3), '" Origin="0.0 0.0 ', -0.5d0*lattice(3), '" Properties=species:S:1:pos:R:3 pbc="T T F"'
        open(newunit=unit_xyz, file=trim(filename), status='replace', action='write')
        write(unit_xyz, '(I0)') n_atoms_tot
        write(unit_xyz, '(A)') trim(header)
        do i_atom = 1, n_atoms_tot
            write(unit_xyz, '(A,3(1X,F12.5))') trim(species(i_atom)), pos_cluster(:, i_atom)
        enddo
        close(unit_xyz)

    endsubroutine


    subroutine save_xyz2

        use variable, only: pos, species, box
        use descriptor, only: id_sim_bis, n_atoms

        integer             ::  i_atom, unit_xyz
        character(len=255)  ::  filename, header
        character(len=255)  ::  id_str

        write(id_str, '(I0)') id_sim_bis
        filename = 'xyz_substrate/' // trim(id_str) // '_original.xyz'

        write(header, '(A,F0.3,A,F0.3,A,F0.3,A)') 'Lattice="', box(1), ' 0.0 0.0 0.0 ', box(2), &
            ' 0.0 0.0 0.0 ', box(3), '" Properties=species:S:1:pos:R:3 pbc="T T F"'
        open(newunit=unit_xyz, file=trim(filename), status='replace', action='write')
        write(unit_xyz, '(I0)') n_atoms
        write(unit_xyz, '(A)') trim(header)
        do i_atom = 1, n_atoms
            write(unit_xyz, '(A,3(1X,F12.5))') trim(species(i_atom)), pos(:, i_atom)
        enddo
        close(unit_xyz)

    endsubroutine

    
    subroutine prepare_hrtem_particle
        
        use constants, only: dbf_ag, dbf_co, box_hrtem
        use descriptor, only: n_atoms
        use variable, only: species, atom_typ1, atom_typ2, pos_cluster, atomic_number, biso, dbf_c, n_atoms_tot

        integer ::  i_atom

        do i_atom = 1, n_atoms_tot
            if (species(i_atom) .eq. atom_typ1) then
                atomic_number(i_atom) = 47
                biso(i_atom) = dbf_ag
            else if (species(i_atom) .eq. atom_typ2) then
                atomic_number(i_atom) = 27
                biso(i_atom) = dbf_co
            else 
                atomic_number(i_atom) = 6
                biso(i_atom) = dbf_c(i_atom - n_atoms)
            endif
            pos_cluster(:, i_atom) = pos_cluster(:, i_atom) / (box_hrtem*10.0d0)
        enddo

    endsubroutine
    
end module