module slc

    use constants, only: nx_max, ny_max, nz_max, n_types_max, n_atoms_tot_max

    implicit none

    integer             ::  n_types = 0, type_index(n_atoms_tot_max), atomic_number_type(n_types_max)
    integer             ::  slice_count(nz_max)
    double precision    ::  volume_slc, f_re, f_im, gthr2, gx(nx_max), gy(ny_max), biso_type(n_types_max), a(2), b(6)
    double complex      ::  uhat(nx_max, ny_max), trans_fft(nx_max, ny_max), ftab(nx_max, ny_max, n_types_max)

    integer             ::  slice(nz_max, n_atoms_tot_max), m_fft(nx_max)

    ! Facteur de structure par blocs de n_blk atomes (slice_potential): phases selon x et y, produit par ZGEMM
    integer, parameter  ::  n_blk = 256
    double complex      ::  phase_x_blk(nx_max, n_blk), phase_y_blk(ny_max, n_blk), sfac(nx_max, ny_max)

    interface
        subroutine zgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
            character, intent(in)           ::  transa, transb
            integer, intent(in)             ::  m, n, k, lda, ldb, ldc
            double complex, intent(in)      ::  alpha, beta, a(lda, *), b(ldb, *)
            double complex, intent(inout)   ::  c(ldc, *)
        endsubroutine
    endinterface


    double precision, parameter ::  v(1:98) = (/ 0.5d0, 0.5d0, 0.5d0, 0.3d0, 0.5d0, 0.5d0, 0.5d0, 0.5d0, 0.5d0, 0.5d0, &
                                                & 0.5d0, 0.5d0, 0.4d0, 0.5d0, 0.5d0, 0.5d0, 0.5d0, 0.5d0, 0.2d0, 0.3d0, &
                                                & 0.5d0, 0.5d0, 0.5d0, 0.5d0, 0.5d0, 0.5d0, 0.5d0, 0.5d0, 0.5d0, 0.5d0, &
                                                & 0.5d0, 0.5d0, 0.5d0, 0.5d0, 0.5d0, 0.5d0, 0.2d0, 0.3d0, 0.5d0, 0.5d0, &
                                                & 0.5d0, 0.5d0, 0.5d0, 0.4d0, 0.5d0, 0.5d0, 0.5d0, 0.3d0, 0.4d0, 0.6d0, &
                                                & 0.6d0, 0.6d0, 0.4d0, 0.4d0, 0.1d0, 0.1d0, 0.3d0, 0.3d0, 0.2d0, 0.2d0, &
                                                & 0.2d0, 0.2d0, 0.1d0, 0.2d0, 0.1d0, 0.2d0, 0.1d0, 0.2d0, 0.1d0, 0.1d0, &
                                                & 0.1d0, 0.1d0, 0.4d0, 0.2d0, 0.5d0, 0.4d0, 0.5d0, 0.5d0, 0.4d0, 0.4d0, &
                                                & 0.4d0, 0.3d0, 0.4d0, 0.4d0, 0.4d0, 0.4d0, 0.1d0, 0.2d0, 0.2d0, 0.3d0, &
                                                & 0.2d0, 0.2d0, 0.2d0, 0.2d0, 0.2d0, 0.3d0, 0.2d0, 0.3d0 /)
    double precision, parameter ::  bb(6, 1:98) = reshape((/ &
            &  48.75740d0,  4.96588d0, 18.24440d0, 18.24440d0, 18.24440d0, 18.24440d0, &
            &  2.54216d0,  8.74302d0, 12.69098d0,  0.43711d0,  5.29446d0, 28.25045d0, &
            &  0.68454d0,  3.06497d0,  6.23974d0,126.17816d0,131.20160d0,131.76538d0, &
            &  0.53996d0,  3.38752d0, 55.62340d0, 50.78098d0, 67.00502d0, 96.36635d0, &
            &  0.33138d0,  2.97485d0, 34.01118d0, 35.98365d0, 36.68364d0, 60.80991d0, &
            &  0.29458d0,  3.93381d0, 24.97836d0, 25.27916d0, 25.46696d0, 46.70328d0, &
            &  0.23925d0,  4.93515d0, 18.11895d0, 15.69698d0, 15.81922d0, 40.24150d0, &
            &  6.37582d0,  8.03744d0, 27.20649d0,  0.11157d0,  0.38686d0, 10.89944d0, &
            &  0.21800d0,  6.76987d0,  7.05056d0,  6.67484d0, 12.38148d0, 28.08398d0, &
            &  0.20055d0,  5.49814d0,  6.28052d0,  7.19211d0,  7.54763d0, 23.26388d0, &
            &  0.21902d0,  5.30022d0,  5.31938d0,  5.28281d0,  5.28546d0,128.18391d0, &
            &  1.97633d0,  2.80902d0, 16.39184d0,  0.05494d0,  2.06121d0,121.70512d0, &
            &  2.29692d0,  2.35822d0, 24.98576d0,  0.07462d0,  0.55953d0,128.50104d0, &
            &  1.73656d0,  3.04329d0, 30.57191d0,  0.05070d0,  0.99181d0, 86.18340d0, &
            &  0.17949d0,  2.63250d0,  2.67559d0, 34.57098d0, 36.77888d0, 54.06180d0, &
            &  1.00609d0,  4.90414d0, 31.34909d0,  0.03699d0,  0.98700d0, 44.94354d0, &
            &  0.18464d0,  1.47963d0,  5.20989d0, 24.79470d0, 32.06184d0, 39.09933d0, &
            &  0.20060d0,  6.53262d0, 22.72092d0,  1.20022d0,  1.27398d0, 36.25907d0, &
            &  0.44424d0,  3.36735d0, 19.63031d0,  0.01824d0, 23.51332d0,212.86819d0, &
            &  0.18274d0,  2.06638d0, 16.99062d0, 11.57795d0, 13.97594d0,186.10446d0, &
            &  0.14245d0,  1.46588d0, 15.46955d0,  4.24287d0,  9.80399d0,121.46864d0, &
            &  0.12782d0,  1.45591d0, 12.09738d0,  4.61747d0, 11.96791d0,105.00546d0, &
            &  0.13126d0,  1.39923d0,  8.00762d0,  7.98129d0, 13.41408d0, 95.30811d0, &
            &  0.12311d0,  2.38386d0,  9.92149d0,  1.64793d0, 11.00035d0, 68.45583d0, &
            &  0.48173d0,  3.78306d0,  8.47337d0,  0.04690d0,  8.74544d0, 77.44405d0, &
            &  0.44704d0,  6.89364d0,  6.90335d0,  0.05691d0,  3.02647d0, 70.86599d0, &
            &  0.10705d0,  3.63573d0,  7.55825d0,  1.27986d0,  5.14045d0, 67.16051d0, &
            &  0.11069d0,  1.61889d0,  6.00325d0,  5.97496d0,  6.06049d0, 59.41419d0, &
            &  0.11293d0,  1.89077d0,  5.08503d0,  5.07335d0,  5.09928d0, 46.38955d0, &
            &  0.10209d0,  1.73365d0,  4.78298d0,  4.80706d0,  5.64485d0, 51.21828d0, &
            &  0.10642d0,  1.53735d0,  5.13798d0,  4.74298d0,  4.99974d0, 61.42872d0, &
            &  0.09583d0,  1.67715d0,  4.70275d0,  2.91198d0,  7.87009d0, 64.93623d0, &
            &  0.09428d0,  2.21409d0,  3.95060d0,  1.52064d0, 15.81446d0, 52.41380d0, &
            &  0.09252d0,  1.60168d0,  3.04917d0,  3.18476d0, 18.93890d0, 47.62742d0, &
            &  0.09246d0,  1.77298d0,  3.48134d0,  1.88354d0, 22.68630d0, 40.69434d0, &
            &  0.49321d0,  2.08254d0, 11.41282d0,  0.03333d0,  2.09673d0, 42.38068d0, &
            &  0.15796d0,  1.71505d0,  9.39164d0,  1.67464d0, 23.58663d0,152.53635d0, &
            &  0.36052d0,  2.12757d0, 12.45815d0,  0.01526d0,  2.10824d0,133.17088d0, &
            &  0.09003d0,  1.41396d0,  2.05348d0, 10.25766d0, 10.74831d0, 90.63555d0, &
            &  0.10094d0,  1.15419d0,  2.34669d0, 10.58145d0, 10.94962d0, 82.82259d0, &
            &  0.09243d0,  1.16977d0,  5.93969d0,  1.30554d0, 13.43475d0, 66.37486d0, &
            &  0.43543d0,  1.24830d0,  7.45369d0,  0.03543d0,  9.91366d0, 61.72203d0, &
            &  0.45943d0,  1.18155d0,  8.31728d0,  0.03226d0,  8.32296d0, 64.97874d0, &
            &  0.08603d0,  1.39552d0, 11.69728d0,  1.39552d0,  3.45200d0, 55.55519d0, &
            &  0.09214d0,  1.11341d0,  7.65767d0,  1.12566d0,  8.32517d0, 48.38017d0, &
            &  0.09005d0,  1.12460d0,  9.69801d0,  1.08539d0,  5.70912d0, 33.48585d0, &
            &  0.08938d0,  3.19060d0,  9.10000d0,  0.80898d0,  0.81439d0, 41.34453d0, &
            &  0.28851d0,  1.61312d0,  8.99691d0,  0.01711d0,  9.46666d0, 58.13256d0, &
            &  0.08948d0,  1.23258d0,  8.23129d0,  1.22390d0,  7.06201d0, 59.69622d0, &
            &  0.07124d0,  0.85532d0,  6.40081d0,  1.33637d0,  6.38240d0, 50.92361d0, &
            &  0.35749d0,  1.32481d0,  6.51696d0,  0.03550d0,  6.51913d0, 50.80984d0, &
            &  0.50089d0,  3.95301d0,  7.62830d0,  0.03005d0,  0.50737d0, 49.62628d0, &
            &  0.08429d0,  1.12959d0,  8.86209d0,  1.12981d0,  9.13243d0, 56.01965d0, &
            &  0.27796d0,  1.62147d0, 11.45200d0,  0.02032d0,  3.27497d0, 51.44078d0, &
            &  0.12045d0,  1.53654d0,  9.81569d0, 41.21656d0, 42.62216d0,224.34816d0, &
            &  0.12230d0,  1.44909d0,  9.50159d0, 49.40860d0, 74.94942d0,217.04485d0, &
            &  0.08930d0,  1.26225d0,  8.09703d0,  1.20293d0, 17.65554d0,116.61481d0, &
            &  0.08504d0,  1.28286d0, 11.22123d0,  1.32741d0,  4.61040d0,112.19678d0, &
            &  0.09805d0,  1.52628d0,  8.58953d0,  1.23893d0, 22.49126d0,140.02856d0, &
            &  0.09413d0,  1.26616d0,  5.98844d0, 17.78775d0, 18.14397d0,132.59305d0, &
            &  0.09447d0,  1.25111d0,  5.91205d0, 16.28675d0, 16.73089d0,127.90916d0, &
            &  0.09061d0,  1.59281d0, 10.64077d0,  1.78861d0,  2.22148d0,124.56328d0, &
            &  0.10485d0,  1.54396d0,  8.65223d0,  7.09290d0, 53.36537d0,183.69014d0, &
            &  0.09338d0,  1.38681d0,  7.35883d0,  1.55122d0, 20.81916d0,111.03201d0, &
            &  0.10190d0,  1.52368d0,  7.16923d0, 20.86269d0, 49.29465d0,166.09206d0, &
            &  0.08402d0,  1.40890d0,  7.14042d0,  1.34848d0, 11.42203d0,108.01204d0, &
            &  0.09441d0,  1.61807d0,  6.27142d0, 40.34946d0, 42.82722d0,130.59616d0, &
            &  0.08211d0,  1.25106d0,  4.81241d0, 10.84493d0, 10.90164d0,100.07855d0, &
            &  0.09662d0,  1.60236d0,  5.67480d0, 30.59014d0, 31.12732d0,138.69682d0, &
            &  0.09493d0,  1.60220d0,  5.43916d0, 28.31076d0, 29.27660d0,138.08665d0, &
            &  0.09658d0,  1.56751d0,  5.32170d0, 34.18217d0, 35.25187d0,121.42893d0, &
            &  0.09294d0,  1.55499d0,  5.25121d0, 37.51883d0, 38.88302d0,105.16978d0, &
            &  0.06298d0,  0.81950d0,  2.89124d0,  5.54290d0,  5.98101d0, 54.42459d0, &
            &  0.07902d0,  1.37096d0,  8.23364d0,  1.38300d0,  1.39219d0, 77.11813d0, &
            &  0.05266d0,  0.90718d0,  4.43830d0,  0.94590d0,  4.37477d0, 43.97909d0, &
            &  0.22700d0,  1.56975d0,  6.34451d0,  0.01564d0,  1.61769d0, 46.15815d0, &
            &  0.05055d0,  0.86775d0,  5.09325d0,  0.88123d0,  3.56919d0, 39.77390d0, &
            &  0.05253d0,  0.83773d0,  3.95899d0,  0.81515d0,  6.44217d0, 34.21146d0, &
            &  0.54927d0,  1.72752d0,  6.71952d0,  0.02637d0,  0.07253d0, 35.45745d0, &
            &  0.21941d0,  1.41611d0,  6.68241d0,  0.01472d0,  1.57578d0, 37.15826d0, &
            &  0.22459d0,  1.12822d0,  4.30289d0,  0.01485d0,  7.15607d0, 43.08737d0, &
            &  0.06432d0,  1.19406d0,  7.39342d0,  1.14160d0,  1.28905d0, 51.13401d0, &
            &  0.05380d0,  0.86719d0,  1.87540d0,  7.64796d0,  7.86794d0, 45.63897d0, &
            &  0.50112d0,  1.63784d0,  6.78551d0,  0.02187d0,  0.08602d0, 46.72951d0, &
            &  0.22321d0,  1.10827d0,  3.59116d0,  0.01011d0, 11.63732d0, 45.06839d0, &
            &  0.21152d0,  1.14015d0,  3.41473d0,  0.01188d0, 13.41211d0, 43.11389d0, &
            &  0.09435d0,  1.02649d0,  6.25480d0, 32.51444d0, 36.29119d0,149.11722d0, &
            &  0.07300d0,  1.01825d0,  5.89629d0,  1.03089d0, 20.37389d0,115.34722d0, &
            &  0.07515d0,  0.94941d0,  3.72527d0, 17.58346d0, 19.75388d0,109.12856d0, &
            &  0.06385d0,  0.90194d0,  4.65715d0,  0.90253d0, 15.70771d0, 83.69695d0, &
            &  0.07557d0,  0.84920d0,  4.00991d0, 16.95003d0, 17.78767d0,100.20415d0, &
            &  0.07142d0,  1.14907d0,  9.21231d0,  0.95923d0,  1.20275d0,104.32746d0, &
            &  0.06918d0,  0.98102d0,  5.95437d0,  0.99086d0, 22.06437d0, 90.98156d0, &
            &  0.07136d0,  0.95772d0,  6.13183d0,  0.97438d0, 15.67499d0, 89.86625d0, &
            &  0.07301d0,  0.93267d0,  6.34836d0,  0.91032d0, 13.26179d0, 86.85986d0, &
            &  0.05778d0,  0.72273d0,  3.01146d0,  9.21882d0,  9.53410d0, 65.86810d0, &
            &  0.07088d0,  0.77587d0,  6.14295d0,  1.79036d0, 15.12379d0, 83.56983d0, &
            &  0.06164d0,  0.81363d0,  6.56165d0,  0.83805d0,  4.18914d0, 61.41408d0 /), (/ 6, 98 /))

    
    contains



    subroutine run_slc

        use constants, only: sigma0, box_hrtem
        use variable, only: nx, ny, nz, lambda, dx, dy, dz, trans, gmax, n_atoms_tot, pos_cluster
        use fft, only: fft_index, fft2

        integer             ::  i_px, mx, my, index_slice, i_atom
        double precision    ::  sigma_lambda
        double complex      ::  pot(nx, ny), trans_filtered(nx, ny)


        sigma_lambda = sigma0 * lambda
        dx = box_hrtem(1) / dble(nx)
        dy = box_hrtem(2) / dble(ny)
        dz = box_hrtem(3) / dble(nz)
        volume_slc = box_hrtem(1) * box_hrtem(2) * dz

        do i_px = 1, nx
            mx = fft_index(i_px, nx) ! indices pour espaces reciproques et grille centrée en 0
            m_fft(i_px) = mx
            gx(i_px) = dble(mx) / box_hrtem(1) ! Frequence spatiale en x
            my = fft_index(i_px, ny)
            gy(i_px) = dble(my) / box_hrtem(2) ! Frequence spatiale en y
        enddo
        gmax = 0.5d0*nx / box_hrtem(1) ! Limite de Nyquist

        call unique_scattering_factors ! Type (Z, B) de chaque atome, facteurs de diffusion calcules une seule fois par type

        slice_count = 0
        do i_atom = 1, n_atoms_tot
            index_slice = int((pos_cluster(3, i_atom) + 0.5d0) * dble(nz)) + 1
            if (index_slice .lt. 1) index_slice = 1
            if (index_slice .gt. nz) index_slice = nz
            slice_count(index_slice) = slice_count(index_slice) + 1
            slice(index_slice, slice_count(index_slice)) = i_atom
        enddo

        do index_slice = 1, nz
            call slice_potential(index_slice)
            call fft2(uhat(1:nx, 1:ny), pot, 1)
            trans(1:nx, 1:ny, index_slice) = exp(dcmplx(0.0d0, sigma_lambda * dz) * pot)
            call fft2(trans(1:nx, 1:ny, index_slice), trans_fft(1:nx, 1:ny), -1)
            gthr2 = (gmax * (2.0d0/3.0d0)) ** 2
            call apply_hard_aperture
            call fft2(trans_fft(1:nx, 1:ny), trans_filtered, 1)
            trans(1:nx, 1:ny, index_slice) = trans_filtered / (nx * ny)
        enddo
        if (sum(slice_count(1:nz)) .ne. n_atoms_tot) &
            write(*,*) 'WARNING: ', n_atoms_tot - sum(slice_count(1:nz)), ' atoms outside the slices'


    endsubroutine



    subroutine apply_hard_aperture

        use constants, only: box_hrtem
        use variable, only: nx, ny
        use fft, only: fft_index

        integer             ::  i_px, j_px, mx, my
        double precision    ::  gx_hrtem, gy_hrtem

        do j_px = 1, ny        
            my = fft_index(j_px, ny)
            gy_hrtem = my / box_hrtem(2)
            do i_px = 1, nx
                mx = fft_index(i_px, nx)
                gx_hrtem = mx / box_hrtem(1)
                if (gx_hrtem**2 + gy_hrtem**2 .gt. gthr2) trans_fft(i_px, j_px) = dcmplx(0.0d0, 0.0d0)
            enddo
        enddo
    
    endsubroutine



    subroutine slice_potential(index_slice)
    
        use constants, only: v0
        use variable, only: pos_cluster, nx, ny

        integer             ::  k, i_atom, i_type, n_in_blk
        logical             ::  first_blk
        integer, intent(in) ::  index_slice

        ! uhat(g) = v0/V * somme sur les types t de ftab_t(g) * S_t(g), ou le facteur de structure du type t dans la
        ! tranche S_t(gx, gy) = somme sur les atomes de exp(-2 pi i gx x) exp(-2 pi i gy y) = (Px . Py^T)(gx, gy)
        ! est un produit matriciel (ZGEMM), accumule par blocs de n_blk atomes
        uhat = dcmplx(0.0d0, 0.0d0)
        do i_type = 1, n_types
            first_blk = .true.
            n_in_blk = 0
            do k = 1, slice_count(index_slice)
                i_atom = slice(index_slice, k)
                if (type_index(i_atom) .ne. i_type) cycle
                n_in_blk = n_in_blk + 1
                call atom_phases(pos_cluster(1, i_atom), nx, phase_x_blk(:, n_in_blk))
                call atom_phases(pos_cluster(2, i_atom), ny, phase_y_blk(:, n_in_blk))
                if (n_in_blk .eq. n_blk) call add_block(n_in_blk, first_blk)
            enddo
            if (n_in_blk .gt. 0) call add_block(n_in_blk, first_blk)
            if (.not. first_blk) uhat(1:nx, 1:ny) = uhat(1:nx, 1:ny) + ftab(1:nx, 1:ny, i_type) * sfac(1:nx, 1:ny)
        enddo
        uhat = (v0 / volume_slc) * uhat

    endsubroutine



    subroutine add_block(n_in_blk, first_blk)

        use variable, only: nx, ny

        integer, intent(inout)  ::  n_in_blk
        logical, intent(inout)  ::  first_blk
        double complex          ::  beta

        ! sfac = Px . Py^T pour le premier bloc du type, sfac = sfac + Px . Py^T ensuite
        beta = dcmplx(1.0d0, 0.0d0)
        if (first_blk) beta = dcmplx(0.0d0, 0.0d0)
        call zgemm('N', 'T', nx, ny, n_in_blk, dcmplx(1.0d0, 0.0d0), phase_x_blk, nx_max, phase_y_blk, ny_max, &
                    beta, sfac, nx_max)
        first_blk = .false.
        n_in_blk = 0

    endsubroutine



    subroutine atom_phases(frac, n, ph)

        use constants, only: pi

        integer                         ::  i, m
        double complex                  ::  w, pw(0:nx_max)
        integer, intent(in)             ::  n
        double precision, intent(in)    ::  frac
        double complex, intent(out)     ::  ph(n)

        ! ph(i) = exp(-2 pi i m frac) pour la frequence m = m_fft(i), par puissances successives de exp(-2 pi i frac)
        w = dcmplx(cos(2.0d0*pi*frac), -sin(2.0d0*pi*frac))
        pw(0) = dcmplx(1.0d0, 0.0d0)
        do m = 1, n/2
            pw(m) = pw(m-1) * w
        enddo
        do i = 1, n
            m = m_fft(i)
            if (m .ge. 0) then
                ph(i) = pw(m)
            else
                ph(i) = conjg(pw(-m))
            endif
        enddo

    endsubroutine



    subroutine unique_scattering_factors

        use variable, only: atomic_number, biso, n_atoms_tot

        integer ::  i_atom, i_type

        ! La table des types et ftab sont conserves d'une image a l'autre (grille, tension et B fixes pendant tout le
        ! calcul): ftab n'est calcule qu'a la premiere apparition d'un couple (Z, B)
        do i_atom = 1, n_atoms_tot
            type_index(i_atom) = 0
            do i_type = 1, n_types
                if (atomic_number(i_atom) .eq. atomic_number_type(i_type) .and. abs(biso(i_atom)-biso_type(i_type)).le.1.0d-12) then
                    type_index(i_atom) = i_type
                    exit
                endif
            enddo
            if (type_index(i_atom) .eq. 0) then
                if (n_types .ge. n_types_max) stop 'error: more than n_types_max distinct (Z, B) pairs, see README.md'
                n_types = n_types + 1
                atomic_number_type(n_types) = atomic_number(i_atom)
                biso_type(n_types) = biso(i_atom)
                type_index(i_atom) = n_types
                call compute_form_factor(n_types)
            endif
        enddo

    endsubroutine



    subroutine compute_form_factor(index_type)

        use variable, only: nx, ny, g2, gmax

        integer             ::  i_px, j_px
        double precision    ::  apod
        integer, intent(in) ::  index_type

        call get_weko(index_type)
        do j_px = 1, ny
            do i_px = 1, nx
                g2 = gx(i_px)**2 + gy(j_px)**2
                apod = 0.5d0 - 0.5d0*tanh((sqrt(g2)/gmax - 0.9d0) * 30.0d0)
                call scattering_factor(index_type)
                ftab(i_px, j_px, index_type) = apod * dcmplx(f_re, f_im)
            enddo
        enddo

    endsubroutine
    


    subroutine scattering_factor(index_type)

        use constants, only: e0, hc, pi, r8pi2
        use variable, only: g2, ht

        double precision    ::  gamma, ua, ga, k0, dwf, g
        integer, intent(in) ::  index_type

        g = sqrt(g2)
        gamma = (e0 + ht) / e0
        dwf = exp(-0.25d0*biso_type(index_type)*g2)
        f_re = 0.1d0 * gamma * 4.0d0*pi * weko_real(0.05d0 * g) * dwf
        f_im = 0.0d0
        if(biso_type(index_type) .gt. 0.0d0 .and. abs(f_re) .gt. tiny(f_re)) then
            ua = 10.0d0 * sqrt(biso_type(index_type) * r8pi2)
            ga = 0.1d0 * g * 2.0d0*pi
            k0 = 2.0d0*pi/hc*0.1d0 * sqrt((2.0d0*e0 + ht) * ht)
            f_im = 0.1d0 * gamma**2 * weko_imag(ga, ua) / k0
        endif

    endsubroutine



    subroutine get_weko(index_type)

        integer             ::  i
        double precision    ::  v_val
        integer, intent(in) ::  index_type
        
        v_val = v(atomic_number_type(index_type))
        a(1) = 0.023933659d0 * dble(atomic_number_type(index_type)) / (3.0d0 * (1.0d0 + v_val))
        a(2) = v_val * a(1)
        do i = 1, 6
            b(i) = bb(i, atomic_number_type(index_type))
        enddo

    endsubroutine


    
    double precision function weko_real(s)

        integer                         ::  i, j
        double precision                ::  argu
        double precision, intent(in)    ::  s

        weko_real = 0.0d0
        do i = 1, 6
            j = 1 + (i-1)/3
            argu = b(i)*s**2
            if(argu.lt.0.1d0) then
                weko_real = weko_real + a(j)*b(i)*(1.0d0 - 0.5d0*argu)
            else if(argu.gt.20.0d0) then
                weko_real = weko_real + a(j)/s**2
            else
                weko_real = weko_real + a(j)*(1.0d0 - exp(-argu))/s**2
            endif
        enddo
    endfunction



    double precision function weko_imag(g, ul)

        use constants, only: pi

        double precision, intent(in)    ::  g, ul
        integer                         ::  i, j, ii, jj
        double precision                ::  b1(6), dwf

        do i = 1, 6
            b1(i) = b(i) / (4.0d0*pi*pi)
        enddo

        dwf = exp(-0.5d0*ul**2*g**2)
        weko_imag = 0.0d0
        do j = 1, 6
            jj = 1 + (j-1)/3
            do i = 1, 6
                ii = 1 + (i-1)/3
                weko_imag = weko_imag + a(jj)*a(ii)*(dwf*ri1(b1(i), b1(j), g) - ri2(b1(i), b1(j), g, ul))
            enddo
        enddo
        weko_imag = 4.0d0*pi**4 * weko_imag
    
    endfunction



    double precision function ri1(bi, bj, g)

        use constants, only: pi, euler

        double precision                ::  big2, bjg2, x1, x2, g2
        double precision, intent(in)    ::  bi, bj, g

        if(abs(g).le.tiny(g)) then
            ri1 = pi * (bi*log((bi+bj)/bi) + bj*log((bi+bj)/bj))
            return
        endif
        g2 = g**2
        big2 = bi*g2
        bjg2 = bj*g2
        ri1 = 2.0d0*euler + log(big2) + log(bjg2) - 2.0d0*ei(-bi*bj*g2/(bi+bj))
        x1 = big2
        x2 = big2*bi/(bi+bj)
        ri1 = ri1 + rih1(x1, x2, x1)
        x1 = bjg2
        x2 = bjg2*bj/(bi+bj)
        ri1 = ri1 + rih1(x1, x2, x1)
        ri1 = pi * ri1 / g2

    endfunction



    double precision function ri2(bi, bj, g, ul)

        use constants, only: pi

        double precision                ::  biu, bju, biuh, bjuh, g2, u2, x1, x2, x3
        double precision, intent(in)    ::  bi, bj, g, ul

        u2 = ul**2
        g2 = g**2
        if(abs(g).le.tiny(g)) then
            ri2 = (bi+u2) * log((bi+bj+u2)/(bi+u2))
            ri2 = ri2 + bj * log((bi+bj+u2)/(bj+u2))
            ri2 = ri2 + u2 * log(u2/(bj+u2))
            ri2 = pi * ri2
            return
        endif
        biuh = bi + 0.5d0*u2
        bjuh = bj + 0.5d0*u2
        biu = bi + u2
        bju = bj + u2
        ri2 = ei(-0.5d0*u2*g2*biuh/biu) + ei(-0.5d0*u2*g2*bjuh/bju)
        ri2 = ri2 - ei(-biuh*bjuh*g2/(biuh+bjuh)) - ei(-0.25d0*u2*g2)
        ri2 = 2.0d0 * ri2
        x1 = 0.5d0*u2*g2
        x2 = 0.25d0*u2*g2
        x3 = 0.25d0*u2*u2*g2/biu
        ri2 = ri2 + rih1(x1, x2, x3)
        x3 = 0.25d0*u2*u2*g2/bju
        ri2 = ri2 + rih1(x1, x2, x3)
        x1 = biuh*g2
        x2 = biuh*biuh*g2/(biuh+bjuh)
        x3 = biuh*biuh*g2/biu
        ri2 = ri2 + rih1(x1, x2, x3)
        x1 = bjuh*g2
        x2 = bjuh*bjuh*g2/(biuh+bjuh)
        x3 = bjuh*bjuh*g2/bju
        ri2 = ri2 + rih1(x1, x2, x3)
        ri2 = pi * ri2 / g2
    
    endfunction



    double precision function rih1(x1, x2, x3)

        double precision, intent(in)    ::  x1, x2, x3

        if(x2.le.20.0d0 .and. x3.le.20.0d0) then
            rih1 = exp(-x1) * (ei(x2) - ei(x3))
            return
        endif
        if(x2.gt.20.0d0) then
            rih1 = exp(x2-x1) * rih2(x2) / x2
        else
            rih1 = exp(-x1) * ei(x2)
        endif
        if(x3.gt.20.0d0) then
            rih1 = rih1 - exp(x3-x1) * rih2(x3) / x3
        else
            rih1 = rih1 - exp(-x1) * ei(x3)
        endif
    endfunction



    double precision function rih2(x)

        integer                         ::  i, i1
        double precision                ::  x1
        double precision, intent(in)    ::  x
        double precision, parameter     ::  f(0:20) = (/ &
            1.000000d0, 1.005051d0, 1.010206d0, 1.015472d0, 1.020852d0, &
            1.026355d0, 1.031985d0, 1.037751d0, 1.043662d0, 1.049726d0, &
            1.055956d0, 1.062364d0, 1.068965d0, 1.075780d0, 1.082830d0, &
            1.090140d0, 1.097737d0, 1.105647d0, 1.113894d0, 1.122497d0, 1.131470d0 /)

        x1 = 1.0d0 / x
        i = int(200.0d0*x1)
        if(i.lt.0) i = 0
        if(i.gt.19) i = 19
        i1 = i + 1
        rih2 = f(i) + 200.0d0*(f(i1)-f(i))*(x1 - 0.005d0*dble(i))
    
    endfunction



    double precision function ei(x)

        use constants, only: euler

        integer                         ::  i
        double precision                ::  xp, si
        double precision, intent(in)    ::  x
        double precision, parameter     ::  a1 = 8.5733287401d0, a2 = 18.0590169730d0, a3 = 8.6347608925d0, a4 = 0.2677737343d0
        double precision, parameter     ::  b1 = 9.5733223454d0, b2 = 25.6329561486d0, b3 = 21.0996530827d0, b4 = 3.9584969228d0

        if(x.gt.60.0d0) stop 'error: ei argument too large'
        if(x.lt.-60.0d0) then
            ei = 0.0d0
            return
        endif
        if(abs(x).lt.1.0d-14) then
            ei = 0.0d0
            return
        endif
        if(x.lt.-1.0d0) then
            xp = abs(x)
            ei = -(a4+xp*(a3+xp*(a2+xp*(a1+xp)))) / (b4+xp*(b3+xp*(b2+xp*(b1+xp)))) * exp(-xp) / xp
            return
        endif

        ei = euler + log(abs(x)) + x
        i = 1
        si = x
        do
            si = si * x * i / (i + 1.0d0)**2
            ei = ei + si
            if (isnan(ei)) then
                write(*,*) 'error: NaN in ei(), x=', x
                stop
            endif
            
            if(abs(si/x).le.1.0d-6) exit
            i = i + 1

            if (i .gt. 1000) then
                write(*,*) 'error: ei() did not converge, x=', x
                stop
            endif
        enddo
    
    endfunction


endmodule