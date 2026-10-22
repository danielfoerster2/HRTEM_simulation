module wavimg

    use variable, only: nx_max, ny_max

    implicit none

    double complex      ::  image_c(nx_max, ny_max)
    double precision    ::  wx, wy, dchix, dchiy, coef(0:8, 24)
    logical             ::  aberr_on(24)
    integer, parameter  ::  am(24) = (/ 1, 2, 2, 3, 3, 4, 4, 4, 5, 5, 5, 6, 6, 6, 6, 7, 7, 7, 7, 8, 8, 8, 8, 8 /)
    integer, parameter  ::  an(24) = (/ 1, 0, 2, 1, 3, 0, 2, 4, 1, 3, 5, 0, 2, 4, 6, 1, 3, 5, 7, 0, 2, 4, 6, 8 /)

    contains

    subroutine run_wavimg
    
        use variable, only: wave, lambda, dx, nx, ny, wave_fft, image, oapr
        use fft, only: fft2


        if (oapr .le. 0.0d0) oapr = 1000.0d0 * lambda * 0.5d0/dx
        call fft2(wave(1:nx, 1:ny), wave_fft(1:nx, 1:ny), -1)
        call explicit_focus_image
        image(1:nx, 1:ny) = dble(image_c(1:nx, 1:ny))
        call apply_vibration
        call apply_detector_model

    endsubroutine


    subroutine apply_vibration
        
        use constants, only: pi
        use variable, only: image, dx, dy, vib1, vib2, vibdir, dovib, nx, ny
        use fft, only: fft2

        integer             ::  nnx, nny, j_px, i_px
        double precision    ::  v12, v22, cd, sd, py, py2, px, px2, vtmp, vmean
        double complex      ::  cvib(nx, ny), image_fft(nx, ny), image_complex(nx, ny)

        if (dovib .eq. 0 .or. vib1 .le. 0.0d0) return
        nnx = nx / 2
        nny = ny / 2
        v12 = vib1**2
        if (dovib .eq. 1) then
            v22 = v12
        else
            v22 = vib2**2
        endif
        if (v22 .le. 0.0d0) v22 = v12
        cd = cos(pi*vibdir/180.0d0)
        sd = sin(pi*vibdir/180.0d0)
        cvib = dcmplx(0.0d0, 0.0d0)
        do j_px = 1, ny
            if (j_px .le. nny) then
                py = dy*(j_px-1)
            else
                py = dy*(j_px-1-ny)
            endif
            py2 = py**2
            do i_px = 1, nx
                if (i_px .le. nnx) then
                    px = dx * (i_px-1)
                else
                    px = dx * (i_px-1-nx)
                endif
                px2 = px**2
                vtmp = exp(-((py2*v12 + px2*v22)*cd*cd + ((px2*v12 + py2*v22)*sd &
                        - 2.0d0*px*py*(v12-v22)*cd)*sd) / (2.0d0*v12*v22))
                cvib(i_px, j_px) = dcmplx(vtmp, 0.0d0)
            enddo
        enddo

        call fft2(cvib, image_fft, -1)
        vmean = dble(image_fft(1, 1))
        cvib = image_fft / vmean
        image_complex = dcmplx(image(1:nx, 1:ny), 0.0d0)
        call fft2(image_complex, image_fft, -1)
        image_fft = image_fft * conjg(cvib)
        call fft2(image_fft, image_complex, 1)
        image_complex = image_complex / (nx * ny)
        image(1:nx, 1:ny) = dble(image_complex)

    endsubroutine

    subroutine apply_detector_model

        use variable, only: nx, ny, image, dx, dy, dose_e_per_a2, readout_noise_e

        integer             ::  i_px, j_px
        double precision    ::  e_per_pixel, pixel_area_a2, scaled

        pixel_area_a2 = 100.0d0 * dx * dy
        e_per_pixel = dose_e_per_a2 * pixel_area_a2
        if(e_per_pixel.le.0.0d0) return
        do j_px = 1, ny
            do i_px = 1, nx
                scaled = max(image(i_px, j_px), 0.0d0) * e_per_pixel
                scaled = dble(sample_poisson(scaled))
                if(readout_noise_e.gt.0.0d0) then
                    scaled = scaled + readout_noise_e * sample_standard_normal()
                endif
                image(i_px, j_px) = scaled / e_per_pixel
            enddo
        enddo

    endsubroutine

    subroutine explicit_focus_image

        use variable, only: wave_fft, doptc, fs, lambda, nx, ny, dx, dy, edge, dopsc, sc_mrad, oapr
        use constants, only: pi
        use fft, only: fft_index, fft2

        integer             ::  nkfs, j_px, i_px, k
        double precision    ::  info_width, df_step, wsum, defocus_prefac, gy, gx, g2, df, wt, chi, sc_rad
        double precision    ::  w2tab(nx, ny), wa_mrad, aperture, chi0(nx, ny), envs, dchix0, dchiy0, scpf
        double complex      ::  transfer0(nx, ny), image_wave_fft(nx, ny), image_wave_rs(nx, ny)

        image_c = dcmplx(0.0d0, 0.0d0)
        nkfs = 0
        if (doptc .ne. 0 .and. fs .gt. 0.0d0) then
            info_width = ((pi*lambda*fs)**2 / 8.0d0)**0.25d0
            df_step = 2.0d0 * info_width**2 / (9.0d0*lambda)
            nkfs = max(1, 1 + int(2.0d0*fs / df_step))
        endif
        sc_rad = sc_mrad * 0.001d0
        scpf = (0.5d0*sc_rad) **2
        call aberration_coefficients
        wx = 0.0d0
        wy = 0.0d0
        call aberration_gradient
        dchix0 = dchix
        dchiy0 = dchiy
        wsum = 0.0d0
        defocus_prefac = pi / lambda
        do j_px = 1, ny
            gy = fft_index(j_px, ny) / (dy*ny)
            wy = lambda * gy
            do i_px = 1, nx
                gx = fft_index(i_px, nx) / (dx*nx)
                wx = lambda * gx
                g2 = gx**2 + gy**2
                w2tab(i_px, j_px) = wx**2 + wy**2
                wa_mrad = 1000.0d0 * lambda * sqrt(g2)
                if (edge .gt. 0.0d0) then
                    aperture = 0.5d0 - 0.5d0*tanh(pi*(wa_mrad - oapr) / (oapr*edge))
                else
                    aperture = 1.0d0
                    if (wa_mrad .gt. oapr) aperture = 0.0d0
                endif
                chi0(i_px, j_px) = aberration_chi()
                envs = 1.0d0
                if (dopsc .ne. 0) then
                    call aberration_gradient
                    envs = exp(-scpf*((dchix-dchix0)**2 + (dchiy-dchiy0)**2))
                endif
                transfer0(i_px, j_px) = aperture * envs * dcmplx(cos(chi0(i_px, j_px)), -sin(chi0(i_px, j_px)))
            enddo
        enddo

        do k = -nkfs, nkfs
            df = 0.0d0
            wt = 1.0d0
            if (nkfs .gt. 0) then
                df = 2.0d0 * fs * k / nkfs
                wt = exp(-(df/fs)**2)
            endif
            do j_px = 1, ny
                do i_px = 1, nx
                    chi = defocus_prefac * df * w2tab(i_px, j_px)
                    image_wave_fft(i_px, j_px) = wave_fft(i_px, j_px) * transfer0(i_px, j_px) * dcmplx(cos(chi), -sin(chi))
                enddo
            enddo
            call fft2(image_wave_fft, image_wave_rs, 1)
            image_wave_rs = image_wave_rs / (nx*ny)
            image_c(1:nx, 1:ny) = image_c(1:nx, 1:ny) + wt * image_wave_rs * conjg(image_wave_rs)
            wsum = wsum + wt
        enddo
        image_c = image_c / wsum

    endsubroutine

    subroutine aberration_gradient
        
        use constants, only: pi
        use variable, only: lambda

        integer             ::  k, l, m, n, j, p
        double precision    ::  w, prefac, poly, dpolyx, dpolyy, radial_x, radial_y
        double precision    ::  wxpow(0:8), wypow(0:8), wpow(0:8)

        w = sqrt(wx*wx + wy*wy)
        wxpow(0) = 1.0d0
        wypow(0) = 1.0d0
        wpow(0) = 1.0d0
        do j = 1, 8
            wxpow(j) = wxpow(j-1) * wx
            wypow(j) = wypow(j-1) * wy
            wpow(j) = wpow(j-1) * w
        enddo
        dchix = 0.0d0
        dchiy = 0.0d0
        do k = 1, 24
            if (.not. aberr_on(k)) cycle
            m = am(k)
            n = an(k)
            p = m - n
            poly = 0.0d0
            dpolyx = 0.0d0
            dpolyy = 0.0d0
            do l = 0, n
                j = n - l
                poly = poly + coef(l, k) * wxpow(j) * wypow(l)
                if(j.gt.0) then
                    dpolyx = dpolyx + coef(l, k) * dble(j) * wxpow(j-1) * wypow(l)
                endif
            enddo
            do l = 1, n
                j = n - l
                dpolyy = dpolyy + coef(l, k) * dble(l) * wxpow(j) * wypow(l-1)
            enddo
            radial_x = 0.0d0
            radial_y = 0.0d0
            if(p.gt.0 .and. w.gt.0.0d0) then
                radial_x = dble(p) * wpow(p-2) * wx * poly
                radial_y = dble(p) * wpow(p-2) * wy * poly
            endif
            dchix = dchix + (radial_x + wpow(p)*dpolyx) / dble(m)
            dchiy = dchiy + (radial_y + wpow(p)*dpolyy) / dble(m)
        enddo
        prefac = 2.0d0*pi / lambda
        dchix = prefac * dchix
        dchiy = prefac * dchiy
    
    endsubroutine


    subroutine aberration_coefficients

        use variable, only: aberr_re, aberr_im

        integer ::  k, l

        ! Coefficients des polynomes d'aberration, independants du pixel: calcules une fois par image
        do k = 1, 24
            aberr_on(k) = .not. (aberr_re(k)*aberr_re(k) + aberr_im(k)*aberr_im(k) .le. 1.0d-30)
            coef(:, k) = 0.0d0
            do l = 0, an(k)
                coef(l, k) = (sign_cycle(l) * aberr_re(k) + sign_cycle(l+3) * aberr_im(k)) * dble(binomial(an(k), l))
            enddo
        enddo

    endsubroutine


    double precision function sign_cycle(i)

        integer, intent(in) ::  i
    
        select case(modulo(i, 4))
        case(0)
            sign_cycle = 1.0d0
        case(2)
            sign_cycle = -1.0d0
        case default
            sign_cycle = 0.0d0
        endselect

    endfunction


    integer function binomial(n, k)

        integer, intent(in) ::  n, k
        integer             ::  i
    
        binomial = 1
        do i = 1, k
            binomial = binomial * (n + 1 - i) / i
        enddo
    
    endfunction


    double precision function aberration_chi()

        use constants, only: pi
        use variable, only: lambda

        integer             ::  k, l, m, n, j
        double precision    ::  w, term, prefac
        double precision    ::  wxpow(0:8), wypow(0:8), wpow(0:8)

        w = sqrt(wx*wx + wy*wy)
        wxpow(0) = 1.0d0
        wypow(0) = 1.0d0
        wpow(0) = 1.0d0
        do j = 1, 8
            wxpow(j) = wxpow(j-1) * wx
            wypow(j) = wypow(j-1) * wy
            wpow(j) = wpow(j-1) * w
        enddo
        aberration_chi = 0.0d0
        do k = 1, 24
            if (.not. aberr_on(k)) cycle
            m = am(k)
            n = an(k)
            term = 0.0d0
            do l = 0, n
                j = n - l
                term = term + coef(l, k) * wxpow(j) * wypow(l)
            enddo
            term = term * wpow(m-n) / dble(m)
            aberration_chi = aberration_chi + term
        enddo
        prefac = 2.0d0*pi / lambda
        aberration_chi = prefac * aberration_chi
    
    endfunction


    integer function sample_poisson(lambda)

        integer                         ::  k
        double precision                ::  limit, p, u
        double precision, intent(in)    ::  lambda

        if(lambda.le.0.0d0) then
            sample_poisson = 0
            return
        endif
        if(lambda.lt.30.0d0) then
            limit = exp(-lambda)
            p = 1.0d0
            k = 0
            do
                call random_number(u)
                p = p * u
                if(p.le.limit) exit
                k = k + 1
            enddo
            sample_poisson = k
            return
        endif
        sample_poisson = max(0, nint(lambda + sqrt(lambda) * sample_standard_normal()))
    
    endfunction

    
    double precision function sample_standard_normal()

        use constants, only: pi 

        double precision    ::  u1, u2

        call random_number(u1)
        call random_number(u2)
        u1 = max(u1, tiny(u1))
        sample_standard_normal = sqrt(-2.0d0*log(u1)) * cos(2.0d0*pi*u2)
    
    endfunction

endmodule