module constants

    implicit none

    integer, parameter              ::  record_unit = 11
    integer, parameter              ::  base_seed = 1357911
    integer, parameter              ::  aberr_periodicity(24) = (/ &
                                            & 1, 0, 2, 1, 3, 0, 2, 4, 1, 3, 5, 0, &
                                            & 2, 4, 6, 1, 3, 5, 7, 0, 2, 4, 6, 8 /)
    integer, parameter              ::  nx_max = 256, ny_max = 256, nz_max = 20
    integer, parameter              ::  n_atoms_max = 5000, n_atoms_substrate_max = 35000, n_types_max = 5
    integer, parameter              ::  n_atoms_tot_max = n_atoms_max + n_atoms_substrate_max
    double precision, parameter     ::  interatomic_distance = 2.889d0
    double precision, parameter     ::  cluster_cutoff = 1.5d0 * interatomic_distance
    double precision, parameter     ::  pi = acos(-1.0d0)
    double precision, parameter     ::  lateral_border_margin = 2.0d0
    double precision, parameter     ::  dbf_ag = 0.019d0, dbf_co = 0.015d0, dbf_c = 0.003d0 ! nm², voir README.md
    double precision, parameter     ::  hc = 1.2398419843320026d0 ! Planck's constant times speed of light
    double precision, parameter     ::  e0 = 510.99895069d0, sigma0 = 2.0886573497d0
    double precision, parameter     ::  r8pi2 = 1.0d0/(8.0d0*pi*pi)
    double precision, parameter     ::  euler = 0.5772156649015328606065120d0
    double precision, parameter     ::  v0 = 0.03809982119d0
    double precision, parameter     ::  box_hrtem(3) = (/ 10.0d0, 10.0d0, 15.0d0 /)
    character(len=255), parameter   ::  file_list = "xyz_file_list.tmp"
    character(len=255), parameter   ::  substrate_list = "substrate_list.tmp"
    
endmodule


module variable

    use constants, only: n_atoms_max, n_atoms_substrate_max, n_atoms_tot_max, n_types_max, nx_max, ny_max, nz_max

    implicit none

    integer             ::  max_files, snapshot_index, augmentation_index, seed, doptc, dopsc, dovib, n_ranks, nx, ny, nz, n_seed
    integer             ::  atomic_number(n_atoms_tot_max), n_atoms_substrate, n_atoms_tot
    integer             ::  n_substrate_files
    double precision    ::  box(3), shift(3), biso(n_atoms_tot_max), image(nx_max, ny_max), pos(3, n_atoms_max), epot(n_atoms_max)
    double precision    ::  pos_cluster(3, n_atoms_tot_max), readout_noise_e, aberr_re(24), aberr_im(24), position_scale
    double precision    ::  ht, fs, edge, sc_mrad, vib1, vib2, vibdir, oapr, dose_e_per_a2, lambda, g2, dx, dy, dz, gmax
    double precision    ::  phi, theta, psi
    double precision    ::  pos_substrate(3, n_atoms_substrate_max)
    character(len=2)    ::  species(n_atoms_tot_max), atom_typ1, atom_typ2, species_substrate(n_atoms_substrate_max)
    character(len=255)  ::  xyz_files(48000), data_file, output_data, substrate_file
    character(len=255)  ::  substrate_dir
    character(len=:), allocatable   ::  substrate_files(:) ! noms des fichiers de substrat (sans le dossier)
    double complex      ::  trans(nx_max, ny_max, nz_max), wave(nx_max, ny_max), wave_fft(nx_max, ny_max)
    logical             ::  placed, found


endmodule





module descriptor

    implicit none

    integer             :: id_sim               ! Simulation ID
    integer             :: id_sim_bis           ! Simulation ID (after rotation)
    integer             :: n_atoms              ! Number of atoms
    double precision    :: composition          ! Composition
    integer             :: n_steps              ! Number of steps
    double precision    :: initial_temperature  ! Initial temperature
    double precision    :: epot_total           ! Total potential energy
    double precision    :: gyration_radius      ! Gyration radius
    integer             :: nat1                 ! Number of atoms of type 1
    integer             :: nat2                 ! Number of atoms of type 2
    integer             :: nat1_out             ! Number of atoms of type 1 outside the virtual sphere
    integer             :: nat2_out             ! Number of atoms of type 2 outside the virtual sphere
    integer             :: nat1_in              ! Number of atoms of type 1 inside the virtual sphere
    integer             :: nat2_in              ! Number of atoms of type 2 inside the virtual sphere
    integer             :: n_clusters           ! Number of clusters
    double precision    :: d_com                ! Distance between centers of mass of the two clusters
    double precision    :: coreshell_index      ! Core-shell index
    double precision    :: mass_center(3)
    double precision    :: mass_center_nat1(3)
    double precision    :: mass_center_nat2(3)

endmodule