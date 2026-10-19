# HRTEM_simulation

Multislice HRTEM image simulation of MD nanoparticles (AgCo) on amorphous-carbon substrates, parallelised with MPI.
Settings and paths are in `run.sh`; `make` builds `main`.

## Requirements

- An MPI Fortran compiler (`mpifort`), FFTW3 (double precision) and OpenBLAS. The library paths are set in `LDLIBS`
  in the `Makefile` (`/lib/x86_64-linux-gnu/libfftw3.so.3` and `/lib/x86_64-linux-gnu/libopenblas.so.0`, from the
  Ubuntu packages `libfftw3-double3` and `libopenblas0`).
- `run.sh` sets `OPENBLAS_NUM_THREADS=1`. Each MPI rank must use a single BLAS thread; otherwise every rank starts one
  thread per core and the machine is oversubscribed.

## Debye-Waller factors

The isotropic Debye-Waller factors are set in `constants.f90`, in nm² (1 nm² = 100 Å²):

```fortran
double precision, parameter :: dbf_ag = 0.019d0, dbf_co = 0.015d0, dbf_c = 0.003d0 ! nm²
```

Convention: B = 8π² ⟨u_x²⟩, where ⟨u_x²⟩ is the mean-square displacement along one axis. The elastic form factor of each
atom is damped by exp(-B g²/4), with g in nm⁻¹ (`scattering_factor` in `slc.f90`). B also sets the absorptive part of the
form factor (thermal diffuse scattering), so a larger B means both weaker high-frequency contrast and more absorption.

| Element | In the code | Bulk value at 300 K (Debye model) | Status |
|---|---|---|---|
| Ag | 0.019 nm² (1.9 Å²) | 0.64 Å² (θ_D = 225 K) | to be discussed |
| Co | 0.015 nm² (1.5 Å²) | 0.31 Å² (θ_D = 445 K) | to be discussed |
| C (a-C substrate) | 0.003 nm² (0.3 Å²) | see below | proposed |

### Carbon (substrate)

Proposed value: **B = 0.3 Å² (0.003 nm²)**, one value for all substrate atoms.

- The a-C slabs come from classical MD. Their per-atom `dwf` column has slab means of 0.14-0.22 Å². Classical MD has no
  zero-point motion, which is large for carbon at 300 K because its vibrational frequencies are far above k_B T.
- A Debye model fitted to the classical slab means (effective θ_D = 1100-1400 K) gives a quantum correction of ×1.4-1.5,
  i.e. 0.21-0.30 Å². This is a lower bound: the high-frequency sp² modes carry most of the zero-point motion, and for
  diamond the same model gives ×2.1 (0.058 Å² classical, 0.12 Å² with zero-point motion).
- Measured room-temperature values are about 0.14-0.2 Å² for diamond and 0.2-0.7 Å² for graphite (anisotropic). An
  a-C film (mostly sp², 1.8-2.3 g/cm³) lies between the two.
- The image is not very sensitive to this value: at the 256 px band limit (8.53 nm⁻¹, 0.117 nm) the damping is 0.947
  for 0.3 Å² and 0.971 for 0.16 Å².

The substrate file reader now reads only species and positions; the per-atom `dwf` column is ignored.

### Silver and cobalt: open questions

The current values are about 3× (Ag) and 5× (Co) the bulk room-temperature values. Two questions decide what they
should be:

1. **Where do 1.9 Å² and 1.5 Å² come from?** They could be meant to include surface softening or the MD temperature.
2. **Are the MD snapshot positions instantaneous or time-averaged?**
   - Instantaneous positions already contain the thermal displacements: each image is then one frozen-phonon
     configuration, and adding a B on top counts the thermal motion twice. B should be close to 0 in that case
     (averaging over several snapshots then plays the role of the Debye-Waller damping).
   - Time-averaged (equilibrium) positions need the full B. Bulk values are the starting point; surface atoms vibrate
     more than bulk atoms, so a somewhat larger value is defensible for small particles.

Effect of the choice on the damping factor exp(-B g²/4):

| Reflection | g (nm⁻¹) | Ag, 1.9 Å² | Ag, 0.64 Å² | Co, 1.5 Å² | Co, 0.31 Å² |
|---|---|---|---|---|---|
| Ag 111 (0.236 nm) | 4.24 | 0.918 | 0.972 | | |
| Co 111, fcc (0.205 nm) | 4.89 | | | 0.914 | 0.982 |
| Ag 220 (0.144 nm) | 6.92 | 0.796 | 0.926 | | |
| Co 220, fcc (0.125 nm) | 7.98 | | | 0.788 | 0.952 |
| band limit, 256 px (0.117 nm) | 8.53 | 0.708 | 0.890 | 0.761 | 0.945 |

### Form factor cache: one B per element

The form factors (Weickenmeier-Kohl, with absorptive part) are computed once per run for each distinct (Z, B) pair and
reused for every image (`unique_scattering_factors` and `compute_form_factor` in `slc.f90`). This is valid because the
grid, the voltage and B are fixed for the whole run. At 256 px the computation takes 1-2 s per pair and per rank,
once per run instead of once per image.

Consequences:

- **B must be a single value per element.** The cache has `n_types_max` = 5 slots. A per-slab or per-atom B would fill
  it and stop the run with `more than n_types_max distinct (Z, B) pairs`. A per-atom B would need a different scheme
  (for example a few B bins per element).
- To change a value, edit `constants.f90` and rebuild. Nothing else depends on it.
