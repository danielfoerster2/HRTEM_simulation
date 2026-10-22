# HRTEM_simulation

Multislice HRTEM image simulation of MD nanoparticles (AgCo) on amorphous-carbon substrates, parallelised with MPI.
Settings and paths are in `run.sh`; `make` builds `main`.

## Requirements

- An MPI Fortran compiler (`mpifort`), FFTW3 (double precision) and OpenBLAS. The library paths are set in `LDLIBS`
  in the `Makefile` (`/lib/x86_64-linux-gnu/libfftw3.so.3` and `/lib/x86_64-linux-gnu/libopenblas.so.0`, from the
  Ubuntu packages `libfftw3-double3` and `libopenblas0`).
- `run.sh` sets `OPENBLAS_NUM_THREADS=1`. Each MPI rank must use a single BLAS thread; otherwise every rank starts one
  thread per core and the machine is oversubscribed.

## Output

A run writes a single file, `output_data` in `run.sh` (`dataset.npy`). It is a NumPy file holding one record per
image, and each record contains the descriptors and the image of that simulation. Both are written by one binary
`write` statement (`save_data_row` in `utils_io.f90`), so a descriptor row cannot be paired with the wrong image.

Fields of a record, in file order:

| Fields | Type |
|---|---|
| `id_sim n_atoms n_steps nat1 nat2 nat1_out nat2_out nat1_in nat2_in` | int32 |
| `initial_temperature epot_total composition gyration_radius d_com coreshell_index` | float64 |
| `phi theta psi position_scale` (rotation angles in rad, scale factor) | float64 |
| `defocus A1_re A1_im B2_re B2_im A2_re A2_im C3 S3_re S3_im A3_re A3_im B4_re B4_im D4_re D4_im A4_re A4_im` (nm) | float64 |
| `image` (ny × nx, first index = y, min-max scaled to 0-255) | uint8 |

`defocus` is C1. `d_com` is NaN when the cluster contains a single species. A record takes 260 + n² bytes: 64.3 KB
at 256 px. The field names and their order are set in `int_fields` and `real_fields` (`utils_io.f90`).

```python
import numpy as np
data = np.load("dataset.npy", mmap_mode="r")      # memory-mapped, nothing is read until used
data["image"][0]                                  # first image, shape (256, 256), dtype uint8
data[0]["phi"], data[0]["defocus"]                # descriptors of the same record
data.dtype.names                                  # all field names

import pandas as pd                               # descriptors as a table
descriptors = pd.DataFrame({name: data[name] for name in data.dtype.names if name != "image"})
```

Reading one descriptor for all images touches every record of the file, which is slow for a large dataset. If the
descriptors are needed often, extract them once (as above) and save that table; `id_sim` links it back to the
records.

`images_to_png.py` writes one `{id_sim}.png` per image if PNG files are needed. For training, reading the
memory-mapped array directly avoids creating one file per image.

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
| Ag | 0.019 nm² (1.9 Å²) | 0.64 Å² (θ_D = 225 K) | 3× too large, to be discussed (see below) |
| Co | 0.015 nm² (1.5 Å²) | 0.31 Å² (θ_D = 445 K) | 3× too large, to be discussed (see below) |
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

#### Substrate positions are time-averaged (verified)

The slabs are made with LAMMPS (Tersoff carbon, melt-quench to 300 K; `in.substrate` / `in.substrate.old-lammps` and
`run.sh` of the `carbon_substrate` generator). Their positions are time averages, so applying B to them does not count
the thermal motion twice:

- During the last 500 fs at 300 K, `fix ave/atom 1 500 500` averages the unwrapped coordinates and their squares at
  every step. The file holds the averaged positions ⟨r⟩ and, per atom, B = 8π²/3 (⟨r²⟩ - ⟨r⟩²), i.e. with the factor
  1/3 that turns the 3D mean-square displacement into ⟨u_x²⟩.
- Checked on a 5 x 5 x 2 nm test run that saved a frame every 10 fs through the averaging window: the positions in
  the output agree with an independent average of the raw frames to 0.001 Å (rms) and differ from the instantaneous
  last frame by 0.092 Å (rms). Mean B: 0.1986 Å² in the file, 0.1987 Å² recomputed, per-atom correlation 1.000.
- Checked on 21 production slabs: each file holds the averaged frame (not the empty step-0 frame LAMMPS also writes),
  the `dwf_mean` in the header matches the column, z lies in [0, thickness], and x, y go up to 0.18 Å outside
  [0, 100) (unwrapped averages; the HRTEM code wraps them when it applies the random lateral shift).
- So the slab means of 0.15-0.20 Å² are correct classical values, and the 0.3 Å² used here adds the zero-point motion.
- Minor: 0.02-0.35 % of the atoms hop between two sites during the window (B > 1 Å²) and their averaged position lies
  between the two; the 500 fs window does not fully average the slowest in-plane vibrations (period about 1 ps at
  10 nm wavelength). Both are negligible at 0.12 nm resolution.

### Silver and cobalt: open questions

1. **The current values are 3× too large.** They come from `debye.py` (repository `vbarrere/Debye-Waller_factor`),
   which computes B = 8π² ⟨|r - ⟨r⟩|²⟩ from a series of MD frames. ⟨|r - ⟨r⟩|²⟩ sums the three axes, whereas
   B = 8π² ⟨u_x²⟩ = 8π² ⟨|r - ⟨r⟩|²⟩ / 3. Divided by 3, Ag gives 0.0063 nm² (0.63 Å²), which matches the classical bulk
   value for Ag at 300 K (0.63 Å²), and Co gives 0.005 nm² (0.5 Å², assuming it came from the same script), somewhat
   above the bulk estimate of about 0.3 Å². (Minor: the loop over frames starts at the second frame but divides by
   the total number of frames.)
2. **Are the MD snapshot positions instantaneous or time-averaged?** Unlike the substrate, everything points to
   instantaneous positions: the MD code writes instantaneous frames (`debye.py` relies on their fluctuations), the
   dataset files have the same layout with a per-atom potential energy, and nothing in the workflow averages positions.
   What is not known here is the temperature at the end of the runs.
   - Instantaneous positions at about 300 K already contain the thermal displacements: each image is then one
     frozen-phonon configuration, and adding a B on top counts the thermal motion twice. B should be close to 0 in
     that case (averaging over several snapshots then plays the role of the Debye-Waller damping).
   - Time-averaged positions, or a quench that ends near 0 K, need the full B: the corrected values, 0.0063 nm² (Ag)
     and 0.005 nm² (Co). Surface atoms vibrate more than bulk atoms, so a somewhat larger value is defensible for
     small particles.
   - A check on one snapshot: for like-species nearest neighbours in the core (away from the surface and the Ag/Co
     interface), the spread of the pair distance is about 0.07-0.1 Å for instantaneous positions at 300 K, and below
     about 0.03 Å for averaged positions or a near-0 K quench.

Effect of the choice on the damping factor exp(-B g²/4) (B = 0 gives 1):

| Reflection | g (nm⁻¹) | Ag, 1.9 Å² (now) | Ag, 0.63 Å² (corrected) | Co, 1.5 Å² (now) | Co, 0.5 Å² (corrected) |
|---|---|---|---|---|---|
| Ag 111 (0.236 nm) | 4.24 | 0.918 | 0.972 | | |
| Co 111, fcc (0.205 nm) | 4.89 | | | 0.914 | 0.971 |
| Ag 220 (0.144 nm) | 6.92 | 0.796 | 0.927 | | |
| Co 220, fcc (0.125 nm) | 7.98 | | | 0.788 | 0.923 |
| band limit, 256 px (0.117 nm) | 8.53 | 0.708 | 0.892 | 0.761 | 0.913 |

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
