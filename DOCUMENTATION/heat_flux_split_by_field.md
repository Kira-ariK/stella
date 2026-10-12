# Heat flux split by field (ExB / flutter / compressional)

This change makes stella write the turbulent heat flux separated into the
contribution from each field, in addition to the total:

| Piece | Field | Physical meaning |
|---|---|---|
| `phi`  | φ   | ExB (electrostatic) heat flux |
| `apar` | A∥  | Magnetic-flutter heat flux |
| `bpar` | B∥  | Compressional (δB∥) heat flux |

At every output time the three pieces add up exactly to the existing total
`qflux_vs_s` (and `qflux_vs_kxkys`):

```
qflux = qflux_phi + qflux_apar + qflux_bpar
```

The total heat flux itself is computed exactly as before. Existing outputs,
including the `.fluxes` text file, are unchanged.

---

## New NetCDF variables

These are written to `<run>.out.nc`.

### Time traces

| Variable | Dimensions (Python order) | Written when |
|---|---|---|
| `qflux_phi_vs_s`  | `(t, species)` | always (flux-tube runs) |
| `qflux_apar_vs_s` | `(t, species)` | `include_apar = .true.` |
| `qflux_bpar_vs_s` | `(t, species)` | `include_bpar = .true.` |

### (kx, ky) spectra

| Variable | Dimensions (Python order) | Written when |
|---|---|---|
| `qflux_phi_vs_kxkys`  | `(t, species, kx, ky)` | `write_fluxes_kxky = .true.` |
| `qflux_apar_vs_kxkys` | `(t, species, kx, ky)` | `write_fluxes_kxky = .true.` and `include_apar = .true.` |
| `qflux_bpar_vs_kxkys` | `(t, species, kx, ky)` | `write_fluxes_kxky = .true.` and `include_bpar = .true.` |

**Units** are the same as `qflux_vs_s`: `n_ref * T_ref * v_ref * (rho_ref/a_ref)^2`,
with `v_ref = sqrt(2 T_ref / m_ref)`.

Normalisation follows the existing total. That means the same `flux_norm`
factor, the same flux-surface or field-line average, the same `ñ_s T̃_s` factor
and the same division by `ntubes`. As with `qflux_vs_kxkys`, the spectra already
include the reality condition, so summing over all (kx, ky) gives the mode
total. The A∥ and B∥ variables are omitted, not written as zeros, when that
field is switched off. In an electrostatic run, therefore, only the `phi`
piece appears, and it equals the total.

---

## Input file

You don't need a new input flag; the split time traces are written whenever
the fluxes are. To get the split (kx, ky) spectra as well, turn on the existing
spectrum output:

```fortran
&diagnostics_fluxes
  write_fluxes_vs_time = .true.   ! qflux_*_vs_s
  write_fluxes_kxky    = .true.   ! qflux_*_vs_kxkys
/
```

Output cadence follows the usual `&diagnostics` settings (`nwrite`, `nc_mult`).

---

## Reading the output

```python
import netCDF4 as nc
import numpy as np

d = nc.Dataset("input.out.nc")
t = d["t"][:]

q_tot  = d["qflux_vs_s"][:]                   # (t, species)
q_exb  = d["qflux_phi_vs_s"][:]
q_flut = d["qflux_apar_vs_s"][:] if "qflux_apar_vs_s" in d.variables else 0 * q_exb
q_bpar = d["qflux_bpar_vs_s"][:] if "qflux_bpar_vs_s" in d.variables else 0 * q_exb

assert np.allclose(q_tot, q_exb + q_flut + q_bpar)

# (kx, ky) spectrum of the ExB heat flux, species 0, last time
Q_exb_kxky = d["qflux_phi_vs_kxkys"][-1, 0]     # (kx, ky)
```

With xarray: `xr.open_dataset("input.out.nc")["qflux_phi_vs_kxkys"]` gives a
labelled `(t, species, kx, ky)` array.

---

## Building

The changes are in three source files:

- `STELLA_CODE/diagnostics/diagnostics_fluxes_fluxtube.f90`: computes the
  heat flux from each field into its own array, and takes the total as their
  sum.
- `STELLA_CODE/diagnostics/diagnostics_fluxes.f90`: stores the split arrays
  and passes them to the NetCDF writer.
- `STELLA_CODE/diagnostics/write_diagnostics_to_netcdf.fpp`: new routines
  `write_qflux_byfield_vs_time_nc` and `write_qflux_byfield_kxkys_nc`.

Build stella as usual. **Watch out:** after linking, `COMPILATION/CMakeLists.txt`
creates symlinks at `<build>/COMPILATION/../../stella` and `../../../stella`.
Those paths assume the build directory is `COMPILATION/build_cmake`. A build
directory placed directly in the repository root (e.g. `build_x/`) therefore
**overwrites the top-level `stella` symlink**. To make a test build without
touching the production executable, nest it two levels deep:

```bash
source COMPILATION/LoadModules/load_modules.pitagora
module load cmake
cmake -S . -B build_test/sub/build -DCMAKE_BUILD_TYPE=RelWithDebInfo \
      -DCMAKE_Fortran_COMPILER=$(which mpifort) -DSTELLA_ENABLE_NETCDF=ON
cmake --build build_test/sub/build -j 16
# executable: build_test/stella  (build_*/ is gitignored)
```

---

## Verification

These tests used the electromagnetic nonlinear test input
(`AUTOMATIC_TESTS/numerical_tests/test_7_electromagnetic/EM_kxky_grid_box_nonlinear.in`,
with `include_apar = include_bpar = .true.` and `beta = 0.04`). I added
`write_fluxes_kxky = .true.` and ran 10 steps on 2 MPI ranks.

| Check | Result |
|---|---|
| `qflux_vs_s` vs `phi + apar + bpar` | identical (relative difference 0) |
| `qflux_vs_kxkys` vs sum of the split spectra | 4 × 10⁻¹⁶ (rounding) |
| New `qflux_vs_s` vs the previous executable | 1 × 10⁻¹⁵ (rounding) |
| New `qflux_vs_kxkys`, `pflux_vs_s`, `vflux_vs_s`, `phi2` vs the previous executable | bit-identical |
| Electrostatic copy of the input | only `qflux_phi_*` written; equals the total exactly |

---

## Limitations

1. **Flux-tube only.** Full-flux-surface runs (`full_flux_surface = .true.`)
   use a different flux routine, and they don't write the split variables.
2. **Radial-variation runs aren't validated.** With `radial_variation = .true.`
   the split variables are still written by the flux-tube routine, but they
   haven't been tested there. They use `phi` without the `phi_corr_QN`
   correction that the radial-variation fluxes (`qflux_x`) include.
3. **Only the heat flux is split.** The particle flux (`pflux_*`) and momentum
   flux (`vflux_*`) are still totals only. They could be split the same way.
4. **No split z-resolved output.** `qflux_vs_kxkyzs` (`write_fluxes_kxkyz`) is
   not split by field.
5. **The text output is unchanged.** `.fluxes` still holds only the totals,
   so the split is available in NetCDF only.
6. **The B∥ convention differs from GENE.** stella reports B∥ as a separate
   piece, whereas GENE folds the B∥ contribution into its "electromagnetic"
   flux. To compare with GENE, compare `qflux_apar + qflux_bpar` to GENE's
   EM heat flux.
7. **Restarts aren't tested.** On a restart that appends to an existing
   `.out.nc` written by an older executable, the new variables are created at
   the first write. Earlier time indices then hold NetCDF fill values.
8. **Post-processing tools don't read it yet.** stellapy and pyrokinetics'
   stella reader don't use the new variables yet. Read them directly as shown
   above.
9. **Testing was limited.** It covers one small electromagnetic case and its
   electrostatic copy, on 2 MPI ranks with the GCC 12.3 / OpenMPI 4.1.6 build
   on Pitagora. Larger parallel layouts and other compilers haven't been
   tested.
10. **There is a small memory overhead.** Each diagnostic call allocates one
    temporary array of size `naky × nakx × nztot × ntubes × nspec × 3`
    (reals), which is negligible compared with the distribution function.
