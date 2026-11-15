# COMSOL models

COMSOL Multiphysics 6.4 files. They need a COMSOL licence to open; the figures and `data/` are the readable record of the results. Solutions were cleared before saving, so each study must be re-run to regenerate results.

| File | Size | Content |
|---|---|---|
| `DE_membrane_axisym.mph` | 26.8 MB | Axisymmetric membrane model, studies 1-7 (see `docs/comsol_validation.md`): checks against closed forms, deflection-controlled pressure continuation, centre control, dynamic attempt |
| `DE_shell_axisym_Model_VHB.mph` | 2.0 MB | Axisymmetric shell model with the VHB 4910 constants: spectrum of the unstable flat state, pressure continuation |

## Notes

* **Thickness setting.** The membrane results quoted as "frozen thickness" use an External Stress with the stretched thickness `hf = 0.5[mm]/lam^2`. The "deformed thickness" runs use `mbrn.ddef`. The setting saved in the file is the one used last and has not been checked against each result, so confirm it before reproducing a result.
* **Control parameters** in the membrane file: `dh` (deflection target), `sw` (switch between deflection and pressure control), `pc` (pressure value in pressure control), as described in `docs/comsol_validation.md`.
* **Size.** GitHub rejects files above 100 MB per push and warns above 50 MB. Use *File > Compact History* and *Remove Solution* before saving, or Git LFS for larger files.
