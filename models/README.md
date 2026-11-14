# COMSOL models

Place the COMSOL Multiphysics 6.4 files here (they need a COMSOL licence to open; the figures and `data/` are the readable record of the results).

| File | Content |
|---|---|
| `DE_membrane_axisym.mph` | membrane model, studies 1-7 (see `docs/comsol_validation.md`) |
| `DE_shell_axisym_Model_VHB.mph` | shell model with the VHB 4910 constants (spectrum of the flat state, pressure continuation) |

**Before you upload:**
* Files above 50 MB trigger a GitHub warning and above 100 MB are rejected; the files used here are about 18 MB. If a file is larger, use Git LFS or *File > Compact History and Save* in COMSOL.
* Save with *File > Compact History*, or *Remove Solution*, to keep the file small and to avoid shipping a stale solution.
* The membrane file must have the External Stress set to the version you describe (`hf` for the frozen thickness, `mbrn.ddef` for the deformed thickness) and the parameters `dh`, `sw`, `pc` as in `docs/comsol_validation.md`.
