# Data

All files are tables copied from COMSOL Multiphysics 6.4 (Derived Values > Global Evaluation, Line Average, Point Evaluation) with no post-processing, except where stated.

| File | Content | Model, setting |
|---|---|---|
| `membrane_frozen_mean_control.csv` | `dh_mm` mean deflection (deformed-area mean, negative = downward), `pu_Pa` holding pressure | membrane, $V_0=9000$ V, frozen thickness. Rows to $-10$ mm (step 0.05) from an earlier run on a coarse mesh (probably 20 elements); rows $-10.1$ to $-11.754$ mm (step 0.1) from a 100-element run. The two runs agree to the displayed digits where they overlap (51.269 Pa at -2.2 mm in both). |
| `membrane_frozen_centre_control.csv` | `dh_centre_mm` centre deflection, `pc_Pa` (0, pressure control off), `mean_deflection_mm`, `pu_Pa` | membrane, 9000 V, frozen thickness, 100 elements, centre-deflection control |
| `membrane_deformed_thickness_mean_control.csv` | `dh_mm`, `pu_Pa`, `line_average_w_mm` | membrane, 9000 V, deformed thickness; solver stops at $-2.2$ mm. The $-2.2$ mm value was 46.816 Pa on 20 elements and 46.830 Pa on 100 elements. |
| `membrane_frozen_pressure_control_tail.csv` | pressure-controlled continuation from 71.5 Pa down, `mean_deflection_mm` | membrane, 9000 V, frozen thickness |
| `shell_pressure_vs_mean_deflection_V8755.csv` | mean deflection and pressure | shell, $V_0=8755$ V |
| `shell_unstable_growth_rates_V8755.csv` | mode number, FE growth rate, membrane-theory rate, difference | shell, flat state, $V_0=8755$ V |
| `membrane_checks_summary.csv` | table of checks against closed forms | membrane (and one shell row) |

The sign convention follows COMSOL: positive pressure pushes towards +z, the membrane deflects opposite to the load (negative $w$), because the flat state has negative stiffness for $V_0>V_c$.
