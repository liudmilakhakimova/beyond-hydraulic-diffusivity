# Input data

Input datasets for the reaction-front diffusivity calculations of [Khakimova, Schmalholz and Podladchikov (2026)](https://doi.org/10.5194/egusphere-2026-3107).

| File | Contents | Used by |
| --- | --- | --- |
| `DATA_Fusseis_plots.mat` | Digitized porosity profiles: `x1,y1`, `x2,y2`, `x3,y3` | `FIG_8.m` |
| `DATA_Liudmilalike_Dehyd.mat` | Saved single-front dehydration calculation and analytical front positions | `FIG_2.m` |
| `DATA_Liudmilalike_Hyd.mat` | Saved single-front hydration calculation and analytical front positions | `FIG_2.m` |
| `DATA_DEHY_Dehyd.mat` | Saved single-front dehydration calculation | `FIG_2.m` |
| `DATA_DEHY_Hyd.mat` | Saved single-front hydration calculation | `FIG_2.m` |

For the Fusseis profiles, distance is in micrometres and porosity is in percent. The pairs correspond to 600, 4200 and 7800 s.

The porosity profiles are redrawn after:

Fusseis, F., Schrank, C., Liu, J., Karrech, A., Llana-Fúnez, S., Xiao, X., and Regenauer-Lieb, K. (2012). Pore formation during dehydration of a polycrystalline gypsum sample observed and quantified in a time-series synchrotron X-ray micro-tomography experiment. *Solid Earth*, 3, 71-86. [doi:10.5194/se-3-71-2012](https://doi.org/10.5194/se-3-71-2012).

The four simulation datasets are MATLAB workspaces loaded by `FIG_2.m`. The dataset-generating solvers and the digitization procedure are not included in this repository.

## Variables in each file

`DATA_Fusseis_plots.mat` contains three column-vector pairs: `x1,y1` (14 points, 600 s), `x2,y2` (25 points, 4200 s), and `x3,y3` (26 points, 7800 s). These are digitized one-dimensional profiles. Figure 8 plots the third profile as `y3-1`; this offset is applied by the script, not stored in the data.

`DATA_Liudmilalike_Dehyd.mat` contains the final pressure and total-density profiles (`Pf`, `rhot`) on the 200-point grid `x`, a four-point constitutive relation (`rhot_lut`, `P_lut`), and 2,310,059 samples of `Time_vec`, `X_front` and `X_front_ana`. The scalar parameters include `rhof`, `k_etaf`, `dP`, `Delta_rhoT`, `dx` and `dt`. The analytical history uses the saved reaction density jump `Delta_rhoT`.

`DATA_Liudmilalike_Hyd.mat` contains the corresponding hydration calculation on a 200-point grid, with the same field names and 1,155,030 time/front samples. Its constitutive relation and boundary values differ from the dehydration case.

`DATA_DEHY_Dehyd.mat` contains final nodal profiles on the 151-point grid `x`, including `Pf`, `rhot` and `phi`; 150 values of Darcy flux `qD`; and 240 samples of `Time_vec` and `X_front`. It also stores the pressure/density and material-property lookup arrays, solver settings and residuals. No `X_front_ana` array is present.

`DATA_DEHY_Hyd.mat` contains the corresponding hydration calculation on a 151-point grid, with 150 flux values and 253 time/front samples. Its saved fields include the same pressure, density, porosity, lookup-table and solver quantities. No `X_front_ana` array is present.

`FIG_2.m` loads the four simulation workspaces in sequence. The plotted quantities use the normalization stored in these files, without unit conversion.
