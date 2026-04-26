# mms-curvature

Python 3 routines for calculating magnetic field curvature, spatial gradients, curl, and divergence
at the mesocenter of the MMS fleet in tetrahedron formation, using data from NASA's Magnetospheric
Multiscale (MMS) mission.

[![DOI](https://zenodo.org/badge/187900473.svg)](https://zenodo.org/badge/latestdoi/187900473)

---

## Version 0.9

---

## Contents

- [Installation](#installation)
- [Quick Start](#quick-start)
- [Core API](#core-api)
  - [mms\_Grad](#mms_grad)
  - [mms\_Curvature](#mms_curvature)
  - [mms\_CurlB](#mms_curlb)
  - [mms\_DivB](#mms_divb)
- [Data Loading Functions](#data-loading-functions)
- [Curvature Runner (utils/Curvature\_Runner\_v6.py)](#curvature-runner)
  - [Command-Line Usage](#command-line-usage)
  - [Runner Functions](#runner-functions)
  - [Output Columns](#output-columns)
- [References](#references)
- [License](#license)

---

## Installation

Clone the repository and install the dependencies:

```bash
git clone https://github.com/unh-mms-rogers/mms_curvature.git
cd mms_curvature
pip install -r requirements.txt
```

**Core dependencies** (`requirements.txt`):

```
numpy
pandas
requests
cdflib
python-dateutil
```

**Optional extras** (`extra_reqs.txt`):

```
matplotlib
pandas-stubs
```

The package is used directly from the repository root (no `pip install` step is required). Ensure
the repository root is on your `PYTHONPATH` or run scripts from within it.

---

## Quick Start

```python
import re
import numpy as np
from mms_curvature import mms_load_fgm, mms_load_mec, mms_Grad, mms_Curvature

trange = ['2017-05-04', '2017-05-05']

# Load FGM magnetic field data
fgmdata, _ = mms_load_fgm(trange=trange, probe=['1','2','3','4'],
                           data_rate='srvy', time_clip=True)

# Extract position and field arrays for each probe
postimes  = [fgmdata['mms'+p+'_fgm_r_gsm_srvy_l2']['x'] for p in '1234']
posvalues = [fgmdata['mms'+p+'_fgm_r_gsm_srvy_l2']['y'] for p in '1234']
magtimes  = [fgmdata['mms'+p+'_fgm_b_gsm_srvy_l2']['x'] for p in '1234']
magvalues = [fgmdata['mms'+p+'_fgm_b_gsm_srvy_l2']['y'] for p in '1234']

# Calculate the spatial gradient and related quantities
grad_Harvey, bm, bmag, rm, t_master = mms_Grad(
    postimes=postimes, posvalues=posvalues,
    magtimes=magtimes, magvalues=magvalues,
    normalize=True, method='RV')

# Calculate the curvature vector k = b · grad(b)
curve_Harvey = mms_Curvature(grad_Harvey, bm)

# Save results
np.savetxt("t_master.csv",      t_master,     delimiter=",")
np.savetxt("curve_Harvey.csv",  curve_Harvey, delimiter=",")
np.save("grad_Harvey.npy", grad_Harvey)
```

---

## Core API

### mms\_Grad

```python
mms_Grad(postimes, posvalues, magtimes, magvalues, normalize=True, method='RV')
```

Dispatcher for the spatial gradient calculation. Selects between the Reciprocal Vectors (`'RV'`)
or Least Squares Minimization (`'LSM'`) methods.

**Parameters**

| Parameter | Type | Description |
|---|---|---|
| `postimes` | list of `np.ndarray` | Unix timestamps for each spacecraft's position data. Order must be consistent (e.g. MMS1 at index 0, MMS4 at index 3). |
| `posvalues` | list of `np.ndarray` | Position vectors for each spacecraft (km, GSM). |
| `magtimes` | list of `np.ndarray` | Unix timestamps for each spacecraft's FGM data. |
| `magvalues` | list of `np.ndarray` | Magnetic field vectors for each spacecraft (nT, GSM). Accepts either `(N, 3)` (Bx, By, Bz) or `(N, 4)` (Bx, By, Bz, \|B\|). |
| `normalize` | bool | If `True` (default), normalizes B before calculating the gradient — required for curvature. If `False`, uses full B vectors — required for curl and divergence. |
| `method` | str | `'RV'` (default) — reciprocal vectors, efficient for 4-spacecraft systems. `'LSM'` — least squares minimization, valid for arbitrary multi-observatory systems. |

**Returns** `(grad_Harvey, bm, bmag, rm, t_master)`

| Variable | Shape | Description |
|---|---|---|
| `grad_Harvey` | `(N, 3, 3)` | Spatial gradient of B (or b = B/\|B\|) in GSM, units of 1/km |
| `bm` | `(N, 3)` | Barycentric mean B unit vector (or full B if `normalize=False`) |
| `bmag` | `(N,)` | Mean \|B\| across spacecraft, in nT |
| `rm` | `(N, 3)` | Mesocenter position in km |
| `t_master` | `(N,)` | Master time series (Unix seconds, float64) |

The RV method additionally accepts `rvecs=True` to append the reciprocal vectors array of shape
`(4, N, 3)` to the return tuple.

---

### mms\_Curvature

```python
mms_Curvature(grad, bm)
```

Calculates the magnetic field line curvature vector **k = b · ∇b**.

**Parameters**

| Parameter | Shape | Description |
|---|---|---|
| `grad` | `(N, 3, 3)` | Gradient of normalized B — output of `mms_Grad(..., normalize=True)` |
| `bm` | `(N, 3)` | Normalized barycentric B vector — second return value of `mms_Grad` |

**Returns** `curve_Harvey` — shape `(N, 3)`, curvature vector in GSM coordinates, units of 1/km.

The radius of curvature in km is `1 / np.linalg.norm(curve_Harvey, axis=1)`.

---

### mms\_CurlB

```python
mms_CurlB(Grad)
```

Calculates ∇ × B from the gradient of the **full** (unnormalized) magnetic field.

**Parameters**

| Parameter | Shape | Description |
|---|---|---|
| `Grad` | `(N, 3, 3)` | Output of `mms_Grad(..., normalize=False)` |

**Returns** `CurlB` — shape `(N, 3)`.

---

### mms\_DivB

```python
mms_DivB(Grad)
```

Calculates ∇ · B as the trace of the gradient of the **full** magnetic field.

**Parameters**

| Parameter | Shape | Description |
|---|---|---|
| `Grad` | `(N, 3, 3)` | Output of `mms_Grad(..., normalize=False)` |

**Returns** `DivB` — shape `(N,)`.

---

## Data Loading Functions

The package exposes shim functions that download and cache MMS data from the Science Data Center.
All return a `(data_dict, metadata_dict)` tuple.

```python
from mms_curvature import (
    mms_load_fgm,        # FGM magnetic field and position
    mms_load_mec,        # MEC ephemeris / attitude
    mms_load_fpi,        # FPI particle moments (dis-moms, des-moms)
    mms_load_hpca,       # HPCA ion moments
    mms_load_scm,        # SCM search-coil magnetometer
    mms_load_edi,        # EDI electron drift instrument
    mms_load_edp,        # EDP electric field
    mms_load_dsp,        # DSP digital signal processor
    mms_load_aspoc,      # ASPOC active spacecraft potential control
    mms_load_fsm,        # FSM fluxgate search-coil merged
    mms_load_ancillary,  # Ancillary products (e.g. DEFERR)
)
```

Common keyword arguments for all loaders:

| Argument | Example | Description |
|---|---|---|
| `trange` | `['2017-05-04', '2017-05-05']` | Start and end time strings |
| `probe` | `['1','2','3','4']` | MMS probe(s) to load |
| `data_rate` | `'srvy'` or `'brst'` | Survey or burst mode |
| `level` | `'l2'` | Data level |
| `time_clip` | `True` | Clip returned data to `trange` |

Data dictionary keys follow the pattern `mms{probe}_{instrument}_{product}_{rate}_{level}`, with
`'x'` holding the time array and `'y'` holding the data array.

---

## Curvature Runner

`utils/Curvature_Runner_v6.py` is a self-contained batch processing script that performs the
complete curvature analysis pipeline including uncertainty quantification, plasma parameters, and
CSV/HDF5 output.

### Command-Line Usage

```bash
# Interactive mode — review and confirm parameters before running
python utils/Curvature_Runner_v6.py

# Non-interactive mode — accept all defaults and run immediately
python utils/Curvature_Runner_v6.py --no-prompt
python utils/Curvature_Runner_v6.py -y
```

Default parameters are set at the top of `main()` and can be changed there or overridden
interactively at runtime:

| Parameter | Default | Description |
|---|---|---|
| `trange` | `['2017-05-01', '2017-06-01']` | Analysis time range |
| `data_rate` | `'srvy'` | FGM/FPI data rate (`'srvy'` or `'brst'`) |
| `prefix` | `"~/Work/v6.2/CurveGSM_"` | Output file path and name prefix |
| `suffix` | `"_v6.2"` | Output filename suffix (before `.csv`) |
| `save_csv` | `True` | Write CSV output |
| `save_h5` | `False` | Also write HDF5 output |

Output files are named automatically from the time range, e.g.:
`CurveGSM_2017-05-01--2017-06-01_v6.2.csv`

For sub-day ranges the hours and minutes are included:
`CurveGSM_2017-05-01_0000--2017-05-01_2359_v6.2.csv`

### Runner Functions

The runner exposes all pipeline stages as importable functions:

#### Data Loading

```python
from utils.Curvature_Runner_v6 import load_fgm_data, load_fpi_data, load_positional_uncertainty

pos_times, b_times, pos_values, b_values, gse_values = load_fgm_data(trange, data_rate)
fpidata, fpirate = load_fpi_data(trange, data_rate)
outRerr = load_positional_uncertainty(trange, num_probes=4, t_master=t_master)
```

| Function | Returns | Description |
|---|---|---|
| `load_fgm_data(trange, data_rate, num_probes=4)` | `pos_times, b_times, pos_values, b_values, gse_values` | FGM position and field in GSM and GSE |
| `load_fpi_data(trange, data_rate, level='l2')` | `fpidata, fpirate` | FPI ion/electron moments for all probes; MMS4 silently dropped on failure |
| `load_positional_uncertainty(trange, num_probes, t_master)` | `outRerr` shape `(probes, N, 4)` | DEFERR ancillary data interpolated to `t_master`, in km |

#### Nominal Products

```python
from utils.Curvature_Runner_v6 import calc_nominal

grad_0n, grad_0f, bm_0, Bmag_0, rm_0, t_master, curve_0, curl_0, div_0 = calc_nominal(
    pos_times, pos_values, b_times, b_values)
```

Returns both normalized (`_0n`) and full-vector (`_0f`) gradients plus curvature, curl, and
divergence in a single call.

#### Uncertainty Quantification

```python
from utils.Curvature_Runner_v6 import (
    calc_positional_uncertainty, calc_magnetometer_uncertainty, combine_uncertainties)

r_unc = calc_positional_uncertainty(
    pos_times, pos_values, b_times, b_values, outRerr, t_master,
    grad_0n, grad_0f, curve_0, curl_0, div_0)

b_unc = calc_magnetometer_uncertainty(
    pos_times, pos_values, b_times, b_values,
    grad_0n, grad_0f, curve_0, curl_0, div_0,
    mag_uncertainty=0.1)   # nT; default ±0.1 nT per component per probe

(sum_uncertainty_grad_n, sum_uncertainty_grad_f,
 sum_uncertainty_curve, sum_uncertainty_curl,
 sum_uncertainty_div, uncertainty_rb_ratio_n) = combine_uncertainties(*r_unc, *b_unc)
```

Uncertainties are propagated by perturbing each input (position component or B component) by its
estimated error and accumulating squared deviations (RSS method). The final combined uncertainty
is `sqrt(sigma_r^2 + sigma_b^2)`.

#### Plasma Parameters

```python
from utils.Curvature_Runner_v6 import mesoGyroradius, calc_plasma_beta, calc_perp_ion_velocity

r_i, r_e = mesoGyroradius(fpidata=fpidata, fpirate=fpirate, t_master=t_master, bmag=Bmag_0)
# r_i, r_e in meters; divide by 1000 for km

beta_i, beta_e, beta_total = calc_plasma_beta(
    fpidata=fpidata, fpirate=fpirate, t_master=t_master, Bmag_0=Bmag_0)

v_perp_mag, v_perp_vec = calc_perp_ion_velocity(
    fpidata=fpidata, fpirate=fpirate, t_master=t_master, bm_0=b_gse)
```

| Function | Returns | Description |
|---|---|---|
| `mesoGyroradius(...)` | `(r_i, r_e)` in meters | Ion and electron thermal gyroradii using FPI T_perp |
| `calc_plasma_beta(...)` | `(beta_i, beta_e, beta_total)` | Plasma beta from FPI number density and T_perp |
| `calc_perp_ion_velocity(...)` | `(v_perp_mag, v_perp_vec)` | Ion bulk velocity component perpendicular to B in km/s |

#### Building and Saving the Output DataFrame

```python
from utils.Curvature_Runner_v6 import build_dataframe, save_results, generate_filename

curvedf = build_dataframe(
    t_master, curve_0, sum_uncertainty_curve, bm_0, Bmag_0, r_i, r_e,
    curl_0, sum_uncertainty_curl, div_0, sum_uncertainty_div, uncertainty_rb_ratio_n,
    beta_i, beta_e, beta_total, v_perp_i)

filename = generate_filename(trange, prefix="./output/CurveGSM_", suffix="_v6.2")
save_results(curvedf, filename, save_csv=True, save_h5=False)
```

### Output Columns

The output DataFrame (index = `Time`, Unix seconds) contains:

| Column | Units | Description |
|---|---|---|
| `Rc(km)` | km | Radius of curvature (1 / \|k\|) |
| `\|curve\|` | 1/km | Magnitude of curvature vector |
| `Curvature_X(GSM)` | 1/km | Curvature vector X component |
| `Curvature_Y(GSM)` | 1/km | Curvature vector Y component |
| `Curvature_Z(GSM)` | 1/km | Curvature vector Z component |
| `error_Kx` | 1/km | Total RSS uncertainty in curvature X |
| `error_Ky` | 1/km | Total RSS uncertainty in curvature Y |
| `error_Kz` | 1/km | Total RSS uncertainty in curvature Z |
| `error_Rc` | km | Propagated uncertainty in radius of curvature |
| `error_\|curve\|` | 1/km | Propagated uncertainty in \|k\| |
| `error_r/b_ratio` | — | Ratio of positional to magnetometer uncertainty norms |
| `b_x` | — | Mean normalized B unit vector X (GSM) |
| `b_y` | — | Mean normalized B unit vector Y (GSM) |
| `b_z` | — | Mean normalized B unit vector Z (GSM) |
| `\|B\|` | nT | Mean field magnitude across the fleet |
| `R_gi(km)` | km | Ion thermal gyroradius |
| `R_ge(km)` | km | Electron thermal gyroradius |
| `curlB_x` | nT/km | Curl of B, X component (GSM) |
| `curlB_y` | nT/km | Curl of B, Y component (GSM) |
| `curlB_z` | nT/km | Curl of B, Z component (GSM) |
| `error_curlx` | nT/km | Total RSS uncertainty in curl X |
| `error_curly` | nT/km | Total RSS uncertainty in curl Y |
| `error_curlz` | nT/km | Total RSS uncertainty in curl Z |
| `div(B)` | nT/km | Divergence of B |
| `error_div(B)` | nT/km | Total RSS uncertainty in div(B) |
| `beta_i` | — | Ion plasma beta |
| `beta_e` | — | Electron plasma beta |
| `beta_total` | — | Total plasma beta |
| `\|v_perp_i\|(km/s)` | km/s | Ion bulk speed perpendicular to B |
| `v_perp_x` | km/s | Perpendicular ion velocity X (GSE) |
| `v_perp_y` | km/s | Perpendicular ion velocity Y (GSE) |
| `v_perp_z` | km/s | Perpendicular ion velocity Z (GSE) |

---

## References

Equations and methods used can be found in:

- *Analysis Methods for Multi-Spacecraft Data* (Paschmann & Daly, Eds.)
  - Ch. 12: "Spatial Gradients and Volumetric Tensor" by Christopher C. Harvey
  - Ch. 14: "Spatial Interpolation for Four Spacecraft: Theory" by G. Chanteur
  - **http://www.issibern.ch/PDF-Files/analysis_methods_1_1a.pdf**

- C. Shen et al. (2003), "Analyses on the geometrical structure of magnetic field in the current
  sheet based on cluster measurements," *Journal of Geophysical Research*.
  **[DOI:10.1029/2002ja009612](https://doi.org/10.1029/2002JA009612)**

---

## License

Released under the Apache 2.0 license (see `LICENSE` file). Click
**[HERE](https://tldrlegal.com/license/apache-license-2.0-(apache-2.0))** for a summary of what
is and is not allowed.

Copyright 2020–2024 Anthony Rogers and Timothy Rogers. All rights reserved.
