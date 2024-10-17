# mms-curvature
Python3 routines for calculating the magnetic field curvature at the mesocenter of the MMS fleet in tetrahedron formation.

Also contains basic functions for loading various instrument data from MMS spacecraft, calculating spacial gradients over time, and deriving the curl and divergence from such gradients.

## Version 0.9

**NIX ALL OF THE BLOW AND UPDATE WITH ACTUALL USAGE GUIDE**

Harvey curvature vector is calculated.  Returns numpy arrays of:

- *t_master* => time series in unix time of the interpolated times for which the curvature is calculated
- *grad_Harvey* => grad(B/|B|) ordered in GSM coordinates
- *curve_Harvey* => k-vector of magnetic field line curvature in GSM coordinates and units of (1/km)

Radius of curvature calculation will be added shortly for greater application in the same way as the IDL version.  

NOTE: the python version DOES NOT contain the time step data in the same array as the gradient or curvature vector as in a tplot variable.  For this reason the *t_master* product is returned.

-----

### References
Equations and methods used can be found at:

- *Analysis Methods for Multi-Spacecraft Data* (Paschmann & Daly, Eds.) Ch.12 'Spatial Gradiants and Volumetric Tensor' by Christopher C. Harvey.  **http://www.issibern.ch/PDF-Files/analysis_methods_1_1a.pdf**
- "Analyses on the geometrical structure of magnetic field in the current sheet based on cluster measurements" (2003) by C. Shen, et al.  **[DOI:10.1029/2002ja009612](https://doi.org/10.1029/2002JA009612)**

-----
### License
Released under Apache 2.0 license (see LICENSE file).  Click **[HERE](https://tldrlegal.com/license/apache-license-2.0-(apache-2.0))** for a summary of what is/isn't allowed.



[![DOI](https://zenodo.org/badge/187900473.svg)](https://zenodo.org/badge/latestdoi/187900473)


