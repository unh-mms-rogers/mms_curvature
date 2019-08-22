# mms-curvature
IDL and python3 routines for calculating the magnetic field curvature at the mesocenter of the MMS fleet in tetrahedron formation.


## Version 0.1
### IDL:
Currently have Harvey method working well.  Shen method will follow as soon as I can get to it.  Currently returns tplot variables of:

- *b_gradient_gsm* => grad(B/|B|) ordered in GSM coordinates
- *b_curvature_vector_gsm* => k-vector of magnetic field line curvature in GSM coordinates and units of (1/km)
- *b_curvature_radius* => radius of curvature of magnetic field lines in units of (km)

A more complete user guide will be forthcoming.

### Python3:
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
