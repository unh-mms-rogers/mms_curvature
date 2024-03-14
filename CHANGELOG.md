# Changelog

## [0.8.0] - 04-09-2022

### Added
Harvey curvature vector is calculated.  Returns numpy arrays of:

- *t_master* => time series in unix time of the interpolated times for which the curvature is calculated
- *grad_Harvey* => grad(B/|B|) ordered in GSM coordinates
- *curve_Harvey* => k-vector of magnetic field line curvature in GSM coordinates and units of (1/km)

Radius of curvature calculation will be added shortly for greater application in the same way as the IDL version. 