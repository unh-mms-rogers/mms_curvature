# Copyright 2020-2022 Anthony Rogers.  All rights reserved.
# Released under the Apache 2.0 license.

import copy
import time
import numpy as np
import pandas as pd
from mms_curvature.mms_curvature import mms_Grad, mms_Curvature, mms_CurlB, mms_DivB
from mms_curvature.mms.mms_load_data_shims import mms_load_fgm, mms_load_fpi, mms_load_ancillary


def mesoGyroradius(fpidata, fpirate, t_master=None, bmag=None):
    '''
    Calculates the average gyroradius in the tetrahedron, assumed to be
    representative of the gyroradius at the mesocenter.  Uses FPI data.

    inputs:
    fpidata:    dict of pre-loaded FPI moments data as returned by load_fpi_data
    fpirate:    FPI data rate string ('fast' or 'brst'), as returned by load_fpi_data
    t_master:   Time series for the magnetic field data, assumed from
                mms_curvature.Curvature
    bmag:       |B| aligned with t_master.  Assumed from mms_curvature.Curvature

    outputs:
    (r_i, r_e) -- 2-element structure aligned with t_master time series where:
        r_i:        average gyroradius of ions (assumed protons) in meters
        r_e:        average gyroradius of electrons in meters
    '''

    me = 9.1094e-31     # electron mass in kilograms
    mp = 1.6726e-27     # proton mass in kg; used as proxy for all ions
    q  = 1.602177e-19   # elementary charge in Coulombs

    # Extract time and T_perp for probes 1-3 (always present)
    ion_times  = []
    ion_tperp  = []
    elec_times = []
    elec_tperp = []
    for probe in ['1', '2', '3']:
        ion_times.append(fpidata['mms' + probe + '_dis_tempperp_' + fpirate]['x'])
        ion_tperp.append(fpidata['mms' + probe + '_dis_tempperp_' + fpirate]['y'])
        elec_times.append(fpidata['mms' + probe + '_des_tempperp_' + fpirate]['x'])
        elec_tperp.append(fpidata['mms' + probe + '_des_tempperp_' + fpirate]['y'])

    mitime = ion_times[0]   # master ion time (mms1)
    metime = elec_times[0]  # master electron time (mms1)

    # Interpolate probes 2 and 3 onto mms1 time grid
    for i in range(1, 3):
        ion_tperp[i]  = np.interp(mitime, ion_times[i],  ion_tperp[i])
        elec_tperp[i] = np.interp(metime, elec_times[i], elec_tperp[i])

    # Include mms4 if it was successfully loaded
    if 'mms4_dis_tempperp_' + fpirate in fpidata:
        ion_tperp4  = np.interp(mitime, fpidata['mms4_dis_tempperp_' + fpirate]['x'],
                                         fpidata['mms4_dis_tempperp_' + fpirate]['y'])
        elec_tperp4 = np.interp(metime, fpidata['mms4_des_tempperp_' + fpirate]['x'],
                                          fpidata['mms4_des_tempperp_' + fpirate]['y'])
        tiperp = np.average([ion_tperp[0],  ion_tperp[1],  ion_tperp[2],  ion_tperp4],  axis=0)
        teperp = np.average([elec_tperp[0], elec_tperp[1], elec_tperp[2], elec_tperp4], axis=0)
    else:
        tiperp = np.average([ion_tperp[0],  ion_tperp[1],  ion_tperp[2]],  axis=0)
        teperp = np.average([elec_tperp[0], elec_tperp[1], elec_tperp[2]], axis=0)

    # Interpolate averaged T_perp to t_master and calculate gyroradii.
    # Assumes tperp in eV and bmag in nT; result is in meters.
    tiperp_m = np.interp(t_master, mitime, tiperp)
    teperp_m = np.interp(t_master, metime, teperp)

    r_i = np.sqrt(2 * mp * 11605 * 1.38649e-23) / (1e-9 * q) * np.sqrt(tiperp_m) / bmag
    r_e = np.sqrt(2 * me * 11605 * 1.38649e-23) / (1e-9 * q) * np.sqrt(teperp_m) / bmag

    return (r_i, r_e)


def calc_plasma_beta(fpidata, fpirate, t_master, Bmag_0):
    '''
    Calculate ion, electron, and total plasma beta averaged across the
    MMS tetrahedron, using pre-loaded FPI moments and FGM field magnitude.

    Beta is the ratio of plasma thermal pressure to magnetic pressure:
        beta = (n * k_B * T) / (B^2 / 2*mu_0)

    Temperatures used are T_perp from FPI moments (eV).  For a more
    rigorous result the full scalar temperature should be substituted.
    Density and temperature are averaged across all available probes;
    mms4 is included automatically if present in fpidata.

    inputs:
    fpidata:    dict of pre-loaded FPI moments data as returned by load_fpi_data
    fpirate:    FPI data rate string ('fast' or 'brst'), as returned by load_fpi_data
    t_master:   master time series aligned with Bmag_0 (from calc_nominal)
    Bmag_0:     mean |B| in nT aligned with t_master (from calc_nominal)

    outputs:
    (beta_i, beta_e, beta_total) -- arrays aligned with t_master
    '''
    mu0 = 4 * np.pi * 1e-7   # permeability of free space (H/m)
    q   = 1.602177e-19        # elementary charge (C)

    # Collect density and T_perp for probes 1-3
    ion_n_t  = [];  ion_n  = []
    ion_tp_t = [];  ion_tp = []
    elec_n_t = [];  elec_n = []
    elec_tp_t= [];  elec_tp= []

    for p in ['1', '2', '3']:
        pref = 'mms' + p + '_'
        ion_n_t.append(  fpidata[pref + 'dis_numberdensity_' + fpirate]['x'])
        ion_n.append(    fpidata[pref + 'dis_numberdensity_' + fpirate]['y'])
        ion_tp_t.append( fpidata[pref + 'dis_tempperp_'      + fpirate]['x'])
        ion_tp.append(   fpidata[pref + 'dis_tempperp_'      + fpirate]['y'])
        elec_n_t.append( fpidata[pref + 'des_numberdensity_' + fpirate]['x'])
        elec_n.append(   fpidata[pref + 'des_numberdensity_' + fpirate]['y'])
        elec_tp_t.append(fpidata[pref + 'des_tempperp_'      + fpirate]['x'])
        elec_tp.append(  fpidata[pref + 'des_tempperp_'      + fpirate]['y'])

    # mms1 as master time for each quantity
    mi_n_t  = ion_n_t[0];   mi_tp_t = ion_tp_t[0]
    me_n_t  = elec_n_t[0];  me_tp_t = elec_tp_t[0]

    # Interpolate probes 2 and 3 onto mms1 time grid
    for i in range(1, 3):
        ion_n[i]  = np.interp(mi_n_t,  ion_n_t[i],   ion_n[i])
        ion_tp[i] = np.interp(mi_tp_t, ion_tp_t[i],  ion_tp[i])
        elec_n[i] = np.interp(me_n_t,  elec_n_t[i],  elec_n[i])
        elec_tp[i]= np.interp(me_tp_t, elec_tp_t[i], elec_tp[i])

    # Include mms4 if it was successfully loaded
    if 'mms4_dis_numberdensity_' + fpirate in fpidata:
        ion_n.append( np.interp(mi_n_t,  fpidata['mms4_dis_numberdensity_' + fpirate]['x'],
                                          fpidata['mms4_dis_numberdensity_' + fpirate]['y']))
        ion_tp.append(np.interp(mi_tp_t, fpidata['mms4_dis_tempperp_'      + fpirate]['x'],
                                          fpidata['mms4_dis_tempperp_'      + fpirate]['y']))
        elec_n.append( np.interp(me_n_t,  fpidata['mms4_des_numberdensity_' + fpirate]['x'],
                                           fpidata['mms4_des_numberdensity_' + fpirate]['y']))
        elec_tp.append(np.interp(me_tp_t, fpidata['mms4_des_tempperp_'      + fpirate]['x'],
                                           fpidata['mms4_des_tempperp_'      + fpirate]['y']))

    # Average across all available probes
    ni_avg = np.average(ion_n,  axis=0)   # cm^-3
    ti_avg = np.average(ion_tp, axis=0)   # eV
    ne_avg = np.average(elec_n, axis=0)   # cm^-3
    te_avg = np.average(elec_tp,axis=0)   # eV

    # Interpolate averaged quantities to t_master
    ni = np.interp(t_master, mi_n_t,  ni_avg)
    ti = np.interp(t_master, mi_tp_t, ti_avg)
    ne = np.interp(t_master, me_n_t,  ne_avg)
    te = np.interp(t_master, me_tp_t, te_avg)

    # Magnetic pressure (Pa): B in nT -> SI
    P_mag = (Bmag_0 * 1e-9)**2 / (2 * mu0)

    # Plasma pressure (Pa): n in cm^-3 -> SI, T in eV -> J
    P_ion  = (ni * 1e6) * (ti * q)
    P_elec = (ne * 1e6) * (te * q)

    beta_i     = P_ion  / P_mag
    beta_e     = P_elec / P_mag
    beta_total = (P_ion + P_elec) / P_mag

    return beta_i, beta_e, beta_total


def calc_perp_ion_velocity(fpidata, fpirate, t_master, bm_0):
    '''
    Calculate the magnitude of the ion bulk velocity component perpendicular
    to the magnetic field, averaged across the MMS tetrahedron.

    The perpendicular velocity is computed as:
        v_perp = v - (v . b_hat) * b_hat
    where b_hat is the mean normalized magnetic field unit vector (bm_0).

    inputs:
    fpidata:    dict of pre-loaded FPI moments data as returned by load_fpi_data
    fpirate:    FPI data rate string ('fast' or 'brst'), as returned by load_fpi_data
    t_master:   master time series aligned with bm_0 (from calc_nominal)
    bm_0:       mean normalized B unit vector array (n, 3) from calc_nominal

    outputs:
    v_perp_mag: magnitude of perpendicular ion bulk velocity in km/s,
                aligned with t_master
    '''

    # Collect bulk velocity vectors for probes 1-3
    v_times = []
    v_vals  = []
    for p in ['1', '2', '3']:
        pref = 'mms' + p + '_'
        v_times.append(fpidata[pref + 'dis_bulkv_' + fpirate]['x'])
        v_vals.append( fpidata[pref + 'dis_bulkv_' + fpirate]['y'])

    mitime = v_times[0]

    # Interpolate probes 2 and 3 onto mms1 time grid (component-wise)
    for i in range(1, 3):
        v_interp = np.zeros_like(v_vals[0])
        for dim in range(3):
            v_interp[:, dim] = np.interp(mitime, v_times[i], v_vals[i][:, dim])
        v_vals[i] = v_interp

    # Include mms4 if it was successfully loaded
    if 'mms4_dis_bulkv_' + fpirate in fpidata:
        v4 = fpidata['mms4_dis_bulkv_' + fpirate]
        v_interp4 = np.zeros_like(v_vals[0])
        for dim in range(3):
            v_interp4[:, dim] = np.interp(mitime, v4['x'], v4['y'][:, dim])
        v_vals.append(v_interp4)

    # Average across all available probes
    v_avg = np.average(v_vals, axis=0)   # shape (n_fpi, 3), km/s

    # Interpolate averaged velocity to t_master (component-wise)
    v_master = np.zeros((len(t_master), 3))
    for dim in range(3):
        v_master[:, dim] = np.interp(t_master, mitime, v_avg[:, dim])

    # Remove parallel component: v_perp = v - (v . b_hat) * b_hat
    v_dot_b    = np.sum(v_master * bm_0, axis=1)     # scalar projection onto B
    v_para     = v_dot_b[:, np.newaxis] * bm_0        # parallel component vector
    v_perp     = v_master - v_para                    # perpendicular component vector
    v_perp_mag = np.linalg.norm(v_perp, axis=1)       # magnitude in km/s

    return v_perp_mag


def generate_filename(trange, prefix, suffix):
    '''
    Generate a sane output filename from trange, prefix, and suffix.

    inputs:
    trange:     2-element list of start/end time strings
    prefix:     path and filename prefix string
    suffix:     filename suffix string (before .csv extension)

    outputs:
    filename:   full output filename string (always ends in .csv)
    '''
    if len(trange[0]) > 10 or len(trange[1]) > 10:
        return (prefix
                + trange[0][:10] + "_" + trange[0][11:13] + trange[0][14:16]
                + "--"
                + trange[1][:10] + "_" + trange[1][11:13] + trange[1][14:16]
                + suffix + ".csv")
    return prefix + trange[0][:10] + "--" + trange[1][:10] + suffix + ".csv"


def load_fgm_data(trange, data_rate, num_probes=4):
    '''
    Load FGM position and magnetic field data for all probes.

    inputs:
    trange:     2-element list of start/end time strings
    data_rate:  'srvy' or 'brst'
    num_probes: number of MMS probes to load (default 4)

    outputs:
    (pos_times, b_times, pos_values, b_values) -- each a list of length
    num_probes containing the time and value arrays for each probe
    '''
    fgmdata = mms_load_fgm(trange=trange, probe=['1', '2', '3', '4'], data_rate=data_rate, time_clip=True)[0]
    pos_times  = [None] * num_probes
    b_times    = [None] * num_probes
    pos_values = [None] * num_probes
    b_values   = [None] * num_probes
    for probe in range(num_probes):
        key = 'mms' + str(probe + 1) + '_fgm_'
        pos_times[probe]  = np.copy(fgmdata[key + 'r_gsm_' + data_rate + '_l2']['x'])
        b_times[probe]    = np.copy(fgmdata[key + 'b_gsm_' + data_rate + '_l2']['x'])
        pos_values[probe] = np.copy(fgmdata[key + 'r_gsm_' + data_rate + '_l2']['y'])
        b_values[probe]   = np.copy(fgmdata[key + 'b_gsm_' + data_rate + '_l2']['y'])
    return pos_times, b_times, pos_values, b_values


def load_fpi_data(trange, data_rate, level='l2'):
    '''
    Load FPI moments (dis-moms, des-moms) for all probes.  Probes 1-3 are
    always loaded; probe 4 is attempted and silently dropped on failure.

    inputs:
    trange:     2-element list of start/end time strings
    data_rate:  'srvy' or 'brst'
    level:      data level (default 'l2')

    outputs:
    (fpidata, fpirate) where:
        fpidata:    merged dict of all successfully loaded FPI variables
        fpirate:    FPI data rate string used ('fast' or 'brst')
    '''
    fpirate = 'fast' if data_rate == 'srvy' else data_rate
    fpidata = {}
    for probe in ['1', '2', '3']:
        data, _ = mms_load_fpi(trange=trange, probe=probe, data_rate=fpirate, level=level,
                               datatype=['dis-moms', 'des-moms'], time_clip=True)
        fpidata.update(data)
    try:
        data4, _ = mms_load_fpi(trange=trange, probe='4', data_rate=fpirate, level=level,
                                datatype=['dis-moms', 'des-moms'], time_clip=True)
        fpidata.update(data4)
    except:
        print('Error loading mms4 FPI data.  Will drop mms4 from this dataset.')
    return fpidata, fpirate


def load_positional_uncertainty(trange, num_probes, t_master):
    '''
    Load DEFERR ancillary positional uncertainty data and interpolate onto
    t_master -- the authoritative time grid produced by calc_nominal --
    following the same component-wise np.interp approach used by mms_Grad_RV
    for magnetic field data.  Each probe's DEFERR array is interpolated
    directly without first converting the per-probe list to a numpy array,
    avoiding inhomogeneous-shape failures when probes have different numbers
    of DEFERR time points.

    inputs:
    trange:     2-element list of start/end time strings
    num_probes: number of MMS probes
    t_master:   master time array as returned by calc_nominal

    outputs:
    outRerr:    array of shape (num_probes, len(t_master), 4) with positional
                uncertainty in kilometers, aligned with t_master
    '''
    deferr_in = mms_load_ancillary(probe=['1', '2', '3', '4'], anc_product='deferr', trange=trange, time_clip=True)
    Rerr_arr = [None] * num_probes
    for probe in range(1, num_probes + 1):
        Rerr_arr[probe - 1] = deferr_in[0]["MMS" + str(probe) + "_DEFERR"].to_numpy()

    # Interpolate each probe's DEFERR component-wise directly onto t_master.
    # Rerr_arr[bird] is kept as a list of per-probe arrays to avoid the
    # inhomogeneous-shape issue that arises from np.asarray() when probe
    # DEFERR arrays have differing lengths.
    outRerr = np.ndarray((num_probes, len(t_master), 4))
    for bird in range(num_probes):
        for dim in range(4):
            outRerr[bird, :, dim] = np.interp(t_master, Rerr_arr[bird][:, 0], Rerr_arr[bird][:, dim])

    # Convert positional uncertainty from meters to kilometers
    return 1e-3 * outRerr


def calc_nominal(pos_times, pos_values, b_times, b_values):
    '''
    Calculate nominal (zero-uncertainty) gradient, curvature, curl, and
    divergence products.

    inputs:
    pos_times:  list of position time arrays (one per probe)
    pos_values: list of position value arrays (one per probe)
    b_times:    list of magnetic field time arrays (one per probe)
    b_values:   list of magnetic field value arrays (one per probe)

    outputs:
    (grad_0n, grad_0f, bm_0, Bmag_0, rm_0, t_master, curve_0, curl_0, div_0)
        grad_0n:    gradient of normalized B
        grad_0f:    gradient of full (unnormalized) B
        bm_0:       mean normalized B vector
        Bmag_0:     mean |B|
        rm_0:       mesocenter position
        t_master:   master time series
        curve_0:    curvature vector
        curl_0:     curl of B
        div_0:      divergence of B
    '''
    grad_0n, bm_0, Bmag_0, rm_0, t_master = mms_Grad(
        postimes=pos_times, posvalues=pos_values, magtimes=b_times, magvalues=b_values, normalize=True)
    grad_0f, Bm_0 = mms_Grad(
        postimes=pos_times, posvalues=pos_values, magtimes=b_times, magvalues=b_values, normalize=False)[:2]
    curve_0 = mms_Curvature(grad_0n, bm_0)
    curl_0  = mms_CurlB(grad_0f)
    div_0   = mms_DivB(grad_0f)
    return grad_0n, grad_0f, bm_0, Bmag_0, rm_0, t_master, curve_0, curl_0, div_0


def calc_positional_uncertainty(pos_times, pos_values, b_times, b_values, outRerr,
                                t_master, grad_0n, grad_0f, curve_0, curl_0, div_0):
    '''
    Compute squared uncertainty contributions from positional errors (DEFERR).
    Each probe's position is perturbed by +/- its uncertainty in each spatial
    dimension and the resulting change in each product is accumulated.

    outRerr is assumed to be aligned with t_master (as returned by
    load_positional_uncertainty).  pos_values are resampled to t_master
    before perturbations are applied so that all arrays share a common grid.

    inputs:
    pos_times:  list of position time arrays (one per probe)
    pos_values: list of nominal position value arrays (one per probe)
    b_times:    list of magnetic field time arrays (one per probe)
    b_values:   list of magnetic field value arrays (one per probe)
    outRerr:    positional uncertainty array (num_probes, len(t_master), 4) in km
    t_master:   master time array from calc_nominal
    grad_0n:    nominal gradient of normalized B
    grad_0f:    nominal gradient of full B
    curve_0:    nominal curvature vector
    curl_0:     nominal curl of B
    div_0:      nominal divergence of B

    outputs:
    (r_uncertainty_grad_n, r_uncertainty_curve, r_uncertainty_grad_f,
     r_uncertainty_curl, r_uncertainty_div) -- summed squared deviations for
     each product due to positional uncertainty
    '''
    num_probes = len(pos_values)
    r_uncertainty_grad_n = np.zeros_like(grad_0n)
    r_uncertainty_curve  = np.zeros_like(curve_0)
    r_uncertainty_grad_f = np.zeros_like(grad_0f)
    r_uncertainty_curl   = np.zeros_like(curl_0)
    r_uncertainty_div    = np.zeros_like(div_0)

    # Resample each probe's position to t_master so it is consistent with outRerr
    master_pos_times = [t_master] * num_probes
    pos_on_master = [np.zeros((len(t_master), 3)) for _ in range(num_probes)]
    for probe in range(num_probes):
        for dim in range(3):
            pos_on_master[probe][:, dim] = np.interp(t_master, pos_times[probe], pos_values[probe][:, dim])

    tpos = copy.deepcopy(pos_on_master)
    for probe in range(num_probes):
        for spatial_dim in range(3):
            for sign in [-1, 1]:
                tpos[probe][:, spatial_dim] = np.add(
                    pos_on_master[probe][:, spatial_dim],
                    np.multiply(outRerr[probe, :, spatial_dim + 1], sign))

                grad_in, bm_i = mms_Grad(postimes=master_pos_times, posvalues=tpos, magtimes=b_times, magvalues=b_values, normalize=True)[:2]
                grad_if = mms_Grad(postimes=master_pos_times, posvalues=tpos, magtimes=b_times, magvalues=b_values, normalize=False)[0]
                curve_i = mms_Curvature(grad_in, bm_i)
                curl_i  = mms_CurlB(grad_if)
                div_i   = mms_DivB(grad_if)

                r_uncertainty_grad_n = np.add(np.power(np.subtract(grad_in, grad_0n), 2), r_uncertainty_grad_n)
                r_uncertainty_curve  = np.add(np.power(np.subtract(curve_i, curve_0), 2), r_uncertainty_curve)
                r_uncertainty_grad_f = np.add(np.power(np.subtract(grad_if, grad_0f), 2), r_uncertainty_grad_f)
                r_uncertainty_curl   = np.add(np.power(np.subtract(curl_i,  curl_0),  2), r_uncertainty_curl)
                r_uncertainty_div    = np.add(np.power(np.subtract(div_i,   div_0),   2), r_uncertainty_div)

                tpos = copy.deepcopy(pos_on_master)

    return r_uncertainty_grad_n, r_uncertainty_curve, r_uncertainty_grad_f, r_uncertainty_curl, r_uncertainty_div


def calc_magnetometer_uncertainty(pos_times, pos_values, b_times, b_values,
                                  grad_0n, grad_0f, curve_0, curl_0, div_0,
                                  mag_uncertainty=0.1):
    '''
    Compute squared uncertainty contributions from magnetometer measurement
    errors.  Each probe's B field is perturbed by +/- mag_uncertainty in each
    spatial dimension and the resulting change in each product is accumulated.

    inputs:
    pos_times:       list of position time arrays (one per probe)
    pos_values:      list of position value arrays (one per probe)
    b_times:         list of magnetic field time arrays (one per probe)
    b_values:        list of nominal magnetic field value arrays (one per probe)
    grad_0n:         nominal gradient of normalized B
    grad_0f:         nominal gradient of full B
    curve_0:         nominal curvature vector
    curl_0:          nominal curl of B
    div_0:           nominal divergence of B
    mag_uncertainty: magnetometer measurement uncertainty in nT (default 0.1)

    outputs:
    (b_uncertainty_grad_n, b_uncertainty_curve, b_uncertainty_grad_f,
     b_uncertainty_curl, b_uncertainty_div) -- summed squared deviations for
     each product due to magnetometer uncertainty
    '''
    num_probes = len(b_values)
    b_uncertainty_grad_n = np.zeros_like(grad_0n)
    b_uncertainty_curve  = np.zeros_like(curve_0)
    b_uncertainty_grad_f = np.zeros_like(grad_0f)
    b_uncertainty_curl   = np.zeros_like(curl_0)
    b_uncertainty_div    = np.zeros_like(div_0)

    tb = copy.deepcopy(b_values)
    for probe in range(num_probes):
        for spatial_dim in range(3):
            for delta in [-mag_uncertainty, mag_uncertainty]:
                tb[probe][:, spatial_dim] = np.add(b_values[probe][:, spatial_dim], delta)

                grad_in, bm_i = mms_Grad(postimes=pos_times, posvalues=pos_values, magtimes=b_times, magvalues=tb, normalize=True)[:2]
                grad_if = mms_Grad(postimes=pos_times, posvalues=pos_values, magtimes=b_times, magvalues=tb, normalize=False)[0]
                curve_i = mms_Curvature(grad_in, bm_i)
                curl_i  = mms_CurlB(grad_if)
                div_i   = mms_DivB(grad_if)

                b_uncertainty_grad_n = np.add(np.power(np.subtract(grad_in, grad_0n), 2), b_uncertainty_grad_n)
                b_uncertainty_curve  = np.add(np.power(np.subtract(curve_i, curve_0), 2), b_uncertainty_curve)
                b_uncertainty_grad_f = np.add(np.power(np.subtract(grad_if, grad_0f), 2), b_uncertainty_grad_f)
                b_uncertainty_curl   = np.add(np.power(np.subtract(curl_i,  curl_0),  2), b_uncertainty_curl)
                b_uncertainty_div    = np.add(np.power(np.subtract(div_i,   div_0),   2), b_uncertainty_div)

                tb = copy.deepcopy(b_values)

    return b_uncertainty_grad_n, b_uncertainty_curve, b_uncertainty_grad_f, b_uncertainty_curl, b_uncertainty_div


def combine_uncertainties(r_uncertainty_grad_n, r_uncertainty_curve, r_uncertainty_grad_f,
                          r_uncertainty_curl, r_uncertainty_div,
                          b_uncertainty_grad_n, b_uncertainty_curve, b_uncertainty_grad_f,
                          b_uncertainty_curl, b_uncertainty_div):
    '''
    Combine positional and magnetometer squared uncertainties into total
    RSS uncertainty estimates for each product.

    inputs:
    r_uncertainty_*: squared positional uncertainty arrays from calc_positional_uncertainty
    b_uncertainty_*: squared magnetometer uncertainty arrays from calc_magnetometer_uncertainty

    outputs:
    (sum_uncertainty_grad_n, sum_uncertainty_grad_f, sum_uncertainty_curve,
     sum_uncertainty_curl, sum_uncertainty_div, uncertainty_rb_ratio_n)
    '''
    sum_uncertainty_grad_n = np.sqrt(np.add(r_uncertainty_grad_n, b_uncertainty_grad_n))
    sum_uncertainty_grad_f = np.sqrt(np.add(r_uncertainty_grad_f, b_uncertainty_grad_f))
    sum_uncertainty_curve  = np.sqrt(np.add(r_uncertainty_curve,  b_uncertainty_curve))
    sum_uncertainty_curl   = np.sqrt(np.add(r_uncertainty_curl,   b_uncertainty_curl))
    sum_uncertainty_div    = np.sqrt(np.add(r_uncertainty_div,    b_uncertainty_div))
    uncertainty_rb_ratio_n = np.divide(
        np.linalg.norm(r_uncertainty_grad_n, axis=(1, 2)),
        np.linalg.norm(b_uncertainty_grad_n, axis=(1, 2)))
    return (sum_uncertainty_grad_n, sum_uncertainty_grad_f, sum_uncertainty_curve,
            sum_uncertainty_curl, sum_uncertainty_div, uncertainty_rb_ratio_n)


def build_dataframe(t_master, curve_0, sum_uncertainty_curve, bm_0, Bmag_0, r_i, r_e,
                    curl_0, sum_uncertainty_curl, div_0, sum_uncertainty_div,
                    uncertainty_rb_ratio_n, beta_i, beta_e, beta_total, v_perp_i):
    '''
    Build the output pandas DataFrame from all computed quantities.

    inputs:
    t_master:               master time series (used as DataFrame index)
    curve_0:                nominal curvature vector array (n, 3)
    sum_uncertainty_curve:  total curvature uncertainty array (n, 3)
    bm_0:                   mean normalized B vector array (n, 3)
    Bmag_0:                 mean |B| array (n,)
    r_i:                    ion gyroradius array in km (n,)
    r_e:                    electron gyroradius array in km (n,)
    curl_0:                 nominal curl of B array (n, 3)
    sum_uncertainty_curl:   total curl uncertainty array (n, 3)
    div_0:                  nominal divergence of B array (n,)
    sum_uncertainty_div:    total divergence uncertainty array (n,)
    uncertainty_rb_ratio_n: ratio of positional to magnetometer uncertainty norm (n,)
    beta_i:                 ion plasma beta array (n,)
    beta_e:                 electron plasma beta array (n,)
    beta_total:             total plasma beta array (n,)
    v_perp_i:               magnitude of perpendicular ion bulk velocity in km/s (n,)

    outputs:
    curvedf:    pandas DataFrame with Time index and all computed columns
    '''
    curve_norm = np.linalg.norm(curve_0, axis=1)
    curve_mag_error = np.sqrt(np.divide(
        np.square(np.add(
            np.multiply(curve_0[:, 0], sum_uncertainty_curve[:, 0]),
            np.add(
                np.multiply(curve_0[:, 1], sum_uncertainty_curve[:, 1]),
                np.multiply(curve_0[:, 2], sum_uncertainty_curve[:, 2])
            )
        )),
        np.square(curve_norm)
    ))
    rc_error = np.divide(curve_mag_error, np.power(curve_norm, 2))

    curvedf = pd.DataFrame({
        'Rc(km)':           1 / curve_norm,
        '|curve|':          curve_norm,
        'Curvature_X(GSM)': curve_0.take(0, axis=1),
        'Curvature_Y(GSM)': curve_0.take(1, axis=1),
        'Curvature_Z(GSM)': curve_0.take(2, axis=1),
        'error_Kx':         sum_uncertainty_curve.take(0, axis=1),
        'error_Ky':         sum_uncertainty_curve.take(1, axis=1),
        'error_Kz':         sum_uncertainty_curve.take(2, axis=1),
        'error_Rc':         rc_error,
        'b_x':              bm_0.take(0, axis=1),
        'b_y':              bm_0.take(1, axis=1),
        'b_z':              bm_0.take(2, axis=1),
        '|B|':              Bmag_0,
        'R_gi(km)':         r_i,
        'R_ge(km)':         r_e,
        'error_|curve|':    curve_mag_error,
        'error_r/b_ratio':  uncertainty_rb_ratio_n,
        'curlB_x':          curl_0.take(0, axis=1),
        'curlB_y':          curl_0.take(1, axis=1),
        'curlB_z':          curl_0.take(2, axis=1),
        'error_curlx':      sum_uncertainty_curl.take(0, axis=1),
        'error_curly':      sum_uncertainty_curl.take(1, axis=1),
        'error_curlz':      sum_uncertainty_curl.take(2, axis=1),
        'div(B)':           div_0,
        'error_div(B)':     sum_uncertainty_div,
        'beta_i':           beta_i,
        'beta_e':           beta_e,
        'beta_total':       beta_total,
        '|v_perp_i|(km/s)': v_perp_i,
    }, index=t_master)
    curvedf.index.name = "Time"
    return curvedf


def save_results(curvedf, filename, save_csv=True, save_h5=False):
    '''
    Save the results DataFrame to CSV and/or HDF5.

    inputs:
    curvedf:    pandas DataFrame to save
    filename:   output filename (expected to end in .csv)
    save_csv:   if True, save as CSV (default True)
    save_h5:    if True, also save as HDF5 with .h5 extension (default False)
    '''
    if save_csv:
        print("Writing File: " + filename)
        curvedf.to_csv(filename)
    if save_h5:
        curvedf.to_hdf(filename[:-3] + 'h5', key='df')


def _prompt(label, current):
    '''Display current value of a parameter and return user input or current value.'''
    response = input(f"  {label} [{current}]: ").strip()
    return response if response else current


def _prompt_bool(label, current):
    '''Display current boolean parameter and return user-confirmed value.'''
    display = 'y' if current else 'n'
    response = input(f"  {label} [{display}] (y/n): ").strip().lower()
    if response == '':
        return current
    return response == 'y'


def _confirm_parameters(trange, data_rate, prefix, suffix, save_csv, save_h5):
    '''
    Display current run parameters and allow the user to confirm or modify
    each one interactively.  Returns the (possibly updated) parameter set.
    '''
    print("\n--- Run Parameters ---")
    print(f"  trange start : {trange[0]}")
    print(f"  trange end   : {trange[1]}")
    print(f"  data_rate    : {data_rate}")
    print(f"  prefix       : {prefix}")
    print(f"  suffix       : {suffix}")
    print(f"  save_csv     : {save_csv}")
    print(f"  save_h5      : {save_h5}")
    print("\nPress Enter to keep each value, or type a new one.")

    trange    = [_prompt("trange start", trange[0]),
                 _prompt("trange end",   trange[1])]
    data_rate = _prompt("data_rate (srvy/brst)", data_rate)
    prefix    = _prompt("prefix",   prefix)
    suffix    = _prompt("suffix",   suffix)
    save_csv  = _prompt_bool("save_csv", save_csv)
    save_h5   = _prompt_bool("save_h5",  save_h5)

    print("\n--- Confirmed Parameters ---")
    print(f"  trange    : {trange}")
    print(f"  data_rate : {data_rate}")
    print(f"  prefix    : {prefix}")
    print(f"  suffix    : {suffix}")
    print(f"  save_csv  : {save_csv}")
    print(f"  save_h5   : {save_h5}")
    print()

    return trange, data_rate, prefix, suffix, save_csv, save_h5


def main():
    ####################################################
    # Default parameters
    trange     = ['2020-08-01', '2020-08-05']
    data_rate  = 'srvy'
    prefix     = "~/Work/MMS_work/CurveGSM_"
    suffix     = "_v6.2"
    save_csv   = True
    save_h5    = False
    num_probes = 4
    ####################################################

    trange, data_rate, prefix, suffix, save_csv, save_h5 = _confirm_parameters(
        trange, data_rate, prefix, suffix, save_csv, save_h5)

    timeStart = time.strftime("%H:%M:%S", time.localtime())
    print("Files Loading:")

    filename = generate_filename(trange, prefix, suffix)

    pos_times, b_times, pos_values, b_values = load_fgm_data(trange, data_rate, num_probes)
    fgm_load_done_time = time.strftime("%H:%M:%S", time.localtime())
    print("Time started: ", timeStart)
    print("Time FGM Loaded: ", fgm_load_done_time)

    fpidata, fpirate = load_fpi_data(trange, data_rate)
    print("Time FPI Loaded: ", time.strftime("%H:%M:%S", time.localtime()))

    calc_start_time = time.strftime("%H:%M:%S", time.localtime())
    print("Calculating Curvature:")

    grad_0n, grad_0f, bm_0, Bmag_0, rm_0, t_master, curve_0, curl_0, div_0 = calc_nominal(
        pos_times, pos_values, b_times, b_values)

    print("Collecting positional uncertainties...")
    outRerr = load_positional_uncertainty(trange, num_probes, t_master)

    r_unc = calc_positional_uncertainty(
        pos_times, pos_values, b_times, b_values, outRerr, t_master,
        grad_0n, grad_0f, curve_0, curl_0, div_0)

    b_unc = calc_magnetometer_uncertainty(
        pos_times, pos_values, b_times, b_values,
        grad_0n, grad_0f, curve_0, curl_0, div_0)

    (sum_uncertainty_grad_n, sum_uncertainty_grad_f,
     sum_uncertainty_curve, sum_uncertainty_curl,
     sum_uncertainty_div, uncertainty_rb_ratio_n) = combine_uncertainties(*r_unc, *b_unc)

    calc_end_time = time.strftime("%H:%M:%S", time.localtime())
    print("Done calculating Curvature.")

    print("Calculating gyroradii...")
    r_i, r_e = mesoGyroradius(fpidata=fpidata, fpirate=fpirate, t_master=t_master, bmag=Bmag_0)
    r_i = r_i / 1000
    r_e = r_e / 1000   # Convert from meters to km

    print("Calculating plasma beta...")
    beta_i, beta_e, beta_total = calc_plasma_beta(fpidata=fpidata, fpirate=fpirate,
                                                   t_master=t_master, Bmag_0=Bmag_0)

    print("Calculating perpendicular ion velocity...")
    v_perp_i = calc_perp_ion_velocity(fpidata=fpidata, fpirate=fpirate,
                                       t_master=t_master, bm_0=bm_0)
    print(f"  |v_perp_i| -- min: {v_perp_i.min():.2f} km/s  "
          f"mean: {v_perp_i.mean():.2f} km/s  "
          f"max: {v_perp_i.max():.2f} km/s")

    curvedf = build_dataframe(
        t_master, curve_0, sum_uncertainty_curve, bm_0, Bmag_0, r_i, r_e,
        curl_0, sum_uncertainty_curl, div_0, sum_uncertainty_div, uncertainty_rb_ratio_n,
        beta_i, beta_e, beta_total, v_perp_i)

    save_results(curvedf, filename, save_csv=save_csv, save_h5=save_h5)

    print("Time started: ", timeStart)
    print("Grad(B) products calculation start: ", calc_start_time)
    print("Grad(B) products calculation end:   ", calc_end_time)
    print("Time finished: ", time.strftime("%H:%M:%S", time.localtime()))


if __name__ == '__main__':
    main()
