# Copyright 2020-2022 Anthony Rogers.  All rights reserved.
# Released under the Apache 2.0 license.

import copy
import time
import numpy as np
import pandas as pd
from mms_curvature.mms_curvature import mms_Grad, mms_Curvature, mms_CurlB, mms_DivB
from mms_curvature.mms_load_data_shims import mms_load_fgm, mms_load_ancillary
from mms_curvature.utils.mms_gyroradius import DataLoadMoments, CalcRadius


def mesoGyroradius(trange=['2017-05-28', '2017-05-29'], data_rate='srvy', level='l2', t_master=None, bmag=None):
    '''
    Calculates the average gyroradius in the tetrahedron, assumed to be
    representative of the gyroradius at the mesocenter.  Uses FPI data.

    inputs:
    trange:     2-element list of strings for start and end times
                ex.['2017-05-28', '2017-05-29/12:02:01']
    data_rate:  Choices are 'srvy' or 'brst'
    level:      data level to use.  Use 'l2' unless you know for sure
    t_master:   Time series for the magnetic field data, assumed from
                mms_curvature.Curvature
    bmag:       |B| aligned with t_master.  Assumed from mms_curvature.Curvature

    outputs:
    (r_i, r_e) -- 2-element structure aligned with t_master time series where:
        r_i:        average gyroradius of ions (assumed protons)
        r_e:        average gyroradius of electrons
    '''

    me = 9.1094e-31     # electron mass in kilograms
    mp = 1.6726e-27     # proton mass in kg; used as proxy for all ions

    distime1, distempperp1, destime1, destempperp1 = DataLoadMoments(trange=trange, data_rate=data_rate, level='l2', probe='1')
    distime2, distempperp2, destime2, destempperp2 = DataLoadMoments(trange=trange, data_rate=data_rate, level='l2', probe='2')
    distime3, distempperp3, destime3, destempperp3 = DataLoadMoments(trange=trange, data_rate=data_rate, level='l2', probe='3')

    mitime = distime1   # master ion time
    distempperp2 = np.interp(mitime, distime2, distempperp2)
    distempperp3 = np.interp(mitime, distime3, distempperp3)

    metime = destime1   # master electron time
    destempperp2 = np.interp(metime, destime2, destempperp2)
    destempperp3 = np.interp(metime, destime3, destempperp3)

    try:
        distime4, distempperp4, destime4, destempperp4 = DataLoadMoments(trange=trange, data_rate=data_rate, level='l2', probe='4')
        distempperp4 = np.interp(mitime, distime4, distempperp4)
        destempperp4 = np.interp(metime, destime4, destempperp4)
        tiperp = np.average([distempperp1, distempperp2, distempperp3, distempperp4], axis=0)
        teperp = np.average([destempperp1, destempperp2, destempperp3, destempperp4], axis=0)
    except:
        print('Error in loading mms4 data.  Will drop mms4 from this dataset.')
        tiperp = np.average([distempperp1, distempperp2, distempperp3], axis=0)
        teperp = np.average([destempperp1, destempperp2, destempperp3], axis=0)

    r_i = CalcRadius(part_time=mitime, part_tempperp=tiperp, b_time=t_master, b_mag=bmag, part_mass=mp, part_q=1.602177e-19)
    r_e = CalcRadius(part_time=metime, part_tempperp=teperp, b_time=t_master, b_mag=bmag, part_mass=me, part_q=1.602177e-19)

    return (r_i, r_e)


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


def load_positional_uncertainty(trange, num_probes, pos_times):
    '''
    Load DEFERR ancillary positional uncertainty data and interpolate onto
    the FGM position time grid.

    inputs:
    trange:     2-element list of start/end time strings
    num_probes: number of MMS probes
    pos_times:  list of position time arrays (one per probe)

    outputs:
    outRerr:    array of shape (num_probes, n_times, 4) with positional
                uncertainty in kilometers
    '''
    deferr_in = mms_load_ancillary(probe=['1', '2', '3', '4'], anc_product='deferr', trange=trange, time_clip=True)
    Rerr_arr = [None] * num_probes
    for probe in range(1, num_probes + 1):
        Rerr_arr[probe - 1] = deferr_in[0]["MMS" + str(probe) + "_DEFERR"].to_numpy()

    tmpRerr = np.asarray(Rerr_arr)
    outRerr = np.ndarray((num_probes, np.asarray(pos_times).shape[1], 4))
    for bird in range(num_probes):
        for dim in range(4):
            outRerr[bird, :, dim] = np.interp(pos_times[bird], tmpRerr[bird][:, 0], tmpRerr[bird][:, dim])

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
                                grad_0n, grad_0f, curve_0, curl_0, div_0):
    '''
    Compute squared uncertainty contributions from positional errors (DEFERR).
    Each probe's position is perturbed by +/- its uncertainty in each spatial
    dimension and the resulting change in each product is accumulated.

    inputs:
    pos_times:  list of position time arrays (one per probe)
    pos_values: list of nominal position value arrays (one per probe)
    b_times:    list of magnetic field time arrays (one per probe)
    b_values:   list of magnetic field value arrays (one per probe)
    outRerr:    positional uncertainty array (num_probes, n_times, 4) in km
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

    tpos = copy.deepcopy(pos_values)
    for probe in range(num_probes):
        for spatial_dim in range(3):
            for sign in [-1, 1]:
                tpos[probe][:, spatial_dim] = np.add(
                    pos_values[probe][:, spatial_dim],
                    np.multiply(outRerr[probe, :, spatial_dim + 1], sign))

                grad_in, bm_i = mms_Grad(postimes=pos_times, posvalues=tpos, magtimes=b_times, magvalues=b_values, normalize=True)[:2]
                grad_if = mms_Grad(postimes=pos_times, posvalues=tpos, magtimes=b_times, magvalues=b_values, normalize=False)[0]
                curve_i = mms_Curvature(grad_in, bm_i)
                curl_i  = mms_CurlB(grad_if)
                div_i   = mms_DivB(grad_if)

                r_uncertainty_grad_n = np.add(np.power(np.subtract(grad_in, grad_0n), 2), r_uncertainty_grad_n)
                r_uncertainty_curve  = np.add(np.power(np.subtract(curve_i, curve_0), 2), r_uncertainty_curve)
                r_uncertainty_grad_f = np.add(np.power(np.subtract(grad_if, grad_0f), 2), r_uncertainty_grad_f)
                r_uncertainty_curl   = np.add(np.power(np.subtract(curl_i,  curl_0),  2), r_uncertainty_curl)
                r_uncertainty_div    = np.add(np.power(np.subtract(div_i,   div_0),   2), r_uncertainty_div)

                tpos = copy.deepcopy(pos_values)

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
                    uncertainty_rb_ratio_n):
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


def main():
    timeStart = time.strftime("%H:%M:%S", time.localtime())
    print("Files Loading:")

    ####################################################
    # Set parameters here
    trange     = ['2020-08-02/16:40', '2020-08-02/17:30']
    data_rate  = 'brst'
    prefix     = "~/Work/Curvature/testruns/CurveGSM_rg_"
    suffix     = "_sigma_v5.2"
    save_csv   = True
    save_h5    = False
    num_probes = 4
    ####################################################

    filename = generate_filename(trange, prefix, suffix)

    pos_times, b_times, pos_values, b_values = load_fgm_data(trange, data_rate, num_probes)
    fgm_load_done_time = time.strftime("%H:%M:%S", time.localtime())
    print("Time started: ", timeStart)
    print("Time FGM Loaded: ", fgm_load_done_time)

    print("Collecting positional uncertainties...")
    outRerr = load_positional_uncertainty(trange, num_probes, pos_times)

    calc_start_time = time.strftime("%H:%M:%S", time.localtime())
    print("Calculating Curvature:")

    grad_0n, grad_0f, bm_0, Bmag_0, rm_0, t_master, curve_0, curl_0, div_0 = calc_nominal(
        pos_times, pos_values, b_times, b_values)

    r_unc = calc_positional_uncertainty(
        pos_times, pos_values, b_times, b_values, outRerr,
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
    r_i, r_e = mesoGyroradius(trange=trange, data_rate=data_rate, level='l2', t_master=t_master, bmag=Bmag_0)
    r_i = r_i / 1000
    r_e = r_e / 1000   # Convert from meters to km

    curvedf = build_dataframe(
        t_master, curve_0, sum_uncertainty_curve, bm_0, Bmag_0, r_i, r_e,
        curl_0, sum_uncertainty_curl, div_0, sum_uncertainty_div, uncertainty_rb_ratio_n)

    save_results(curvedf, filename, save_csv=save_csv, save_h5=save_h5)

    print("Time started: ", timeStart)
    print("Grad(B) products calculation start: ", calc_start_time)
    print("Grad(B) products calculation end:   ", calc_end_time)
    print("Time finished: ", time.strftime("%H:%M:%S", time.localtime()))


if __name__ == '__main__':
    main()
