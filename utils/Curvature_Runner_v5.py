# Copyright 2020-2022 Anthony Rogers.  All rights reserved.
# Released under the Apache 2.0 license.

# This file/script serves as a reference use case for the library it
# is contained in.  Expect this file to be replaced/rewritten in the
# future to more accurately and completely reflect library usage.
# This file is currently retained only for historical reference.

import copy
import time
import numpy as np
import pandas as pd
from mms_curvature.mms_curvature import mms_Grad, mms_Curvature, mms_CurlB, mms_DivB
from mms_curvature.mms_load_data_shims import mms_load_fgm, mms_load_ancillary
from mms_curvature.utils.mms_gyroradius import DataLoadMoments, CalcRadius
# from mms_curvature.utils.mms_TQF import get_TQF
#from mms_curvature.utils.mms_Rerr import get_Rerr

# Minimum verbosity for status; time-keeping
timeStart =  time.strftime("%H:%M:%S", time.localtime())
print("Files Loading:")

####################################################
# Set parameters here
trange=['2020-08-02/16:40', '2020-08-02/17:30']
data_rate='brst'
prefix="~/Work/Curvature/testruns/CurveGSM_rg_"
suffix="_sigma_v5.2"
save_csv=True
save_h5=False

####################################################

## Calculate some additional attributes for later use
#
#
####################################################
#  Update this to use data shims!!! - ajr 15.09.20
####################################################
#
#
#
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

    # constants for future use
    me = 9.1094e-31     # electron mass in kilograms
    mp = 1.6726e-27     # proton mass in kg; used as proxy for all ions

    # load all the needed plasma data
    distime1, distempperp1, destime1, destempperp1 = DataLoadMoments(trange=trange, data_rate=data_rate, level='l2', probe='1')
    distime2, distempperp2, destime2, destempperp2 = DataLoadMoments(trange=trange, data_rate=data_rate, level='l2', probe='2')
    distime3, distempperp3, destime3, destempperp3 = DataLoadMoments(trange=trange, data_rate=data_rate, level='l2', probe='3')

    # assume time_clip worked and use mms1 as the master clock.  Might swap this for a search later
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
        tiperp = np.average([distempperp1, distempperp2, distempperp3, distempperp4], axis=0)   # average ion T_perp
        teperp = np.average([destempperp1, destempperp2, destempperp3, destempperp4], axis=0)   # average electron T_perp
    except:
        print('Error in loading mms4 data.  Will drop mms4 from this dataset.')
        tiperp = np.average([distempperp1, distempperp2, distempperp3], axis=0)   # average ion T_perp
        teperp = np.average([destempperp1, destempperp2, destempperp3], axis=0)   # average electron T_perp

    r_i = CalcRadius(part_time=mitime, part_tempperp=tiperp, b_time=t_master, b_mag=bmag, part_mass=mp, part_q=1.602177e-19)
    r_e = CalcRadius(part_time=metime, part_tempperp=teperp, b_time=t_master, b_mag=bmag, part_mass=me, part_q=1.602177e-19)

    return (r_i, r_e)








# def mesoGyroradius(trange=['2017-05-28', '2017-05-29'], data_rate='srvy', level='l2', t_master=None, bmag=None):
#     '''
#     Calculates the average gyroradius in the tetrahedron, assumed to be 
#     representative of the gyroradius at the mesocenter.  Uses FPI data.
# 
#     inputs:
#     trange:     2-element list of strings for start and end times 
#                 ex.['2017-05-28', '2017-05-29/12:02:01']
#     data_rate:  Choices are 'srvy' or 'brst'
#     level:      data level to use.  Use 'l2' unless you know for sure
#     t_master:   Time series for the magnetic field data, assumed from mms_curvature.Curvature
#     bmag:       |B| aligned with t_master.  Assumed from mms_curvature.Curvature
# 
#     outputs:
#     (r_i, r_e) -- 2-element structure aligned with t_master time series where:
#     r_i:        average gyroradius of ions (assumed protons) 
#     r_e:        average gyroradius of electrons
# 
#     '''
# 
#     # constants for future use
#     me = 9.1094e-31     # electron mass in kilograms
#     mp = 1.6726e-27     # proton mass in kg; used as proxy for all ions
# 
#     # load all the needed plasma data
#     distime1, distempperp1, destime1, destempperp1 = DataLoadMoments(trange=trange, data_rate=data_rate, level='l2', probe='1')
#     distime2, distempperp2, destime2, destempperp2 = DataLoadMoments(trange=trange, data_rate=data_rate, level='l2', probe='2')
#     distime3, distempperp3, destime3, destempperp3 = DataLoadMoments(trange=trange, data_rate=data_rate, level='l2', probe='3')
#     distime4, distempperp4, destime4, destempperp4 = DataLoadMoments(trange=trange, data_rate=data_rate, level='l2', probe='4')
# 
#     # assume time_clip worked and use mms1 as the master clock.  Might swap this for a search later
# 
#     mitime = distime1   # master ion time
#     distempperp2 = np.interp(mitime, distime2, distempperp2)
#     distempperp3 = np.interp(mitime, distime3, distempperp3)
#     distempperp4 = np.interp(mitime, distime4, distempperp4)
# 
#     metime = destime1   # master electron time
#     destempperp2 = np.interp(metime, destime2, destempperp2)
#     destempperp3 = np.interp(metime, destime3, destempperp3)
#     destempperp4 = np.interp(metime, destime4, destempperp4)
# 
#     tiperp = np.average([distempperp1, distempperp2, distempperp3, distempperp4], axis=0)   # average ion T_perp
#     teperp = np.average([destempperp1, destempperp2, destempperp3, destempperp4], axis=0)   # average electron T_perp
# 
#     r_i = CalcRadius(part_time=mitime, part_tempperp=tiperp, b_time=t_master, b_mag=bmag, part_mass=mp, part_q=1.602177e-19)
#     r_e = CalcRadius(part_time=metime, part_tempperp=teperp, b_time=t_master, b_mag=bmag, part_mass=me, part_q=1.602177e-19)
# 
#     return (r_i, r_e)
# 
############################################################################################################################


# Generate sane filename
if len(trange[0]) > 10 or len(trange[1]) > 10:
    filename = prefix+trange[0][:10]+"_"+trange[0][11:13]+trange[0][14:16]+"--"+trange[1][:10]+"_"+trange[1][11:13]+trange[1][14:16]+suffix+".csv"
else:
    filename = prefix+trange[0][:10]+"--"+trange[1][:10]+suffix+".csv"


# Get B-field data and calculate curvature
numProbes = 4

fgmdata = mms_load_fgm(trange=trange, probe=['1', '2', '3', '4'], data_rate=data_rate, time_clip=True)[0] 


pos_times = [None]*numProbes
b_times = [None]*numProbes
pos_values = [None]*numProbes
b_values = [None]*numProbes

#populate master arrays for reuse
for probe in range(numProbes):
    pos_times[probe] = np.copy(fgmdata['mms'+str(probe+1)+'_fgm_r_gsm_'+str(data_rate)+'_l2']['x'])
    b_times[probe] = np.copy(fgmdata['mms'+str(probe+1)+'_fgm_b_gsm_'+str(data_rate)+'_l2']['x'])
    pos_values[probe] = np.copy(fgmdata['mms'+str(probe+1)+'_fgm_r_gsm_'+str(data_rate)+'_l2']['y'])
    b_values[probe] = np.copy(fgmdata['mms'+str(probe+1)+'_fgm_b_gsm_'+str(data_rate)+'_l2']['y'])

fgm_load_done_time = time.strftime("%H:%M:%S", time.localtime())
print("Time started: ", timeStart)
print("Time FGM Loaded: ", fgm_load_done_time)

# Get DEFERR data fro positional uncertainty
print("Collecting positional uncertainties...")
deferr_in = mms_load_ancillary(probe=['1','2','3','4'], anc_product='deferr', trange=trange, time_clip=True)
Rerr_arr = [None]*numProbes
for probe in range(1,numProbes+1):
    Rerr_arr[(probe-1)] = (deferr_in[0]["MMS"+str(probe)+"_DEFERR"]).to_numpy()
# Returns positional uncertainty in METERS


tmpRerr = np.asarray(Rerr_arr)
outRerr = np.ndarray((numProbes,np.asarray(pos_times).shape[1],4)) 
for bird in range(numProbes): 
    for dim in range(4):
        outRerr[bird,:,dim] = np.interp(pos_times[bird],tmpRerr[bird][:,0], tmpRerr[bird][:,dim])

# Convert positional uncertainty in METERS to KILOMETERS
outRerr = 1e-3 * outRerr

calc_start_time = time.strftime("%H:%M:%S", time.localtime())
print("Calculating Curvature:")

# Calculate nominal values for grad(B) products; no uncertainty
grad_0n, bm_0, Bmag_0, rm_0, t_master = mms_Grad(postimes=pos_times, posvalues=pos_values, magtimes=b_times, magvalues=b_values, normalize=True)
grad_0f, Bm_0 = mms_Grad(postimes=pos_times, posvalues=pos_values, magtimes=b_times, magvalues=b_values, normalize=False)[:2]

curve_0 = mms_Curvature(grad_0n, bm_0)  # uses gradiant of normalized magnetic field

curl_0 = mms_CurlB(grad_0f)             # both Curl and Div use gradiant of 
div_0 = mms_DivB(grad_0f)               # NOT-normalized magnetic field


# Explicitly initialize, because it makes the later loops easier.
r_uncertainty_grad_n = np.zeros_like(grad_0n)
r_uncertainty_curve = np.zeros_like(curve_0)
r_uncertainty_grad_f = np.zeros_like(grad_0f)
r_uncertainty_curl = np.zeros_like(curl_0)
r_uncertainty_div = np.zeros_like(div_0)

b_uncertainty_grad_n = np.zeros_like(grad_0n)
b_uncertainty_curve = np.zeros_like(curve_0)
b_uncertainty_grad_f = np.zeros_like(grad_0f)
b_uncertainty_curl = np.zeros_like(curl_0)
b_uncertainty_div = np.zeros_like(div_0)

sum_uncertainty_grad_n = np.zeros_like(grad_0n)
sum_uncertainty_curve = np.zeros_like(curve_0)
sum_uncertainty_grad_f = np.zeros_like(grad_0f)
sum_uncertainty_curl = np.zeros_like(curl_0)
sum_uncertainty_div = np.zeros_like(div_0)

# positional uncertainties
tpos = copy.deepcopy(pos_values)
for probe in range(numProbes):
    for spatial_dim in range(3):
        for sign in [-1, 1]:
            # Apply uncertainty
            tpos[probe][:,spatial_dim] = np.add(pos_values[probe][:,spatial_dim], np.multiply(outRerr[probe,:,spatial_dim+1],sign))
        
            # Run curvature for current uncertainty
            grad_in, bm_i = mms_Grad(postimes=pos_times, posvalues=tpos, magtimes=b_times, magvalues=b_values, normalize=True)[:2]
            grad_if = mms_Grad(postimes=pos_times, posvalues=tpos, magtimes=b_times, magvalues=b_values, normalize=False)[0]
            
            curve_i = mms_Curvature(grad_in, bm_i)  # uses gradiant of normalized magnetic field
            
            curl_i = mms_CurlB(grad_if)             # both Curl and Div use gradiant of 
            div_i = mms_DivB(grad_if)               # NOT-normalized magnetic field
            
            # Add to uncertainty summation
            ### As I'm typing this, I'm suddenly unsure if the uncertainty sum should be element-wise or something else?
            ###  Assuming element-wise until told otherwise.
            r_uncertainty_grad_n = np.add(np.power(np.subtract(grad_in,grad_0n), 2), r_uncertainty_grad_n)
            r_uncertainty_curve = np.add(np.power(np.subtract(curve_i,curve_0), 2), r_uncertainty_curve)
            r_uncertainty_grad_f = np.add(np.power(np.subtract(grad_if,grad_0f), 2), r_uncertainty_grad_f)
            r_uncertainty_curl = np.add(np.power(np.subtract(curl_i,curl_0), 2), r_uncertainty_curl)
            r_uncertainty_div = np.add(np.power(np.subtract(div_i,div_0), 2), r_uncertainty_div)
            
            # Reset changed input
            tpos = copy.deepcopy(pos_values)

# magnetometer measurement uncertainties
tb = copy.deepcopy(b_values)
for probe in range(numProbes):
    for spatial_dim in range(3):
        # I think you said magnetometer uncertainty was +/- 0.1?  If I'm wrong, adjust the following loop values.
        for mag_uncertainty in [-0.1, 0.1]:
            # Apply uncertainty
            tb[probe][:,spatial_dim] = np.add(b_values[probe][:,spatial_dim], mag_uncertainty)
            
            # Run curvature for current uncertainty
            grad_in, bm_i = mms_Grad(postimes=pos_times, posvalues=pos_values, magtimes=b_times, magvalues=tb, normalize=True)[:2]
            grad_if = mms_Grad(postimes=pos_times, posvalues=pos_values, magtimes=b_times, magvalues=tb, normalize=False)[0]
            
            curve_i = mms_Curvature(grad_in, bm_i)  # uses gradiant of normalized magnetic field
            
            curl_i = mms_CurlB(grad_if)             # both Curl and Div use gradiant of 
            div_i = mms_DivB(grad_if)               # NOT-normalized magnetic field
            
            # Add to uncertainty summation
            ### As I'm typing this, I'm suddenly unsure if the uncertainty sum should be element-wise or something else?
            ###  Assuming element-wise until told otherwise.
            b_uncertainty_grad_n = np.add(np.power(np.subtract(grad_in,grad_0n), 2), b_uncertainty_grad_n)
            b_uncertainty_curve = np.add(np.power(np.subtract(curve_i,curve_0), 2), b_uncertainty_curve)
            b_uncertainty_grad_f = np.add(np.power(np.subtract(grad_if,grad_0f), 2), b_uncertainty_grad_f)
            b_uncertainty_curl = np.add(np.power(np.subtract(curl_i,curl_0), 2), b_uncertainty_curl)
            b_uncertainty_div = np.add(np.power(np.subtract(div_i,div_0), 2), b_uncertainty_div)
        
            # Reset changed input
            tb = copy.deepcopy(b_values)

# sum total uncertainties so far
sum_uncertainty_grad_n = np.sqrt(np.add(r_uncertainty_grad_n, b_uncertainty_grad_n))
sum_uncertainty_grad_f = np.sqrt(np.add(r_uncertainty_grad_f, b_uncertainty_grad_f))
sum_uncertainty_curve = np.sqrt(np.add(r_uncertainty_curve, b_uncertainty_curve))
sum_uncertainty_curl = np.sqrt(np.add(r_uncertainty_curl, b_uncertainty_curl))
sum_uncertainty_div = np.sqrt(np.add(r_uncertainty_div, b_uncertainty_div))

uncertainty_rb_ratio_n = np.divide(np.linalg.norm(r_uncertainty_grad_n, axis=(1,2)), np.linalg.norm(b_uncertainty_grad_n, axis=(1,2)))


# # Square-root the uncertainty arrays to get final magnitudes.
# r_uncertainty_grad  = np.sqrt(r_uncertainty_grad)
# r_uncertainty_curve = np.sqrt(r_uncertainty_curve)
# b_uncertainty_grad  = np.sqrt(b_uncertainty_grad)
# b_uncertainty_curve = np.sqrt(b_uncertainty_curve)
# sum_uncertainty_grad = np.sqrt(sum_uncertainty_grad)
# sum_uncertainty_curve = np.sqrt(sum_uncertainty_curve)

# t_master = f_0[0]
# curve_0 = f_0[2] 
# rm = f_0[5]
# bm = f_0[6]
# bmag = f_0[7]   
#    curldata = mms_bcurl(fields=[fgmdata['mms1_fgm_b_gsm_srvy_l2'], fgmdata['mms2_fgm_b_gsm_srvy_l2'], fgmdata['mms3_fgm_b_gsm_srvy_l2'], fgmdata['mms4_fgm_b_gsm_srvy_l2']], positions=[fgmdata['mms1_fgm_r_gsm_srvy_l2'], fgmdata['mms2_fgm_r_gsm_srvy_l2'], fgmdata['mms3_fgm_r_gsm_srvy_l2'], fgmdata['mms4_fgm_r_gsm_srvy_l2']])

calc_end_time = time.strftime("%H:%M:%S", time.localtime())
print("Done calculating Curvature.")


# Calculate gyroradii
print("Calculating gyroradii...")
r_i, r_e = mesoGyroradius(trange=trange, data_rate=data_rate, level='l2', t_master=t_master, bmag=Bmag_0)
r_i = r_i/1000
r_e = r_e/1000  # Convert from meters to km

# # Get fleet position
# posmltm, posrm = FleetPos(postime1, postime2, postime3, postime4, rm, t_master)
# 
# 
# # Get Tetrahedron Quality Factor (TQF)
# TQFarr = get_TQF(trange=trange)[1]
# bcnt=0
# while TQFarr[bcnt,0] < t_master[0]: bcnt = bcnt + 1     # find beginning of trange in TQFarr
# ecnt = -1
# while TQFarr[ecnt,0] > t_master[-1]: ecnt = ecnt -1     # find end of trange in TQFarr
# print("\n\n*** bcnt= "+str(bcnt)+"   ecnt= "+str(ecnt)+"   ***\n\n")
# TQF = np.interp(t_master, TQFarr[bcnt:ecnt,0].astype('float64'), TQFarr[bcnt:ecnt,1].astype('float64')) # interpolate to match t_master

# m0 = 4.0 * np.pi * 1e-7
# curlB_mag = np.multiply(jmag, (m0 * 1e12))


#e_indicator = np.divide(curldata['divB'], np.linalg.norm(curldata['curlB']))



# Calculate the magnitude of the curvature to find radius
curve_norm = np.linalg.norm(curve_0, axis=1)

curve_mag_error = np.sqrt(np.divide(np.square(np.add(np.multiply(curve_0[:,0], sum_uncertainty_curve[:,0]), np.add(np.multiply(curve_0[:,1], sum_uncertainty_curve[:,   1]), np.multiply(curve_0[:,2], sum_uncertainty_curve[:,2])))), np.square(curve_norm)))

rc_error = np.divide(curve_mag_error, np.power(curve_norm, 2))

# Construct Pandas data frame of all calculated data for export 

curvedf = pd.DataFrame({'Rc(km)': 1/curve_norm, '|curve|': curve_norm, 'Curvature_X(GSM)': curve_0.take(0,axis=1), 'Curvature_Y(GSM)': curve_0.take(1,axis=1),          'Curvature_Z(GSM)': curve_0.take(2,axis=1), 'error_Kx': sum_uncertainty_curve.take(0,axis=1), 'error_Ky':sum_uncertainty_curve.take(1,axis=1), 'error_Kz':sum_uncertainty_curve.take(2,axis=1), 'error_Rc':rc_error, 'b_x':bm_0.take(0, axis=1), 'b_y':bm_0.take(1, axis=1), 'b_z':bm_0.take(2, axis=1), '|B|':Bmag_0, 'R_gi(km)': r_i, 'R_ge(km)': r_e, 'error_|curve|': curve_mag_error, 'error_r/b_ratio':uncertainty_rb_ratio_n, 'curlB_x':curl_0.take(0,axis=1), 'curlB_y':curl_0.take(1, axis=1), 'curlB_z':curl_0.take(2, axis=1), 'error_curlx':sum_uncertainty_curl.take(0, axis=1), 'error_curly':sum_uncertainty_curl.take(1, axis=1), 'error_curlz':sum_uncertainty_curl.take(2, axis=1), 'div(B)':div_0, 'error_div(B)':sum_uncertainty_div}, index=t_master)


print("Writing File: "+filename)

curvedf.index.name="Time"
if save_csv: curvedf.to_csv(filename)
if save_h5: curvedf.to_hdf(filename[:-3]+'h5', key='df')




print("Time started: ", timeStart)
print("Grad(B) products calculation start: ", calc_start_time)
print("Grad(B) products calculation end:   ", calc_end_time)
print("Time finished: ", time.strftime("%H:%M:%S", time.localtime()))
# end
