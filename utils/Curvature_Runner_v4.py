import copy
import time
import numpy as np
import pandas as pd
# from mms_curvature import Curvature
from mms_curvature.mms_curvature import Curvature
# from mms_load_data_shims import mms_load_fgm
from mms_curvature.mms_load_data_shims import mms_load_fgm
#from mms_curvature.utils.mms_bcurl import mms_bcurl
from mms_curvature.utils.mms_gyroradius import DataLoadMoments, CalcRadius
# from utils.mms_gyroradius import DataLoadMoments, CalcRadius
#from mms_TQF import get_TQF
#from pytplot import get_data
from mms_curvature.utils.mms_Rerr import get_Rerr
# from utils.mms_Rerr import get_Rerr

# Minimum verbosity for status; time-keeping
timeStart =  time.strftime("%H:%M:%S", time.localtime())
print("Files Loading:")

####################################################
# Set parameters here
trange=['2017-06-17/20:16', '2017-06-17/20:34']
data_rate='brst'
prefix="~/Work/Curvature/testruns/CurveGSM_rg_"
suffix="_sigma_v4.1"
save_csv=True
save_h5=False

calc_curl=False
####################################################

## Calculate some additional attributes for later use
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
    t_master:   Time series for the magnetic field data, assumed from mms_curvature.Curvature
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
    distime4, distempperp4, destime4, destempperp4 = DataLoadMoments(trange=trange, data_rate=data_rate, level='l2', probe='4')

    # assume time_clip worked and use mms1 as the master clock.  Might swap this for a search later

    mitime = distime1   # master ion time
    distempperp2 = np.interp(mitime, distime2, distempperp2)
    distempperp3 = np.interp(mitime, distime3, distempperp3)
    distempperp4 = np.interp(mitime, distime4, distempperp4)

    metime = destime1   # master electron time
    destempperp2 = np.interp(metime, destime2, destempperp2)
    destempperp3 = np.interp(metime, destime3, destempperp3)
    destempperp4 = np.interp(metime, destime4, destempperp4)

    tiperp = np.average([distempperp1, distempperp2, distempperp3, distempperp4], axis=0)   # average ion T_perp
    teperp = np.average([destempperp1, destempperp2, destempperp3, destempperp4], axis=0)   # average electron T_perp

    r_i = CalcRadius(part_time=mitime, part_tempperp=tiperp, b_time=t_master, b_mag=bmag, part_mass=mp, part_q=1.602177e-19)
    r_e = CalcRadius(part_time=metime, part_tempperp=teperp, b_time=t_master, b_mag=bmag, part_mass=me, part_q=1.602177e-19)

    return (r_i, r_e)

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
    pos_values[probe] = copy.deepcopy(fgmdata['mms'+str(probe+1)+'_fgm_r_gsm_'+str(data_rate)+'_l2']['y'])
    b_values[probe] = copy.deepcopy(fgmdata['mms'+str(probe+1)+'_fgm_b_gsm_'+str(data_rate)+'_l2']['y'])


print("Time started: ", timeStart)
print("Time FGM Loaded: ", time.strftime("%H:%M:%S", time.localtime()))

# Get DEFERR data fro positional uncertainty
print("Collecting positional uncertainties...")
Rerr_arr = [None]*numProbes
for probe in range(1,numProbes+1):
    Rerr_arr[(probe-1)] = get_Rerr(trange=trange, probe=str(probe), datadir="~/data/mms/ancillary/mms"+str(probe)+"/deferr/")[1]

tmpRerr = np.asarray(Rerr_arr)
outRerr = np.ndarray((numProbes,np.asarray(pos_times).shape[1],4)) 
for bird in range(numProbes): 
    for dim in range(4):
        outRerr[bird,:,dim] = np.interp(pos_times[bird],tmpRerr[bird,:,0], tmpRerr[bird,:,dim])



print("Calculating Curvature:")
f_0 = Curvature( \
        pos_times[0], pos_values[0], b_times[0], b_values[0], \
        pos_times[1], pos_values[1], b_times[1], b_values[1], \
        pos_times[2], pos_values[2], b_times[2], b_values[2], \
        pos_times[3], pos_values[3], b_times[3], b_values[3], \
        report_all=True)

# Explicitly initialize, because it makes the later loops easier.
r_uncertainty_grad = np.zeros_like(f_0[1])
r_uncertainty_curve = np.zeros_like(f_0[2])
b_uncertainty_grad = np.zeros_like(f_0[1])
b_uncertainty_curve = np.zeros_like(f_0[2])
sum_uncertainty_grad = np.zeros_like(f_0[1])
sum_uncertainty_curve = np.zeros_like(f_0[2])

# positional uncertainties
tpos = copy.deepcopy(pos_values)
for probe in range(numProbes):
    for spatial_dim in range(3):
        for sign in [-1, 1]:
            # Apply uncertainty
            tpos[probe][:,spatial_dim] = np.add(pos_values[probe][:,spatial_dim], np.multiply(outRerr[probe,:,spatial_dim+1],sign))
        
        # Run curvature for current uncertainty
        f_i = Curvature( \
                pos_times[0], tpos[0], b_times[0], b_values[0], \
                pos_times[1], tpos[1], b_times[1], b_values[1], \
                pos_times[2], tpos[2], b_times[2], b_values[2], \
                pos_times[3], tpos[3], b_times[3], b_values[3], \
                )
        
        # Add to uncertainty summation
        ### As I'm typing this, I'm suddenly unsure if the uncertainty sum should be element-wise or something else?
        ###  Assuming element-wise until told otherwise.
        r_uncertainty_grad = np.add(np.power(np.subtract(f_i[1],f_0[1]), 2), r_uncertainty_grad)
        r_uncertainty_curve = np.add(np.power(np.subtract(f_i[2],f_0[2]), 2), r_uncertainty_curve)
        
        # Reset changed input
        tpos = copy.deepcopy(pos_values)
        # pos_values[probe][:,spatial_dim] = copy.deepcopy(fgmdata['mms'+str(probe+1)+'_fgm_r_gsm_'+str(data_rate)+'_l2']['y'][:,spatial_dim])

# magnetometer measurement uncertainties
tb = copy.deepcopy(b_values)
for probe in range(numProbes):
    for spatial_dim in range(3):
        # I think you said magnetometer uncertainty was +/- 0.1?  If I'm wrong, adjust the following loop values.
        for mag_uncertainty in [-0.1, 0.1]:
            # Apply uncertainty
            tb[probe][:,spatial_dim] = np.add(b_values[probe][:,spatial_dim], mag_uncertainty)
        
        # Run curvature for current uncertainty
        f_i = Curvature( \
                pos_times[0], pos_values[0], b_times[0], tb[0], \
                pos_times[1], pos_values[1], b_times[1], tb[1], \
                pos_times[2], pos_values[2], b_times[2], tb[2], \
                pos_times[3], pos_values[3], b_times[3], tb[3], \
                )
        
        # Add to uncertainty summation
        ###  Assuming element-wise until told otherwise.
        b_uncertainty_grad = np.add(np.power(np.subtract(f_i[1],f_0[1]), 2), b_uncertainty_grad)
        b_uncertainty_curve = np.add(np.power(np.subtract(f_i[2],f_0[2]), 2), b_uncertainty_curve)
        
        # Reset changed input
        tb = copy.deepcopy(b_values)
        # b_values[probe][:,spatial_dim] = copy.deepcopy(fgmdata['mms'+str(probe+1)+'_fgm_b_gsm_'+str(data_rate)+'_l2']['y'][:,spatial_dim])

# sum total uncertainties so far
sum_uncertainty_grad = np.add(r_uncertainty_grad, b_uncertainty_grad)
sum_uncertainty_curve = np.add(r_uncertainty_curve, b_uncertainty_curve)

# Square-root the uncertainty arrays to get final magnitudes.
r_uncertainty_grad  = np.sqrt(r_uncertainty_grad)
r_uncertainty_curve = np.sqrt(r_uncertainty_curve)
b_uncertainty_grad  = np.sqrt(b_uncertainty_grad)
b_uncertainty_curve = np.sqrt(b_uncertainty_curve)
sum_uncertainty_grad = np.sqrt(sum_uncertainty_grad)
sum_uncertainty_curve = np.sqrt(sum_uncertainty_curve)

t_master = f_0[0]
curve_Harvey = f_0[2] 
rm = f_0[5]
bm = f_0[6]
bmag = f_0[7]   
#    curldata = mms_bcurl(fields=[fgmdata['mms1_fgm_b_gsm_srvy_l2'], fgmdata['mms2_fgm_b_gsm_srvy_l2'], fgmdata['mms3_fgm_b_gsm_srvy_l2'], fgmdata['mms4_fgm_b_gsm_srvy_l2']], positions=[fgmdata['mms1_fgm_r_gsm_srvy_l2'], fgmdata['mms2_fgm_r_gsm_srvy_l2'], fgmdata['mms3_fgm_r_gsm_srvy_l2'], fgmdata['mms4_fgm_r_gsm_srvy_l2']])


print("Done calculating Curvature.")


# Calculate gyroradii
r_i, r_e = mesoGyroradius(trange=trange, data_rate=data_rate, level='l2', t_master=t_master, bmag=bmag)
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
curve_norm = np.linalg.norm(curve_Harvey, axis=1)

# curve_mag_error = np.multiply(np.linalg.norm(uk, axis=1), curve_norm)
# curve_mag_error = ( (uk_x * k_x) + (uk_y * k_y) + (uk_z * k_z) )  / |k| 

curve_mag_error = np.sqrt(np.divide(np.square(np.add(np.multiply(curve_Harvey[:,0], sum_uncertainty_curve[:,0]), np.add(np.multiply(curve_Harvey[:,1], sum_uncertainty_curve[:,   1]), np.multiply(curve_Harvey[:,2], sum_uncertainty_curve[:,2])))), np.square(curve_norm)))

rc_error = np.divide(curve_mag_error, np.power(curve_norm, 2))

# Construct Pandas data frame of all calculated data for export 

curvedf = pd.DataFrame({'Rc(km)': 1/curve_norm, '|curve|': curve_norm, 'Curvature_X(GSM)': curve_Harvey.take(0,axis=1), 'Curvature_Y(GSM)': curve_Harvey.take(1,axis=1),          'Curvature_Z(GSM)': curve_Harvey.take(2,axis=1), 'error_x': sum_uncertainty_curve.take(0,axis=1), 'error_y':sum_uncertainty_curve.take(1,axis=1), 'error_z':sum_uncertainty_curve.take(2,axis=1), 'error_Rc':rc_error, 'b_x':bm.take(0, axis=1), 'b_y':bm.take(1, axis=1), 'b_z':bm.take(2, axis=1), '|B|':bmag, 'R_gi(km)': r_i, 'R_ge(km)': r_e, 'error_|curve|': curve_mag_error}, index=t_master)

# curvedf = pd.DataFrame({'Rc(km)': 1/curve_norm, '|curve|': curve_norm, 'Curvature_X(GSM)': curve_Harvey.take(0,axis=1), 'Curvature_Y(GSM)': curve_Harvey.take(1,axis=1),'Curvature_Z(GSM)': curve_Harvey.take(2,axis=1), 'error_x': np.multiply(uk.take(0,axis=1), curve_Harvey.take(0,axis=1)), 'error_y':np.multiply(uk.take(1,axis=1), curve_Harvey.take(1,axis=1)), 'error_z':np.multiply(uk.take(2,axis=1), curve_Harvey.take(2,axis=1)), 'error_Rc':1/curve_mag_error, 'b_x':bm.take(0, axis=1), 'b_y':bm.take(1, axis=1), 'b_z':bm.take(2, axis=1), '|B|':bmag, 'R_gi(km)': r_i, 'R_ge(km)': r_e}, index=t_master)                    

# curvedf = pd.DataFrame({'Rc(km)': 1/curve_norm, '|curve|': curve_norm, 'Curvature_X(GSM)': curve_Harvey.take(0,axis=1), 'Curvature_Y(GSM)': curve_Harvey.take(1,axis=1), 'Curvature_Z(GSM)': curve_Harvey.take(2,axis=1), '|B|': bmag, 'R_gi(km)': r_i, 'R_ge(km)': r_e, 'Position(MLT)': posmltm, 'Position_radius(km)': posrm, 'Position_Xgsm(km)': rm.take(0,axis=1), 'Position_Ygsm(km)': rm.take(1,axis=1), 'Position_Zgsm(km)': rm.take(2,axis=1), 'TQF': TQF}, index=t_master)

print("Writing File: "+filename)

curvedf.index.name="Time"
if save_csv: curvedf.to_csv(filename)
if save_h5: curvedf.to_hdf(filename[:-3]+'h5', key='df')

if calc_curl:
    from mms_curvature.utils.mms_bcurl import mms_bcurl
    curlfields = [fgmdata['mms1_fgm_b_gsm_'+data_rate+'_l2'],fgmdata['mms2_fgm_b_gsm_'+data_rate+'_l2'],fgmdata['mms3_fgm_b_gsm_'+data_rate+'_l2'],fgmdata['mms4_fgm_b_gsm_'+data_rate+'_l2']]
    curlpos = [fgmdata['mms1_fgm_r_gsm_'+data_rate+'_l2'],fgmdata['mms2_fgm_r_gsm_'+data_rate+'_l2'],fgmdata['mms3_fgm_r_gsm_'+data_rate+'_l2'],fgmdata['mms4_fgm_r_gsm_'+data_rate+'_l2']]
    curldict = mms_bcurl(fields=curlfields, positions=curlpos)
    curldf = pd.DataFrame({'bx(nT)':curldict['barcentre'].take(0,axis=1),'by(nT)':curldict['barcentre'].take(1,axis=1), 'bz(nT)':curldict['barcentre'].take(2,axis=1), 'curlBx(nT/m)':curldict['curlB'].take(0,axis=1), 'curlBy(nT/m)':curldict['curlB'].take(1,axis=1), 'curlBz(nT/m)':curldict['curlB'].take(2,axis=1), '|J|(A/m^2)':curldict['jmag'], 'divB(nT/m)':curldict['divB']}, index=curldict['timeseries'])
    curldf.index.name="Time"
    if save_csv: curldf.to_csv("Curl_Products_"+filename)
    if save_h5: curldf.to_hdf("Curl_Products_"+filename[:-3]+'h5', key='df')



print("Time started: ", timeStart)
print("Time finished: ", time.strftime("%H:%M:%S", time.localtime()))
#mag_curve_frame.to_csv(filename)


#np.savetxt("t_master.csv", t_master, delimiter=",")
#np.savetxt("grad_Harvey.csv", grad_Harvey, delimiter=",")
#np.savetxt("curve_Harvey.csv", curve_Harvey, delimiter=",")
#np.save("grad_Harvey.npy", grad_Harvey)
# end
