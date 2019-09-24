'''
python3 + numpy routine to calculate magnetic curvature from MMS data.

Uses pyspedas for MMS data file loading

*****************************************************************************
Example script for loading required data and calculating the curvature vector:
    
import time
import numpy as np
from mms_curvature import DataLoad, Curvature

timeStart =  time.strftime("%H:%M:%S", time.localtime())
print("Files Loading:")

postime1, pos1, magtime1, mag1, postime2, pos2, magtime2, mag2, postime3, pos3, magtime3, mag3, postime4, pos4, magtime4, mag4 = DataLoad(trange=['2017-05-04', '2017-05-05'])

print("Time started: ", timeStart)
print("Time Loaded: ", time.strftime("%H:%M:%S", time.localtime()))

t_master, grad_Harvey, curve_Harvey = Curvature(postime1, pos1, magtime1, mag1, postime2, pos2, magtime2, mag2, postime3, pos3, magtime3, mag3, postime4, pos4, magtime4, mag4)

np.savetxt("t_master.csv", t_master, delimiter=",")
np.savetxt("curve_Harvey.csv", curve_Harvey, delimiter=",")
np.save("grad_Harvey.npy", grad_Harvey)

# end
****************************************************************************

---AJR 08.22.2019
'''
import numpy as np
import pyspedas
from pytplot import get_data

def DataLoad(trange=['2017-05-01', '2017-05-02/15:30:02'], data_rate='srvy', level='l2'):
    '''
    Loads all data needed for calculating magnnetic field curvature from MMS FGM data.  
    Uses pyspedas and pytplot.get_data for accessing the SDC API, file downloading, 
    data file version control, and CDF unpacking.

    Parameters:
    trange:     A list with two strings for the date range [tstart, tend]
                e.g. trange=['2017-05-01', '2017-05-02/15:30:02']

    data_rate:  The cadance of data which should be loaded.
                Options are 'srvy', 'brst'

    level:      The data level which will be loaded.  Use 'l2' unless you're sure otherwise.

    '''
    # load data files from SDC/local storage into tplot variables
    pyspedas.mms_load_mec(trange=trange, probe=['1', '2', '3', '4'], data_rate='srvy', level=level, time_clip=True)
    pyspedas.mms_load_fgm(trange=trange, probe=['1', '2', '3', '4'], data_rate=data_rate, level=level, time_clip=True)
    # extract data from tplot variables to numpy arrays.  NOTE: all done in GSM.
    postime1, pos1 = get_data('mms1_mec_r_gsm')
    postime2, pos2 = get_data('mms2_mec_r_gsm')
    postime3, pos3 = get_data('mms3_mec_r_gsm')
    postime4, pos4 = get_data('mms4_mec_r_gsm')
    magtime1, mag1 = get_data('mms1_fgm_b_gsm_'+data_rate+'_l2')
    magtime2, mag2 = get_data('mms2_fgm_b_gsm_'+data_rate+'_l2')
    magtime3, mag3 = get_data('mms3_fgm_b_gsm_'+data_rate+'_l2')
    magtime4, mag4 = get_data('mms4_fgm_b_gsm_'+data_rate+'_l2')
    # return all arrays
    return (postime1, pos1, magtime1, mag1, postime2, pos2, magtime2, mag2, postime3, pos3, magtime3, mag3, postime4, pos4, magtime4, mag4)


def Curvature(postime1, pos1, magtime1, mag1, postime2, pos2, magtime2, mag2, postime3, pos3, magtime3, mag3, postime4, pos4, magtime4, mag4, report_all=False): 
    '''
    Calculates spacial gradient and curvature vector of the magnetic field.  
    Returns those and the master time (interpolated from FGM data) as numpy arrays
    '''
    
    # normalize magnetic fields
#    bn1 = mag1[:,0:3]/np.linalg.norm(mag1[:,0:3], ord=2, axis=1, keepdims=True)
#    bn2 = mag2[:,0:3]/np.linalg.norm(mag2[:,0:3], ord=2, axis=1, keepdims=True)
#    bn3 = mag3[:,0:3]/np.linalg.norm(mag3[:,0:3], ord=2, axis=1, keepdims=True)
#    bn4 = mag4[:,0:3]/np.linalg.norm(mag4[:,0:3], ord=2, axis=1, keepdims=True)
    
    bn1 = mag1[:,0:3]/mag1[:,3,np.newaxis]
    bn2 = mag2[:,0:3]/mag2[:,3,np.newaxis]
    bn3 = mag3[:,0:3]/mag3[:,3,np.newaxis]
    bn4 = mag4[:,0:3]/mag4[:,3,np.newaxis]

    # find probe with latest beginning point for magnetic field data
    mastersc = np.argmax([magtime1[0], magtime2[0], magtime3[0], magtime4[0]])
    # find earliest end time for all space craft in range
    tend = min(magtime1[-1], magtime2[-1], magtime3[-1], magtime4[-1])
    tend_i = 0  # initialize counting index for finding ending index for master S/C
    
    # initialize master time sequence by trimming time from last S/C to start (mastersc)
    if mastersc == 0:
        while magtime1[tend_i] < tend: tend_i += 1
        t_master = magtime1[0:(tend_i+1)]
    elif mastersc == 1:
        while magtime2[tend_i] < tend: tend_i += 1
        t_master = magtime2[0:(tend_i+1)]
    elif mastersc == 2:
        while magtime3[tend_i] < tend: tend_i += 1
        t_master = magtime3[0:(tend_i+1)]
    elif mastersc == 3:
        while magtime4[tend_i] < tend: tend_i += 1
        t_master = magtime4[0:(tend_i+1)]
    
    # master mag field data arr, with interpolated values
    barr=np.ndarray((4,t_master.shape[0],3))
    
    barr[0,:,0] = np.interp(t_master, magtime1, bn1[:,0]) # MMS1
    barr[0,:,1] = np.interp(t_master, magtime1, bn1[:,1]) 
    barr[0,:,2] = np.interp(t_master, magtime1, bn1[:,2]) 
    
    barr[1,:,0] = np.interp(t_master, magtime2, bn2[:,0]) # MMS2
    barr[1,:,1] = np.interp(t_master, magtime2, bn2[:,1]) 
    barr[1,:,2] = np.interp(t_master, magtime2, bn2[:,2]) 
    
    barr[2,:,0] = np.interp(t_master, magtime3, bn3[:,0]) # MMS3
    barr[2,:,1] = np.interp(t_master, magtime3, bn3[:,1]) 
    barr[2,:,2] = np.interp(t_master, magtime3, bn3[:,2]) 
    
    barr[3,:,0] = np.interp(t_master, magtime4, bn4[:,0]) # MMS4
    barr[3,:,1] = np.interp(t_master, magtime4, bn4[:,1]) 
    barr[3,:,2] = np.interp(t_master, magtime4, bn4[:,2]) 
    
    # The below now done using time_clip=True in DataLoad
    '''
    # need to clear extranious data points from positional data for
    #   np.interp to work as intended.
    bcnt = 0        # MMS1
    while postime1[bcnt] < t_master[0] and bcnt: bcnt = bcnt + 1
    ecnt = bcnt
    while postime1[ecnt] <= t_master[-1]: ecnt = ecnt + 1
    pos1 = pos1[(bcnt-1):(ecnt+1),:]
    postime1 = postime1[(bcnt-1):(ecnt+1),:]
   
    bcnt = 0        # MMS2
    while postime2[bcnt] < t_master[0]: bcnt = bcnt + 1
    ecnt = bcnt
    while postime2[ecnt] <= t_master[-1]: ecnt = ecnt + 1
    pos2 = pos2[(bcnt-1):(ecnt+1),:]
    postime2 = postime2[(bcnt-1):(ecnt+1),:]

    bcnt = 0        # MMS3
    while postime3[bcnt] < t_master[0]: bcnt = bcnt + 1
    ecnt = bcnt
    while postime3[ecnt] <= t_master[-1]: ecnt = ecnt + 1
    pos3 = pos3[(bcnt-1):(ecnt+1),:]
    postime3 = postime3[(bcnt-1):(ecnt+1),:]
    
    bcnt = 0        # MMS4
    while postime4[bcnt] < t_master[0]: bcnt = bcnt + 1
    ecnt = bcnt
    while postime4[ecnt] <= t_master[-1]: ecnt = ecnt + 1
    pos4 = pos4[(bcnt-1):(ecnt+1),:]
    postime4 = postime4[(bcnt-1):(ecnt+1),:]
    ''' 
    
    # master position data array, with interpolated value
    rarr = np.ndarray((4,t_master.shape[0],3))
    
    rarr[0,:,0] = np.interp(t_master, postime1, pos1[:,0]) # MMS1
    rarr[0,:,1] = np.interp(t_master, postime1, pos1[:,1])
    rarr[0,:,2] = np.interp(t_master, postime1, pos1[:,2])
    
    rarr[1,:,0] = np.interp(t_master, postime2, pos2[:,0]) # MMS2
    rarr[1,:,1] = np.interp(t_master, postime2, pos2[:,1])
    rarr[1,:,2] = np.interp(t_master, postime2, pos2[:,2])
    
    rarr[2,:,0] = np.interp(t_master, postime3, pos3[:,0]) # MMS3
    rarr[2,:,1] = np.interp(t_master, postime3, pos3[:,1])
    rarr[2,:,2] = np.interp(t_master, postime3, pos3[:,2])
    
    rarr[3,:,0] = np.interp(t_master, postime4, pos4[:,0]) # MMS4
    rarr[3,:,1] = np.interp(t_master, postime4, pos4[:,1])
    rarr[3,:,2] = np.interp(t_master, postime4, pos4[:,2])
    
    # Now all magnetic fields and positional dataof of the same cadance and at the same times for each index
    # Indicies are: [s/c(0=mms1, 1=mms2, 2=mms3, 3=mms4), time_step, vector(0=x, 1=y, 2=z)]
    
    # calculate position and normalized magneic field at mesocenter of the fleet
    
    

    rm = np.ndarray((t_master.shape[0], 3))
    for i in range(t_master.shape[0]):      # populate each element of the mesocenter with the average across all four s/c
        for j in range(3):                  # there's got to be a way to vctorize this
            #print("rm:", rm.shape, "rarr:", rarr.shape)
            rm[i,j] = (1./4.)*(rarr[0,i,j]+rarr[1,i,j]+rarr[2,i,j]+rarr[3,i,j])
    
    bm = np.ndarray((t_master.shape[0], 3))
    for i in range(t_master.shape[0]):      # populate each element of the mesocenter with the average across all four s/c
        for j in range(3):                  # there's got to be a way to vctorize this
            bm[i,j] = (1./4.)*(barr[0,i,j]+barr[1,i,j]+barr[2,i,j]+barr[3,i,j])
    
    # Calculate volumetric tensor (Harvey, Ch 12.4, Eq 12.23, from "Analysis Methods for Multi-Spacecraft Data" Paschmann ed.)
    
    Rvol = np.ndarray((t_master.shape[0], 3, 3))
    for i in range(t_master.shape[0]):
        #rvol_step = np.zeros([3,3])    # Stepwise method gives same result as explicit below
        #for sc in range(4):
        #    rvol_step = rvol_step + np.outer(rarr[sc,i,:], rarr[sc,i,:])
        ## endfor
        #Rvol[i,:,:] = (1./4.) * (rvol_step - np.outer(rm[i,:], rm[i,:]))
        
        Rvol[i,:,:] = (1./4.)*((np.outer(rarr[0,i,:], rarr[0,i,:]) + np.outer(rarr[1,i,:], rarr[1,i,:]) + np.outer(rarr[2,i,:], rarr[2,i,:]) + np.outer(rarr[3,i,:], rarr[3,i,:])) - np.outer(rm[i,:], rm[i,:]))     # give same result as stepwise above
    
    Rinv = np.linalg.inv(Rvol)
    a_ne_b_list=[[1,2,3],[2,3],[3]]

    grad_Harvey = np.ndarray((t_master.shape[0], 3, 3))
    curve_Harvey = np.ndarray((t_master.shape[0], 3))    
    
    for t in range(t_master.shape[0]):      # steps through each time step
        for i in range(3):                  # for final i-component of the gradient 
            for j in range(3):              # for final j-component of the gradient
                dbdr = np.zeros(3)          # a != b summation row vector from Harvey.  Re-initialized as zeros for each i,j
                for k in range(3):          # step through k-index.  May be able to eliminate this with vectorization later
                    for a in range(3):      # step through spacecraft MMS1-3; MMS4 done implicitly
                        for b in a_ne_b_list[a]:
                            dbdr[k] = dbdr[k] + ((barr[a,t,i] - barr[b,t,i]) * (rarr[a,t,k] - rarr[b,t,k]))
                        # endfor
                    # endfor
                # endfor
                grad_Harvey[t,i,j] = (1./16.) * np.matmul(dbdr, Rinv[t,:,j])     # Gives the same result as below
                #grad_Harvey[t,i,j] = (1./16.) * np.matmul(dbdr, np.linalg.inv(Rvol[t,:,:])[:,j])    # Maybe linalg.inv doesn't vectorize the way I think?
            # endfor
        curve_Harvey[t,:] = np.matmul(bm[t,:], grad_Harvey[t,:,:])               # same thing, probably just should be matrix mult.
    # endfor
    
    if report_all:
        return (t_master, grad_Harvey, curve_Harvey, rarr, barr, rm, bm, Rvol, Rinv)  #used for troubleshooting
    else:
        return (t_master, grad_Harvey, curve_Harvey)
    





# end
