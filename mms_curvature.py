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

    # find probe with latest beginning point for magnetic field data and what time
    # tstart = max(magtime1[0], magtime2[0], magtime3[0], magtime4[0])
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
        rvol_step = np.zeros([3,3])
        #for sc in range(4):
        #    rvol_step = rvol_step + np.outer(rarr[sc,i,:], rarr[sc,i,:])
        # endfor
        #Rvol[i,:,:] = (1./4.) * (rvol_step - np.outer(rm[i,:], rm[i,:]))
        
        Rvol[i,:,:] = (1./4.)*((np.outer(rarr[0,i,:], rarr[0,i,:]) + np.outer(rarr[1,i,:], rarr[1,i,:]) + np.outer(rarr[2,i,:], rarr[2,i,:]) + np.outer(rarr[3,i,:], rarr[3,i,:])) - np.outer(rm[i,:], rm[i,:]))
    
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
                grad_Harvey[t,i,j] = (1./16.) * np.matmul(dbdr, Rinv[t,:,j])     # possibly shouldn't be np.outer... just matrix mult?
            # endfor
        curve_Harvey[t,:] = np.matmul(bm[t,:], grad_Harvey[t,:,:])               # same thing, probably just should be matrix mult.
    # endfor
    
    if report_all:
        return (t_master, grad_Harvey, curve_Harvey, rarr, barr, rm, bm, Rvol, Rinv)  #used for troubleshooting
    else:
        return (t_master, grad_Harvey, curve_Harvey)
    





# end
'''












;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro mms_curvature, trange=trange, data_rate=data_rate, makeplot=makeplot, level=level

if undefined(data_rate) then data_rate='srvy'
if undefined(makeplot) then makeplot=1
if undefined(level) then level='l2'


; Load MMS Ephemeris and Coordinate data
mms_load_mec, trange=trange, probe = [1,2,3,4], level=level, data_rate=date_rate, /time_clip	
;mms_qcotrans, 'mms1_mec_r_eci', 'mms1_mec_r_gsm', out_coord='gsm'	; Convert position to GSM for consistancy
;mms_qcotrans, 'mms2_mec_r_eci', 'mms2_mec_r_gsm', out_coord='gsm'
;mms_qcotrans, 'mms3_mec_r_eci', 'mms3_mec_r_gsm', out_coord='gsm'
;mms_qcotrans, 'mms4_mec_r_eci', 'mms4_mec_r_gsm', out_coord='gsm'

get_data, 'mms1_mec_r_gsm', data=rr1	; Extract data from tplot variable into structure
get_data, 'mms2_mec_r_gsm', data=rr2
get_data, 'mms3_mec_r_gsm', data=rr3
get_data, 'mms4_mec_r_gsm', data=rr4

; Load magnetic field data from the combined fluxgate data product.
mms_load_fgm, trange=trange, probe=[1,2,3,4], level=level, data_rate=date_rate, /time_clip

get_data, 'mms1_fgm_b_gsm_'+data_rate+'_'+level, data=bt1	; Base product provide component and magnitude in the structure
get_data, 'mms2_fgm_b_gsm_'+data_rate+'_'+level, data=bt2
get_data, 'mms3_fgm_b_gsm_'+data_rate+'_'+level, data=bt3
get_data, 'mms4_fgm_b_gsm_'+data_rate+'_'+level, data=bt4


; Generate normalized B-field vectors (b = B/|B|  ->  b = bt/|bt|)
b1 = {x:bt1.x, y:dblarr(n_elements(bt1.x), 3, /nozero)}
for i=0,n_elements(bt1.x)-1 do begin			; Step through all time steps in B data
	for j=0,2 do begin				; Step through all 3 components in each time step
		b1.y[i,j] = bt1.y[i,j]/bt1.y[i,3]	; divide value of component by magnitude
	endfor
endfor

b2 = {x:bt2.x, y:dblarr(n_elements(bt2.x), 3, /nozero)}
for i=0,n_elements(bt2.x)-1 do begin			; Step through all time steps in B data
	for j=0,2 do begin				; Step through all 3 components in each time step
		b2.y[i,j] = bt2.y[i,j]/bt2.y[i,3]	; divide value of component by magnitude
	endfor
endfor

b3 = {x:bt3.x, y:dblarr(n_elements(bt3.x), 3, /nozero)}
for i=0,n_elements(bt3.x)-1 do begin			; Step through all time steps in B data
	for j=0,2 do begin				; Step through all 3 components in each time step
		b3.y[i,j] = bt3.y[i,j]/bt3.y[i,3]	; divide value of component by magnitude
	endfor
endfor

b4 = {x:bt4.x, y:dblarr(n_elements(bt4.x), 3, /nozero)}
for i=0,n_elements(bt4.x)-1 do begin			; Step through all time steps in B data
	for j=0,2 do begin				; Step through all 3 components in each time step
		b4.y[i,j] = bt4.y[i,j]/bt4.y[i,3]	; divide value of component by magnitude
	endfor
endfor

; Find master time series from magnetic field data

;	; Find last beginning time --> master start time
;	; Find first end time --> master end time
;	; Use time from s/c with last start time as base, truncate all after master end time
;	; truncate and interpolate data from all other s/c to match master time base
;	; interpolate positional data from all s/c to master time base

t_start = max([b1.x[0], b2.x[0], b3.x[0], b4.x[0]], master_sc)	; find latest start time and which probe has it
t_end = min([b1.x[-1], b2.x[-1], b3.x[-1], b4.x[-1]])		; find first end time
t_end_index = 0

case master_sc of		; Use the time coords of the last start probe and then truncate at the first end time
	0: begin
		while b1.x[t_end_index] lt t_end do t_end_index += 1
		t_master = b1.x[0:t_end_index]
	end
	1: begin
		while b2.x[t_end_index] lt t_end do t_end_index += 1
		t_master = b2.x[0:t_end_index]
	end
	2: begin
		while b3.x[t_end_index] lt t_end do t_end_index += 1
		t_master = b3.x[0:t_end_index]
	end
	3: begin
		while b4.x[t_end_index] lt t_end do t_end_index += 1
		t_master = b4.x[0:t_end_index]
	end
endcase

; interpolated b timeseries, linear interpolation used


barr = dblarr(4,n_elements(t_master),3) 	; Combined array for later stepping.  Could be condensed for memory performance later
;bi1 = [[interpol(b1.y[*,0], b1.x, t_master)], [interpol(b1.y[*,1], b1.x, t_master)], [interpol(b1.y[*,2], b1.x, t_master)]]
barr[0,*,*] = [[interpol(b1.y[*,0], b1.x, t_master)], [interpol(b1.y[*,1], b1.x, t_master)], [interpol(b1.y[*,2], b1.x, t_master)]]
;bi2 = [[interpol(b2.y[*,0], b2.x, t_master)], [interpol(b2.y[*,1], b2.x, t_master)], [interpol(b2.y[*,2], b2.x, t_master)]]
barr[1,*,*] = [[interpol(b2.y[*,0], b2.x, t_master)], [interpol(b2.y[*,1], b2.x, t_master)], [interpol(b2.y[*,2], b2.x, t_master)]]
;bi3 = [[interpol(b3.y[*,0], b3.x, t_master)], [interpol(b3.y[*,1], b3.x, t_master)], [interpol(b3.y[*,2], b3.x, t_master)]]
barr[2,*,*] = [[interpol(b3.y[*,0], b3.x, t_master)], [interpol(b3.y[*,1], b3.x, t_master)], [interpol(b3.y[*,2], b3.x, t_master)]]
;bi4 = [[interpol(b4.y[*,0], b4.x, t_master)], [interpol(b4.y[*,1], b4.x, t_master)], [interpol(b4.y[*,2], b4.x, t_master)]]
barr[3,*,*] = [[interpol(b4.y[*,0], b4.x, t_master)], [interpol(b4.y[*,1], b4.x, t_master)], [interpol(b4.y[*,2], b4.x, t_master)]]

;stop

; Interpolate position data to increase cadance to that of the FGM data.  Linear interpolation used.
rarr = dblarr(4,n_elements(t_master), 3) 	; Combined array for later stepping.  Could be condensed for memory performance later
;r1 = [[interpol(rr1.y[*,0], rr1.x, t_master)], [interpol(rr1.y[*,1], rr1.x, t_master)], [interpol(rr1.y[*,2], rr1.x, t_master)]]
rarr[0,*,*] = [[interpol(rr1.y[*,0], rr1.x, t_master)], [interpol(rr1.y[*,1], rr1.x, t_master)], [interpol(rr1.y[*,2], rr1.x, t_master)]]
;r2 = [[interpol(rr2.y[*,0], rr2.x, t_master)], [interpol(rr1.y[*,1], rr2.x, t_master)], [interpol(rr2.y[*,2], rr2.x, t_master)]]
rarr[1,*,*] = [[interpol(rr2.y[*,0], rr2.x, t_master)], [interpol(rr2.y[*,1], rr2.x, t_master)], [interpol(rr2.y[*,2], rr2.x, t_master)]]
;r3 = [[interpol(rr3.y[*,0], rr3.x, t_master)], [interpol(rr1.y[*,1], rr3.x, t_master)], [interpol(rr3.y[*,2], rr3.x, t_master)]]
rarr[2,*,*] = [[interpol(rr3.y[*,0], rr3.x, t_master)], [interpol(rr3.y[*,1], rr3.x, t_master)], [interpol(rr3.y[*,2], rr3.x, t_master)]]
;r4 = [[interpol(rr4.y[*,0], rr4.x, t_master)], [interpol(rr1.y[*,1], rr4.x, t_master)], [interpol(rr4.y[*,2], rr4.x, t_master)]]
rarr[3,*,*] = [[interpol(rr4.y[*,0], rr4.x, t_master)], [interpol(rr4.y[*,1], rr4.x, t_master)], [interpol(rr4.y[*,2], rr4.x, t_master)]]

; Now all of the magnetic field and positional data are of the same cadance and at the same times for each time index
; Indicies are: [s/c(0=mms1, ... ,3=mms4), time_step, vector(x=0,y=1,z=2)]

; Calculate position and normalized magnetic field at mesocenter of fleet

rm = dblarr(n_elements(t_master), 3, /nozero)
for i=0,n_elements(t_master)-1 do begin
	rm[i,*] = (1.0/4.0) * (rarr[0,i,*] + rarr[1,i,*] + rarr[2,i,*] + rarr[3,i,*])
	for j=0,2 do begin
		rm[i,j] = (1./4.)*(rarr[0,i,j] + rarr[1,i,j] + rarr[2,i,j] + rarr[3,i,j])
	endfor
endfor

bm = dblarr(n_elements(t_master), 3, /nozero)
for i=0,n_elements(t_master)-1 do begin
	for j=0,2 do begin
		bm[i,j] = (1./4.)*(barr[0,i,j] + barr[1,i,j] + barr[2,i,j] + barr[3,i,j])
	endfor
endfor

; Calculate Volumetric Tensor (Harvey, Ch 12.4, Eq 12.23)

Rvol = dblarr(n_elements(t_master), 3, 3, /nozero)
for i=0,n_elements(t_master)-1 do begin
	Rvol[i,*,*] = (1./4.)*((reform(rarr[0,i,*]) # reform(rarr[0,i,*])) + (reform(rarr[1,i,*]) # reform(rarr[1,i,*])) + (reform(rarr[2,i,*]) # reform(rarr[2,i,*])) + (reform(rarr[3,i,*]) # reform(rarr[3,i,*]))) - (reform(rm[i,*]) # reform(rm[i,*]))	
	; reform() necessary to do matrix math on vectors in time series
	; # operator works on vectors A,B such that A # B = A # transpose(B), removing the need for an explicit transpose
endfor

Rinv = dblarr(n_elements(t_master), 3, 3, /nozero)
for i=0,n_elements(t_master)-1 do Rinv[i,*,*] = la_invert(reform(Rvol[i,*,*]), /double)
; Inverted R for later use

;stop

; Calculate grad(b) using both Harvey and Shen methods

grad_Harvey = {x:t_master, y:dblarr(n_elements(t_master), 3, 3, /nozero)}	; initialize grad(b) variable for Harvey, Shen methods
grad_Shen = {x:t_master, y:dblarr(n_elements(t_master), 3, 3, /nozero)}

k_harvey = {x:t_master, y:dblarr(n_elements(t_master), 3, /nozero)}		; initialize curvature vectors for Harvey, Shen methods
k_shen = {x:t_master, y:dblarr(n_elements(t_master), 3, /nozero)}

rc_harvey = {x:t_master, y:dblarr(n_elements(t_master), /nozero)}		; initialize radius of curvature for Harvey, Shen methods
rc_shen = {x:t_master, y:dblarr(n_elements(t_master), /nozero)}

a_ne_b_list = list([1,2,3],[2,3],[3])		; list of all unique combonations of 4 spacecraft for sc_a != sc_b

dbdr_tik = dblarr(n_elements(t_master), 3, 3)	; initialize the db dr central a != b summation in Harvey method for this time step and i-component

for t = 0,n_elements(t_master)-1 do begin		; Steps thorugh each time step
	for i =0,2 do begin			; for the final i-component of the gradient
		for j = 0,2 do begin		; for the final j-component of the gradient
			for k = 0,2 do begin	; step through the k-index.  May not be needed if can be left for vector operation
				for a = 0,2 do begin	; step through spacecraft MMS1-3
					foreach b, a_ne_b_list[a] do begin	; With previous for loop covers all possible a != b combonaitons
						dbdr_tik[t,i,k] = dbdr_tik[t,i,k] + ((barr[a,t,i] - barr[b,t,i]) * (rarr[a,t,k] - rarr[b,t,k]))
					endforeach
				endfor
			endfor
			; Need to put rest of grad operation in here
;			grad_Harvey.y[t,i,j] = (1./16.) * transpose(dbdr_tik) # reform(Rinv[t,*,j])	; original index order
			grad_Harvey.y[t,i,j] = (1./16.) * transpose(reform(dbdr_tik[t,i,*])) ## reform(Rinv[t,*,j])	; reversing indicies used as trouble-shooting attempt
		endfor
	endfor
	; Put in any other time-step calculations here, e.g. b.grad(b)
	k_harvey.y[t,*] = transpose(bm[t,*]) # reform(grad_Harvey.y[t,*,*])
	rc_harvey.y[t] = 1./sqrt(k_harvey.y[t,0]^2 + k_harvey.y[t,1]^2 + k_harvey.y[t,2]^2)
endfor

store_data, 'b_curvature_vector_gsm', data=k_harvey
store_data, 'b_curvature_radius', data=rc_harvey
store_data, 'b_gradient_gsm', data=grad_Harvey


stop
						; Balls.  I need to expand out the whole of the 6 permutations.  This will wait until tomorrow.

						; Need to endfor all of these.

; need to solve the implicit summation over k as well as every bloody other thing.  Might be able to do this operation vectorized over all of the time axis, but might not be able to.  In any case, this needs to get ironed out for both the Harvey and Shen methods for grad(b)





;for i=0,2 do begin
;	for j=0,2 do begin
;		for k = 0,2 do begin
;			grad_Harvey.y[*,i,j] = (bi1.y[*,i] - bi2.y[*,i])*(r1.y[*,k]
;			grad_Shen.y[i,j,k] = 
;		endfor
;	endfor
;endfor

; Calculate b curvature vectors and radii, as well as normal of the osculating plane using Harvey and Shen methods


end
'''
