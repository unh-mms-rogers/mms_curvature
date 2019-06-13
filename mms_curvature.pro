;
;
;
;

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
for t=0,n_elements(t_master)-1 do begin
	rvol_step = dblarr(3,3)
	for sc = 0,3 do begin		; making summation more explicit 
		rvol_step = rvol_step + (reform(rarr[sc,t,*]) # reform(rarr[sc,t,*]))
	endfor
	Rvol[t,*,*] = (1./4.) * (rvol_step - (reform(rm[t,*]) # reform(rm[t,*])))
;	Rvol[i,*,*] = (1./4.)*((reform(rarr[0,i,*]) # reform(rarr[0,i,*])) + (reform(rarr[1,i,*]) # reform(rarr[1,i,*])) + (reform(rarr[2,i,*]) # reform(rarr[2,i,*])) + (reform(rarr[3,i,*]) # reform(rarr[3,i,*]))) - (reform(rm[i,*]) # reform(rm[i,*]))	
	; reform() necessary to do matrix math on vectors in time series
	; # operator works on vectors A,B such that A # B = A # transpose(B), removing the need for an explicit transpose
endfor

Rinv = dblarr(n_elements(t_master), 3, 3, /nozero)
for i=0,n_elements(t_master)-1 do Rinv[i,*,*] = la_invert(reform(Rvol[i,*,*]), /double)
; Inverted R for later use

;stop

; Calculate grad(b) using both Harvey and Shen methods

grad_Harvey = {x:t_master, y:dblarr(n_elements(t_master), 3, 3, /nozero)}	; initialize grad(b) variable for Harvey, Shen methods
alt_grad_Harvey = {x:t_master, y:dblarr(n_elements(t_master), 3, 3, /nozero)}	; initialize grad(b) variable for Harvey, Shen methods
grad_Shen = {x:t_master, y:dblarr(n_elements(t_master), 3, 3, /nozero)}

k_harvey = {x:t_master, y:dblarr(n_elements(t_master), 3, /nozero)}		; initialize curvature vectors for Harvey, Shen methods
alt_k_harvey = {x:t_master, y:dblarr(n_elements(t_master), 3, /nozero)}		; initialize curvature vectors for Harvey, Shen methods
k_shen = {x:t_master, y:dblarr(n_elements(t_master), 3, /nozero)}

rc_harvey = {x:t_master, y:dblarr(n_elements(t_master), /nozero)}		; initialize radius of curvature for Harvey, Shen methods
alt_rc_harvey = {x:t_master, y:dblarr(n_elements(t_master), /nozero)}		; initialize radius of curvature for Harvey, Shen methods
rc_shen = {x:t_master, y:dblarr(n_elements(t_master), /nozero)}

a_ne_b_list = list([1,2,3],[2,3],[3])		; list of all unique combonations of 4 spacecraft for sc_a != sc_b

dbdr_tik = dblarr(n_elements(t_master), 3, 3)	; initialize the db dr central a != b summation in Harvey method for this time step and i-component

for t = 0,n_elements(t_master)-1 do begin		; Steps thorugh each time step
	for i =0,2 do begin			; for the final i-component of the gradient
		for j = 0,2 do begin		; for the final j-component of the gradient
			dbdr = dblarr(3)	; temp version to test
			for k = 0,2 do begin	; step through the k-index.  May not be needed if can be left for vector operation
				for a = 0,2 do begin	; step through spacecraft MMS1-3
					foreach b, a_ne_b_list[a] do begin	; With previous for loop covers all possible a != b combonaitons
						dbdr_tik[t,i,k] = dbdr_tik[t,i,k] + ((barr[a,t,i] - barr[b,t,i]) * (rarr[a,t,k] - rarr[b,t,k]))
						dbdr[k] = dbdr[k] + ((barr[a,t,i] - barr[b,t,i]) * (rarr[a,t,k] - rarr[b,t,k]))
					endforeach
				endfor
			endfor
			; Need to put rest of grad operation in here
			grad_Harvey.y[t,i,j] = (1./16.) * transpose(dbdr) # reform(Rinv[t,*,j])	; original index order
			alt_grad_Harvey.y[t,i,j] = (1./16.) * transpose(reform(dbdr_tik[t,i,*])) # reform(Rinv[t,*,j])	; reversing indicies used as trouble-shooting attempt
		endfor
	endfor
	; Put in any other time-step calculations here, e.g. b.grad(b)
	k_harvey.y[t,*] = transpose(bm[t,*]) # reform(grad_Harvey.y[t,*,*])
	rc_harvey.y[t] = 1./sqrt(k_harvey.y[t,0]^2 + k_harvey.y[t,1]^2 + k_harvey.y[t,2]^2)
	alt_k_harvey.y[t,*] = transpose(bm[t,*]) # reform(alt_grad_Harvey.y[t,*,*])
	alt_rc_harvey.y[t] = 1./sqrt(alt_k_harvey.y[t,0]^2 + alt_k_harvey.y[t,1]^2 + alt_k_harvey.y[t,2]^2)
endfor

store_data, 'b_curvature_vector_gsm', data=k_harvey
store_data, 'b_curvature_radius', data=rc_harvey
store_data, 'b_gradient_gsm', data=grad_Harvey
store_data, 'alt_curvature', data=alt_k_harvey
store_data, 'alt_radius', data=alt_rc_harvey
store_data, 'alt_gradiant', data=alt_grad_harvey

options, 'b_curvature_vector_gsm', 'labels', ['x_GSM', 'y_GSM', 'z_GSM']
options, 'b_curvature_vector_gsm', 'colors', [2, 4, 6]
options, 'b_curvature_vector_gsm', 'labflag', -1

options, 'b_curvature_radius', ytitle, "R_C (km)"
ylim, 'b_curvature_radius', 1.,1e7,1
ylim, 'alt_radius', 1., 1e7, 1

if makeplot ne 0 then begin
	tplot, trange=trange, ['b_curvature_vector_gsm', 'b_curvature_radius']
endif


print, "Stopping for troubleshooting manipulation of variables."
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
