;
;
;
;

pro mms_curvature, trange=trange, data_rate=data_rate, makeplot=makeplot, level=level

if undefined(data_rate) then data_rate='srvy'
if undefined(makeplot) then makeplot=1
if undefined(level) then level='l2'


; Load MMS Ephemeris and Coordinate data
mms_load_mec, trange=trange, probe = [1,2,3,4], level=level, data_rate=data_rate, /time_clip	

get_data, 'mms1_mec_r_gsm', data=rr1	; Extract data from tplot variable into structure
get_data, 'mms2_mec_r_gsm', data=rr2
get_data, 'mms3_mec_r_gsm', data=rr3
get_data, 'mms4_mec_r_gsm', data=rr4

; Load magnetic field data from the combined fluxgate data product.
mms_load_fgm, trange=trange, probe=[1,2,3,4], level=level, data_rate=data_rate, /time_clip

get_data, 'mms1_fgm_b_gsm_'+data_rate+'_'+level, data=bt1	; Base product provide component and magnitude in the structure
get_data, 'mms2_fgm_b_gsm_'+data_rate+'_'+level, data=bt2
get_data, 'mms3_fgm_b_gsm_'+data_rate+'_'+level, data=bt3
get_data, 'mms4_fgm_b_gsm_'+data_rate+'_'+level, data=bt4


; Generate normalized B-field vectors (b = B/|B|  ->  b = bt/|bt|)
; 	NOTE: The fgm 'b' data product contains 4 components: Bx,By,Bz,|B| so we may normalize by dividing by the last (3rd) index
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
; Interpolation is used to homogonize time steps across all four spacecraft in the MMS fleet.

; barr == B-Array.  This is a single array with all the normalized magnetic field data from all four spacecraft within a single structure 
;		    for easy manipulation via indicies later on.

barr = dblarr(4,n_elements(t_master),3) 	; Combined array for later stepping.  Could be condensed for memory performance later
barr[0,*,*] = [[interpol(b1.y[*,0], b1.x, t_master)], [interpol(b1.y[*,1], b1.x, t_master)], [interpol(b1.y[*,2], b1.x, t_master)]]
barr[1,*,*] = [[interpol(b2.y[*,0], b2.x, t_master)], [interpol(b2.y[*,1], b2.x, t_master)], [interpol(b2.y[*,2], b2.x, t_master)]]
barr[2,*,*] = [[interpol(b3.y[*,0], b3.x, t_master)], [interpol(b3.y[*,1], b3.x, t_master)], [interpol(b3.y[*,2], b3.x, t_master)]]
barr[3,*,*] = [[interpol(b4.y[*,0], b4.x, t_master)], [interpol(b4.y[*,1], b4.x, t_master)], [interpol(b4.y[*,2], b4.x, t_master)]]


; Interpolate position data to increase cadance to that of the FGM data.  Linear interpolation used.

; rarr == Position(r) Array.	This is a single array with all position information (in km) from all four spacecraft within a single
;				for easy manipulation vi indicies later on.

rarr = dblarr(4,n_elements(t_master), 3) 	; Combined array for later stepping.  Could be condensed for memory performance later
rarr[0,*,*] = [[interpol(rr1.y[*,0], rr1.x, t_master)], [interpol(rr1.y[*,1], rr1.x, t_master)], [interpol(rr1.y[*,2], rr1.x, t_master)]]
rarr[1,*,*] = [[interpol(rr2.y[*,0], rr2.x, t_master)], [interpol(rr2.y[*,1], rr2.x, t_master)], [interpol(rr2.y[*,2], rr2.x, t_master)]]
rarr[2,*,*] = [[interpol(rr3.y[*,0], rr3.x, t_master)], [interpol(rr3.y[*,1], rr3.x, t_master)], [interpol(rr3.y[*,2], rr3.x, t_master)]]
rarr[3,*,*] = [[interpol(rr4.y[*,0], rr4.x, t_master)], [interpol(rr4.y[*,1], rr4.x, t_master)], [interpol(rr4.y[*,2], rr4.x, t_master)]]

; Now all of the magnetic field and positional data are of the same cadance and at the same times for each time index
; Indicies for barr and rarr are: [s/c(0=mms1, ... ,3=mms4), time_step, vector(x=0,y=1,z=2)]

; Calculate position and normalized magnetic field at mesocenter of fleet

rm = dblarr(n_elements(t_master), 3, /nozero)	; rm == Position of mesocenter of fleet
for t=0,n_elements(t_master)-1 do begin		; time steps
;	rm[t,*] = (1.0/4.0) * (rarr[0,t,*] + rarr[1,t,*] + rarr[2,t,*] + rarr[3,t,*])
	for j=0,2 do begin
		rm[t,j] = (1./4.)*(rarr[0,t,j] + rarr[1,t,j] + rarr[2,t,j] + rarr[3,t,j])
	endfor
endfor

bm = dblarr(n_elements(t_master), 3, /nozero)	; bm == normalized magnetic field at mesocenter of fleet
for t=0,n_elements(t_master)-1 do begin		; time steps
	for j=0,2 do begin
		bm[t,j] = (1./4.)*(barr[0,t,j] + barr[1,t,j] + barr[2,t,j] + barr[3,t,j])
	endfor
endfor

; Calculate Volumetric Tensor (Harvey, Ch 12.4, Eq 12.22)

Rvol = dblarr(n_elements(t_master), 3, 3, /nozero)	; Rvol == Volumetric Tensor
for t=0,n_elements(t_master)-1 do begin			; time steps
	rvol_step = dblarr(3,3)		; Could probably be condensed into Rvol[] expression below. Made for easier development troubleshooting
	for sc = 0,3 do begin		; making summation more explicit 
		rvol_step = rvol_step + (reform(rarr[sc,t,*]) # reform(rarr[sc,t,*]))
	endfor
	Rvol[t,*,*] = (1./4.) * (rvol_step - (reform(rm[t,*]) # reform(rm[t,*])))
	; reform() necessary to do matrix math on vectors in time series
	; # operator works on vectors A,B such that A # B = A # transpose(B), removing the need for an explicit transpose
endfor

; Calculate Inverted R tensor for later use
; In future might directly calculate Rinv as I don't ever need Rvol for calculation
Rinv = dblarr(n_elements(t_master), 3, 3, /nozero)
for i=0,n_elements(t_master)-1 do Rinv[i,*,*] = la_invert(reform(Rvol[i,*,*]), /double)


; Calculate grad(b) using both Harvey and Shen methods

grad_Harvey = {x:t_master, y:dblarr(n_elements(t_master), 3, 3, /nozero)}	; initialize grad(b) variable for Harvey, Shen methods
grad_Shen = {x:t_master, y:dblarr(n_elements(t_master), 3, 3, /nozero)}

k_harvey = {x:t_master, y:dblarr(n_elements(t_master), 3, /nozero)}		; initialize curvature vectors for Harvey, Shen methods
k_shen = {x:t_master, y:dblarr(n_elements(t_master), 3, /nozero)}

rc_harvey = {x:t_master, y:dblarr(n_elements(t_master), /nozero)}		; initialize radius of curvature for Harvey, Shen methods
rc_shen = {x:t_master, y:dblarr(n_elements(t_master), /nozero)}

a_ne_b_list = list([1,2,3],[2,3],[3])		; list of all unique combonations of 4 spacecraft for sc_a != sc_b


for t = 0,n_elements(t_master)-1 do begin		; Steps thorugh each time step
	for i =0,2 do begin			; for the final i-component of the gradient
		for j = 0,2 do begin		; for the final j-component of the gradient
			dbdr = dblarr(3)	; a != b summation row vector from Harvey.  Re-initialized at [0,0,0] with each i,j
			for k = 0,2 do begin	; step through the k-index.  May not be needed if can be left for vector operation
				for a = 0,2 do begin	; step through spacecraft MMS1-3; MMS4 covered implicitly
					foreach b, a_ne_b_list[a] do begin	; With previous for loop covers all unique a != b combonaitons
						dbdr[k] = dbdr[k] + ((barr[a,t,i] - barr[b,t,i]) * (rarr[a,t,k] - rarr[b,t,k]))
					endforeach
				endfor
			endfor
			grad_Harvey.y[t,i,j] = (1./16.) * transpose(dbdr) # reform(Rinv[t,*,j])	; original index order
		endfor
	endfor
	; Put in any other time-step calculations here, e.g. b.grad(b)
	k_harvey.y[t,*] = transpose(bm[t,*]) # reform(grad_Harvey.y[t,*,*])
	rc_harvey.y[t] = 1./sqrt(k_harvey.y[t,0]^2 + k_harvey.y[t,1]^2 + k_harvey.y[t,2]^2)
endfor

; Next step: Need to calculate the osculating plane normal 


op_harvey = dblarr(n_elements(t_master), 3)	; initialize variable for osculating plane

for t=0,n_elements(t_master)-1 do begin
	op_harvey[t,*] = crossp(reform(bm[t,*]), reform(k_harvey.y[t,*]))	; Calculate cross product for Osc. Plane
	op_norm = sqrt(op_harvey[t,0]^2 + op_harvey[t,1]^2 + op_harvey[t,2]^2)	; Calculate magnitude of cross product
	op_harvey[t,*] = op_harvey[t,*] * (1./op_norm)		; Normalize Osc. Plane calculation
endfor


; clock angle doesn't work well.  Might be worth looking at when error analysis in place

;rot_clock = dblarr(n_elements(t_master), 4)
;for t=0,n_elements(t_master)-1 do begin
;	rot_clock[t,0] = (180/!pi)*atan(k_harvey.y[t,2], k_harvey.y[t,1])	; clock angle about x_gsm (twisting in tail)
;	rot_clock[t,1] = sqrt(k_harvey.y[t,2]^2 + k_harvey.y[t,1]^2)		; magnitude of the curvature about x_gsm  
;	rot_clock[t,2] = (180/!pi)*atan(k_harvey.y[t,2], k_harvey.y[t,0])	; clock angle about y_gsm (tilt or flap in tail)
;	rot_clock[t,3] = sqrt(k_harvey.y[t,2]^2 + k_harvey.y[t,0]^2)		; magnitude of the curvature about y_gsm
;endfor
;
;store_data, 'b_curvature_clock_all', data={x:t_master, y:rot_clock}
;store_data, 'b_curvature_clock_about_xgsm', data = {x:t_master, y:rot_clock[*,0:1]}
;store_data, 'b_curvature_clock_about_ygsm', data={x:t_master, y:rot_clock[*,2:3]}


; Store all the tplot variables and do a first pass at reasonable options and ylimits

store_data, 'b_osc_plane_normal_gsm', data={x:t_master, y:op_harvey}
options, 'b_osc_plane_normal_gsm', 'labels', ['x_GSM', 'y_GSM', 'z_GSM']
options, 'b_osc_plane_normal_gsm', 'colors', [2, 4, 6]
options, 'b_osc_plane_normal_gsm', 'labflag', -1
options, 'b_osc_plane_normal_gsm', 'databar', 0.

store_data, 'b_curvature_vector_gsm', data=k_harvey
store_data, 'b_curvature_radius', data=rc_harvey
store_data, 'b_gradient_gsm', data=grad_Harvey

options, 'b_curvature_vector_gsm', 'labels', ['x_GSM', 'y_GSM', 'z_GSM']
options, 'b_curvature_vector_gsm', 'colors', [2, 4, 6]
options, 'b_curvature_vector_gsm', 'labflag', -1
options, 'b_curvature_vector_gsm', 'databar', 0.

options, 'b_curvature_radius', 'ytitle', 'R_C (km)'
options, 'b_curvature_radius', 'ylog', 1

if makeplot ne 0 then begin
	tplot, trange=trange, ['b_curvature_vector_gsm', 'b_curvature_radius']
endif



; Thanks are owed to Tim Rogers (virmitio) for useful conversations regarding the row-indexed vs column-indexed nature of programming languages as well as his assistance in cleaning up the markdown of the README.md file.

end
