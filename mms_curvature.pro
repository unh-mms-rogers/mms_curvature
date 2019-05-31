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
b1 = {x:bt1.x, y:fltarr(n_elements(bt1.x), 3, /nozero)}
for i=0,n_elements(bt1.x)-1 do begin			; Step through all time steps in B data
	for j=0,2 do begin				; Step through all 3 components in each time step
		b1.y[i,j] = bt1.y[i,j]/bt1.y[i,3]	; divide value of component by magnitude
	endfor
endfor

b2 = {x:bt2.x, y:fltarr(n_elements(bt2.x), 3, /nozero)}
for i=0,n_elements(bt2.x)-1 do begin			; Step through all time steps in B data
	for j=0,2 do begin				; Step through all 3 components in each time step
		b2.y[i,j] = bt2.y[i,j]/bt2.y[i,3]	; divide value of component by magnitude
	endfor
endfor

b3 = {x:bt3.x, y:fltarr(n_elements(bt3.x), 3, /nozero)}
for i=0,n_elements(bt3.x)-1 do begin			; Step through all time steps in B data
	for j=0,2 do begin				; Step through all 3 components in each time step
		b3.y[i,j] = bt3.y[i,j]/bt3.y[i,3]	; divide value of component by magnitude
	endfor
endfor

b4 = {x:bt4.x, y:fltarr(n_elements(bt4.x), 3, /nozero)}
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
		while b1.x[t_end_index] le t_end do t_end_index += 1
		t_master = b1.x[0:t_end_index]
	end
	1: begin
		while b2.x[t_end_index] le t_end do t_end_index += 1
		t_master = b2.x[0:t_end_index]
	end
	2: begin
		while b3.x[t_end_index] le t_end do t_end_index += 1
		t_master = b3.x[0:t_end_index]
	end
	3: begin
		while b4.x[t_end_index] le t_end do t_end_index += 1
		t_master = b4.x[0:t_end_index]
	end
endcase

; interpolated b timeseries, linear interpolation used
bi1 = {x:t_master, y:[[interpol(b1.y[*,0], b1.x, t_master)], [interpol(b1.y[*,1], b1.x, t_master)], [interpol(b1.y[*,2], b1.x, t_master)]]}
bi2 = {x:t_master, y:[[interpol(b2.y[*,0], b2.x, t_master)], [interpol(b2.y[*,1], b2.x, t_master)], [interpol(b2.y[*,2], b2.x, t_master)]]}
bi3 = {x:t_master, y:[[interpol(b3.y[*,0], b3.x, t_master)], [interpol(b3.y[*,1], b3.x, t_master)], [interpol(b3.y[*,2], b3.x, t_master)]]}
bi4 = {x:t_master, y:[[interpol(b4.y[*,0], b4.x, t_master)], [interpol(b4.y[*,1], b4.x, t_master)], [interpol(b4.y[*,2], b4.x, t_master)]]}

; Interpolate position data to increase cadance to that of the FGM data.  Linear interpolation used.
r1 = {x:t_master, y:[[interpol(rr1.y[*,0], rr1.x, t_master)], [interpol(rr1.y[*,1], rr1.x, t_master)], [interpol(rr1.y[*,2], rr1.x, t_master)]]}
r2 = {x:t_master, y:[[interpol(rr2.y[*,0], rr2.x, t_master)], [interpol(rr1.y[*,1], rr2.x, t_master)], [interpol(rr2.y[*,2], rr2.x, t_master)]]}
r3 = {x:t_master, y:[[interpol(rr3.y[*,0], rr3.x, t_master)], [interpol(rr1.y[*,1], rr3.x, t_master)], [interpol(rr3.y[*,2], rr3.x, t_master)]]}
r4 = {x:t_master, y:[[interpol(rr4.y[*,0], rr4.x, t_master)], [interpol(rr1.y[*,1], rr4.x, t_master)], [interpol(rr4.y[*,2], rr4.x, t_master)]]}

; Now all of the magnetic field and positional data are of the same cadance and at the same times for each time index

;t_master = b1.x	; temp until interpolation fleshed out

; Calculate position and normalized magnetic field at mesocenter of fleet
rm = {x:t_master, y:fltarr(n_elements(t_master), 3, /nozero)}
for i=0,n_elements(t_master)-1 do begin
	for j=0,2 do begin
		rm.y[i,j] = (1/4)*total([r1.y[i,j], r2.y[i,j], r3.y[i,j], r4.y[i,j]] )
	endfor
endfor

bm = {x:t_master, y:fltarr(n_elements(t_master), 3, /nozero)}
for i=0,n_elements(t_master)-1 do begin
	for j=0,2 do begin
		bm.y[i,j] = (1/4)*total([bi1.y[i,j], bi2.y[i,j], bi3.y[i,j], bi4.y[i,j]] )
	endfor
endfor

; Calculate Volumetric Tensor (Harvey, Ch 12.4, Eq 12.23)

Rvol = {x:t_master, y:fltarr(n_elements(t_master), 3, 3, /nozero)}
for i=0,n_elements(t_master)-1 do begin
	Rvol.y[i,*,*] = (1/4)*total([(reform(r3.y[i,*]) # reform(r3.y[i,*])), (reform(r3.y[i,*]) # reform(r3.y[i,*])), (reform(r3.y[i,*]) # reform(r3.y[i,*])), (reform(r4.y[i,*]) # reform(r4.y[i,*]))]) - (reform(rm.y[i,*]) # reform(rm.y[i,*]))	
	; reform() necessary to do matrix math on vectors in time series
	; # operator works on vectors A,B such that A # B = A # transpose(B), removing the need for an explicit transpose
endfor

Rinv = {x:t_master, y:fltarr(n_elements(t_master), 3, 3, /nozero)}
for i=0,n_elements(t_master)-1 do Rinv.y[i,*,*] = la_invert(reform(Rvol.y[i,*,*]))
; Inverted R for later use


; Calculate grad(b) using both Harvey and Shen methods
grad_Harvey = {x:t_master, y:dblarr(n_elements(t_master), 3, 3, /nozero)
grad_Shen = {x:t_master, y:dblarr(n_elements(t_master), 3, 3, /nozero)

; need to solve the implicit summation over k as well as every bloody other thing.  Might be able to do this operation vectorized over all of the time axis, but might not be able to.  In any case, this needs to get ironed out for both the Harvey and Shen methods for grad(b)





for i=0,2 do begin
	for j=0,2 do begin
		for k = 0,2 do begin
			grad_Harvey.y[*,i,j] = (bi1.y[*,i] - bi2.y[*,i])*(r1.y[*,k]
			grad_Shen.y[i,j,k] = 
		endfor
	endfor
endfor

; Calculate b curvature vectors and radii, as well as normal of the osculating plane using Harvey and Shen methods


end
