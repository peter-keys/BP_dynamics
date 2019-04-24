;--------------------------------------------------------------------
;	TRIM INVERTED TEMPORAL SUBFIELDS BACK TO ORIGINAL
;
;--------------------------------------------------------------------
; 02/08/2018
;
; > To speed up the inversions of the data I ran 
;   them as a horizontal stack of images in time.
;   Also, I had to split them in half as the 
;   inverisons wouldn't run otherwise. 
;	_1 is 1st batch
;	_2 is 2nd batch
;
;   I also added in 2 frames at the start/end 
;   as it could be useful.
;   *********************** 
;
;	Therefore this code takes the original 
;	values used to trim the data up (after 
; 	the input of the MBP number) and then 
;	loads in the inverted data and returns 
;	a set of time sereis of inverted data
;	Useful data to be returned is:
;		- Observed I
;		- Synthetic I
;		- Observed V
;		- Synthetic V
;		- Temperature
;		- LOS Vel
;		- B LOS Z
;		- GasP
;		- Rho
;		- ELP
;	At:	
;		- tau = 0
;		- tau = -1
;		- tau = -2
;		- tau = -3
;   *********************** 

PRO return_trimmed, MBP_num

;--------------------------------------
; READ IN THE MBP LOCATION INFO
;--------------------------------------
readcol,'/data/solarstore3/phk/data/27Jul2014/final_tracking_file',f='D,D,D,D,D,x', btpt_num, frame, xloc, yloc, area

;MBP_num = 279.

xcoords = xloc[WHERE (btpt_num eq MBP_num)]  
ycoords = yloc[WHERE (btpt_num eq MBP_num)]   				
frame_num = frame[WHERE (btpt_num eq MBP_num)]  
initial_frame = frame_num[0] 
final_frame = frame_num[N_ELEMENTS(frame_num) - 1]   
n_frames = N_ELEMENTS(frame_num) 

	; > Correct values for trimming	
xcoords = xcoords + 60. 
ycoords = ycoords + 65. 
frame_num = frame_num + 3. 
initial_frame = initial_frame + 3. 
final_frame = final_frame +3. 

	; > Add 2 frames either side of the MBP to get info prior to formation/after 
	; > IFs stop crashing later if MBP appears at start/end of the scans
initial_frame = initial_frame - 2. 
IF (initial_frame LT 0) THEN initial_frame = 0. 
final_frame = final_frame +2 
IF (final_frame GT 91) THEN final_frame = 91. 	

; 27Jul2014 QS 6302 Data:
;file = '/data/solarstore4/vh/SSTreduc/2014.07.27_6302dc_remod_2018/crispex/14:18:38/crispex.stokes.6302.14:18:38.time_corrected.fcube'
file = '/data/solarstore4/vh/SSTreduc/2014.07.27_6302dc_remod_2018/crispex/14:18:38/crispex.stokes.6302.14:18:38.time_corrected_cmapcorr.fcube'
file1 = '/data/solarstore4/vh/SSTreduc/2014.07.27_6302dc_remod_2018/crispex/14:18:38/crispex.stokes.6302.14:18:38.time_corrected_sp.fcube'
restore,'/data/solarstore4/vh/SSTreduc/2014.07.27_6302dc/crispex/14:18:38/wav.idl',/ver
;wav = f0('/data/solarstore4/vh/SSTreduc/2014.07.27_6302dc/crispex/14:18:38/wav.6302.f0')	; f0 requires just IDL (not ssw) and ANA enabled

;;
;; Get dimensions from header
;;
lp_header, file, nx = nx, ny = ny
lp_header, file1,nx = nw, ny = nt


;;
;; read one snapshot only the first time it is executed		<--- Might need to augment this bit to do all scans...
;;
if(n_elements(d) eq 0) then begin &$
   openr, lun, file, /get_lun &$
   dat = assoc(lun, fltarr(nx, ny, nw, 4, /nozero), 512) &$
   ;dat = assoc(lun, intarr(nx, ny, nw, 4, /nozero), 512) &$ 	; Careful with the INT/FLT here - may mess things up
   d = float(dat[initial_frame]) &$ ;; for example snapshot 10 		; USE 30 here as it is the 1st frame of this MBP (174) 
   free_lun, lun &$	
endif
im_dims = SIZE(d)

	; > Calculate the biggest/lowest x and y values to get dimensiosn of sub field to invert
	; > Use IFs to make sure the low and high values aren't outside image dimensions
x0 = MIN(xcoords)-50 
x1 = MAX(xcoords)+50 
y0 = MIN(ycoords)-50 
y1 = MAX(ycoords)+50 	
IF (x0 LT 0) THEN x0 = 0. 
IF (x1 GT im_dims[1]) THEN x1 = im_dims[1] - 1. 
IF (y0 LT 0) THEN y0 = 0. 
IF (y1 GT im_dims[2]) THEN y1 = im_dims[2] - 1. 	

; > Dimensions of your data set
nx1 = x1 - x0 + 1
ny1 = y1 - y0 + 1
nf1 = final_frame - initial_frame + 1

;--------------------------------------
; LOAD IN THE DATA TO BE TRIMMED
;--------------------------------------
mbp_as_str = STRCOMPRESS(STRING(FIX(MBP_num)),/REMOVE_ALL)
IF (MBP_num LT 100) THEN mbp_as_Str = '0'+mbp_as_str
dir = '/data/solarstore3/phk/data/27Jul2014/Inversions/new_2018_results/MBP'+mbp_as_str+'/'

; > Set up arrays as:
;	x,y,scan,stokes,time
; > ONly have 2 stokes as Q + U ignored..

obs = FLTARR(nx1,ny1,92,2,nf1)
syn = FLTARR(nx1,ny1,92,2,nf1)
syn_r1 = FLTARR(nx1,ny1,92,2,nf1)

; > Read in your profiles
oi1 = READ_PROFILE(dir+'observed_MBP'+STRCOMPRESS(STRING(FIX(MBP_num)),/REMOVE_ALL)+'_1.nic',oq1,ou1,ov1)
dims1 = SIZE(oi1)
oi2 = READ_PROFILE(dir+'observed_MBP'+STRCOMPRESS(STRING(FIX(MBP_num)),/REMOVE_ALL)+'_2.nic',oq2,ou2,ov2)
dims2 = SIZE(oi2)
si1 = READ_PROFILE(dir+'synthetic_MBP'+STRCOMPRESS(STRING(FIX(MBP_num)),/REMOVE_ALL)+'_1.pro',sq1,su1,sv1)
si2 = READ_PROFILE(dir+'synthetic_MBP'+STRCOMPRESS(STRING(FIX(MBP_num)),/REMOVE_ALL)+'_2.pro',sq2,su2,sv2)
si1_r1 = READ_PROFILE(dir+'r1synthetic_MBP'+STRCOMPRESS(STRING(FIX(MBP_num)),/REMOVE_ALL)+'_1.pro',sq1_r1,su1_r1,sv1_r1)
si2_r1 = READ_PROFILE(dir+'r1synthetic_MBP'+STRCOMPRESS(STRING(FIX(MBP_num)),/REMOVE_ALL)+'_2.pro',sq2_r1,su2_r1,sv2_r1)


; > Number of frames in your data...
frames1 = N_ELEMENTS(oi1[*,0,0])/nx1
frames2 = N_ELEMENTS(oi2[*,0,0])/nx1

IF (frames1+frames2) NE nf1 THEN print,'Your frames dont match... Something is up'

; > Find the divisions for each frame in x
;divx1 = dims1[1]/nx1


; > Loop for frist half of data...
FOR t = 0, frames1-1 DO BEGIN &$
	obs[*,*,*,0,t] = oi1[t*nx1:((t+1)*nx1)-1,*,*] &$
	obs[*,*,*,1,t] = ov1[t*nx1:((t+1)*nx1)-1,*,*] &$
	syn[*,*,*,0,t] = si1[t*nx1:((t+1)*nx1)-1,*,*] &$
	syn[*,*,*,1,t] = sv1[t*nx1:((t+1)*nx1)-1,*,*] &$
	syn_r1[*,*,*,0,t] = si1_r1[t*nx1:((t+1)*nx1)-1,*,*] &$
	syn_r1[*,*,*,1,t] = sv1_r1[t*nx1:((t+1)*nx1)-1,*,*] &$	
	FOR x = 0, dims1[1]-1 DO BEGIN &$
	ENDFOR &$
	print,'Running batch 1, frame ',t,' of ',nf1 &$
ENDFOR
; > Loop for 2nd half of data...
FOR t = frames1, frames1+(frames2)-1 DO BEGIN &$
	obs[*,*,*,0,t] = oi2[(t-frames1)*nx1:(((t-frames1)+1)*nx1)-1,*,*] &$
	obs[*,*,*,1,t] = ov2[(t-frames1)*nx1:(((t-frames1)+1)*nx1)-1,*,*] &$
	syn[*,*,*,0,t] = si2[(t-frames1)*nx1:(((t-frames1)+1)*nx1)-1,*,*] &$
	syn[*,*,*,1,t] = sv2[(t-frames1)*nx1:(((t-frames1)+1)*nx1)-1,*,*] &$
	syn_r1[*,*,*,0,t] = si2_r1[(t-frames1)*nx1:(((t-frames1)+1)*nx1)-1,*,*] &$
	syn_r1[*,*,*,1,t] = sv2_r1[(t-frames1)*nx1:(((t-frames1)+1)*nx1)-1,*,*] &$	
	print,'Running batch 2, frame ',t,' of ',nf1 &$
ENDFOR

info_obs = 'Array is [x,y,lambda,stokes(I&V only),scan_number]'
SAVE,FILENAME=dir+'MBP'+STRCOMPRESS(STRING(FIX(MBP_num)),/REMOVE_ALL)+'_observed_timeseries.sav',obs,info_obs
SAVE,FILENAME=dir+'MBP'+STRCOMPRESS(STRING(FIX(MBP_num)),/REMOVE_ALL)+'_synthetic_timeseries.sav',syn,info_obs
SAVE,FILENAME=dir+'MBP'+STRCOMPRESS(STRING(FIX(MBP_num)),/REMOVE_ALL)+'_synthetic_r1_timeseries.sav',syn_r1,info_obs


;--------------------------------------
; FIND THE ELEMENT CORRESPONDING TO 
; SPECIFIC TAU VALUES IN THE MODELS
;--------------------------------------
; > Want to find the point in the model that is closest to the tau value we want

;> Load in the models and make tau a single dimension array
m = READ_MODEL(dir+'modelout_MBP'+STRCOMPRESS(STRING(FIX(MBP_num)),/REMOVE_ALL)+'_1.mod')
tau1 = reform(m.tau[0,0,*]) 
m2 = READ_MODEL(dir+'modelout_MBP'+STRCOMPRESS(STRING(FIX(MBP_num)),/REMOVE_ALL)+'_2.mod')
tau2 = reform(m2.tau[0,0,*]) 

; > Set the tau valeus you want
zero = 0.
mone = -1.
mtwo = -2.
mthree = -3.

;> Work out the element in the tau array that is closest to these values
;  Here 'near' = the nearest value to the value you're searching for
;       'index' is the number in the array that corresponds to the nearest value
neart0 = MIN(ABS(tau1 - zero),indext0)
neartm1 = MIN(ABS(tau1 - mone),indextm1)
neartm2 = MIN(ABS(tau1 - mtwo),indextm2)
neartm3 = MIN(ABS(tau1 - mthree),indextm3)

; > Now create the relvant models

Bfield = FLTARR(nx1,ny1,4,nf1)
Temp = FLTARR(nx1,ny1,4,nf1)
GasP = FLTARR(nx1,ny1,4,nf1)
rho = FLTARR(nx1,ny1,4,nf1)
ElP = FLTARR(nx1,ny1,4,nf1)
LOSv = FLTARR(nx1,ny1,4,nf1)

Bfield[*] = 0.
Temp[*] = 0.
GasP[*] = 0.
rho[*] = 0.
ElP[*] = 0.
LOSv[*] = 0.

info_mod = 'Array is [x,y,t,tau] while 0:3 in tau is 0, -1, -2, -3, respectively'

FOR t = 0, frames1-1 DO BEGIN &$
	Bfield[*,*,0,t] = m.B_LOS_z[t*nx1:((t+1)*nx1)-1,*,indext0] &$
	Bfield[*,*,1,t] = m.B_LOS_z[t*nx1:((t+1)*nx1)-1,*,indextm1] &$
	Bfield[*,*,2,t] = m.B_LOS_z[t*nx1:((t+1)*nx1)-1,*,indextm2] &$
	Bfield[*,*,3,t] = m.B_LOS_z[t*nx1:((t+1)*nx1)-1,*,indextm3] &$

	Temp[*,*,0,t] = m.T[t*nx1:((t+1)*nx1)-1,*,indext0] &$
	Temp[*,*,1,t] = m.T[t*nx1:((t+1)*nx1)-1,*,indextm1] &$
	Temp[*,*,2,t] = m.T[t*nx1:((t+1)*nx1)-1,*,indextm2] &$
	Temp[*,*,3,t] = m.T[t*nx1:((t+1)*nx1)-1,*,indextm3] &$

	GasP[*,*,0,t] = m.GAS_P[t*nx1:((t+1)*nx1)-1,*,indext0] &$
	GasP[*,*,1,t] = m.GAS_P[t*nx1:((t+1)*nx1)-1,*,indextm1] &$
	GasP[*,*,2,t] = m.GAS_P[t*nx1:((t+1)*nx1)-1,*,indextm2] &$
	GasP[*,*,3,t] = m.GAS_P[t*nx1:((t+1)*nx1)-1,*,indextm3] &$

	Rho[*,*,0,t] = m.RHO[t*nx1:((t+1)*nx1)-1,*,indext0] &$
	Rho[*,*,1,t] = m.RHO[t*nx1:((t+1)*nx1)-1,*,indextm1] &$
	Rho[*,*,2,t] = m.RHO[t*nx1:((t+1)*nx1)-1,*,indextm2] &$
	Rho[*,*,3,t] = m.RHO[t*nx1:((t+1)*nx1)-1,*,indextm3] &$

	ELP[*,*,0,t] = m.EL_P[t*nx1:((t+1)*nx1)-1,*,indext0] &$
	ELP[*,*,1,t] = m.EL_P[t*nx1:((t+1)*nx1)-1,*,indextm1] &$
	ELP[*,*,2,t] = m.EL_P[t*nx1:((t+1)*nx1)-1,*,indextm2] &$
	ELP[*,*,3,t] = m.EL_P[t*nx1:((t+1)*nx1)-1,*,indextm3] &$

	LOSv[*,*,0,t] = m.V_LOS[t*nx1:((t+1)*nx1)-1,*,indext0] &$
	LOSv[*,*,1,t] = m.V_LOS[t*nx1:((t+1)*nx1)-1,*,indextm1] &$
	LOSv[*,*,2,t] = m.V_LOS[t*nx1:((t+1)*nx1)-1,*,indextm2] &$
	LOSv[*,*,3,t] = m.V_LOS[t*nx1:((t+1)*nx1)-1,*,indextm3] &$

	print,'Running Model batch 1, frame ',t,' of ',nf1 &$
ENDFOR
; > Loop for 2nd half of data...
FOR t = frames1, frames1+(frames2)-1 DO BEGIN &$	
	Bfield[*,*,0,t] = m2.B_LOS_z[(t-frames1)*nx1:(((t-frames1)+1)*nx1)-1,*,indext0] &$
	Bfield[*,*,1,t] = m2.B_LOS_z[(t-frames1)*nx1:(((t-frames1)+1)*nx1)-1,*,indextm1] &$
	Bfield[*,*,2,t] = m2.B_LOS_z[(t-frames1)*nx1:(((t-frames1)+1)*nx1)-1,*,indextm2] &$
	Bfield[*,*,3,t] = m2.B_LOS_z[(t-frames1)*nx1:(((t-frames1)+1)*nx1)-1,*,indextm3] &$

	Temp[*,*,0,t] = m2.T[(t-frames1)*nx1:(((t-frames1)+1)*nx1)-1,*,indext0] &$
	Temp[*,*,1,t] = m2.T[(t-frames1)*nx1:(((t-frames1)+1)*nx1)-1,*,indextm1] &$
	Temp[*,*,2,t] = m2.T[(t-frames1)*nx1:(((t-frames1)+1)*nx1)-1,*,indextm2] &$
	Temp[*,*,3,t] = m2.T[(t-frames1)*nx1:(((t-frames1)+1)*nx1)-1,*,indextm3] &$

	GasP[*,*,0,t] = m2.GAS_P[(t-frames1)*nx1:(((t-frames1)+1)*nx1)-1,*,indext0] &$
	GasP[*,*,1,t] = m2.GAS_P[(t-frames1)*nx1:(((t-frames1)+1)*nx1)-1,*,indextm1] &$
	GasP[*,*,2,t] = m2.GAS_P[(t-frames1)*nx1:(((t-frames1)+1)*nx1)-1,*,indextm2] &$
	GasP[*,*,3,t] = m2.GAS_P[(t-frames1)*nx1:(((t-frames1)+1)*nx1)-1,*,indextm3] &$

	Rho[*,*,0,t] = m2.RHO[(t-frames1)*nx1:(((t-frames1)+1)*nx1)-1,*,indext0] &$
	Rho[*,*,1,t] = m2.RHO[(t-frames1)*nx1:(((t-frames1)+1)*nx1)-1,*,indextm1] &$
	Rho[*,*,2,t] = m2.RHO[(t-frames1)*nx1:(((t-frames1)+1)*nx1)-1,*,indextm2] &$
	Rho[*,*,3,t] = m2.RHO[(t-frames1)*nx1:(((t-frames1)+1)*nx1)-1,*,indextm3] &$

	ELP[*,*,0,t] = m2.EL_P[(t-frames1)*nx1:(((t-frames1)+1)*nx1)-1,*,indext0] &$
	ELP[*,*,1,t] = m2.EL_P[(t-frames1)*nx1:(((t-frames1)+1)*nx1)-1,*,indextm1] &$
	ELP[*,*,2,t] = m2.EL_P[(t-frames1)*nx1:(((t-frames1)+1)*nx1)-1,*,indextm2] &$
	ELP[*,*,3,t] = m2.EL_P[(t-frames1)*nx1:(((t-frames1)+1)*nx1)-1,*,indextm3] &$

	LOSv[*,*,0,t] = m2.V_LOS[(t-frames1)*nx1:(((t-frames1)+1)*nx1)-1,*,indext0] &$
	LOSv[*,*,1,t] = m2.V_LOS[(t-frames1)*nx1:(((t-frames1)+1)*nx1)-1,*,indextm1] &$
	LOSv[*,*,2,t] = m2.V_LOS[(t-frames1)*nx1:(((t-frames1)+1)*nx1)-1,*,indextm2] &$
	LOSv[*,*,3,t] = m2.V_LOS[(t-frames1)*nx1:(((t-frames1)+1)*nx1)-1,*,indextm3] &$
	
	print,'Running Model batch 2, frame ',t,' of ',nf1 &$
ENDFOR

SAVE,FILENAME=dir+'MBP'+STRCOMPRESS(STRING(FIX(MBP_num)),/REMOVE_ALL)+'_model_timeseries.sav',Bfield,LOSv,Temp,GasP,Rho,ELP,info_mod

; > Do the same for r1 as it sometimes looks better for talks etc...

m = READ_MODEL(dir+'r1modelout_MBP'+STRCOMPRESS(STRING(FIX(MBP_num)),/REMOVE_ALL)+'_1.mod')
m2 = READ_MODEL(dir+'r1modelout_MBP'+STRCOMPRESS(STRING(FIX(MBP_num)),/REMOVE_ALL)+'_2.mod')


Bfield_r1 = FLTARR(nx1,ny1,4,nf1)
Temp_r1 = FLTARR(nx1,ny1,4,nf1)
GasP_r1 = FLTARR(nx1,ny1,4,nf1)
rho_r1 = FLTARR(nx1,ny1,4,nf1)
ElP_r1 = FLTARR(nx1,ny1,4,nf1)
LOSv_r1 = FLTARR(nx1,ny1,4,nf1)

Bfield_r1[*] = 0.
Temp_r1[*] = 0.
GasP_r1[*] = 0.
rho_r1[*] = 0.
ElP_r1[*] = 0.
LOSv_r1[*] = 0.

info_mod = 'Array is [x,y,t,tau] while 0:3 in tau is 0, -1, -2, -3, respectively'

FOR t = 0, frames1-1 DO BEGIN &$
	Bfield_r1[*,*,0,t] = m.B_LOS_z[t*nx1:((t+1)*nx1)-1,*,indext0] &$
	Bfield_r1[*,*,1,t] = m.B_LOS_z[t*nx1:((t+1)*nx1)-1,*,indextm1] &$
	Bfield_r1[*,*,2,t] = m.B_LOS_z[t*nx1:((t+1)*nx1)-1,*,indextm2] &$
	Bfield_r1[*,*,3,t] = m.B_LOS_z[t*nx1:((t+1)*nx1)-1,*,indextm3] &$

	Temp_r1[*,*,0,t] = m.T[t*nx1:((t+1)*nx1)-1,*,indext0] &$
	Temp_r1[*,*,1,t] = m.T[t*nx1:((t+1)*nx1)-1,*,indextm1] &$
	Temp_r1[*,*,2,t] = m.T[t*nx1:((t+1)*nx1)-1,*,indextm2] &$
	Temp_r1[*,*,3,t] = m.T[t*nx1:((t+1)*nx1)-1,*,indextm3] &$

	GasP_r1[*,*,0,t] = m.GAS_P[t*nx1:((t+1)*nx1)-1,*,indext0] &$
	GasP_r1[*,*,1,t] = m.GAS_P[t*nx1:((t+1)*nx1)-1,*,indextm1] &$
	GasP_r1[*,*,2,t] = m.GAS_P[t*nx1:((t+1)*nx1)-1,*,indextm2] &$
	GasP_r1[*,*,3,t] = m.GAS_P[t*nx1:((t+1)*nx1)-1,*,indextm3] &$

	Rho_r1[*,*,0,t] = m.RHO[t*nx1:((t+1)*nx1)-1,*,indext0] &$
	Rho_r1[*,*,1,t] = m.RHO[t*nx1:((t+1)*nx1)-1,*,indextm1] &$
	Rho_r1[*,*,2,t] = m.RHO[t*nx1:((t+1)*nx1)-1,*,indextm2] &$
	Rho_r1[*,*,3,t] = m.RHO[t*nx1:((t+1)*nx1)-1,*,indextm3] &$

	ELP_r1[*,*,0,t] = m.EL_P[t*nx1:((t+1)*nx1)-1,*,indext0] &$
	ELP_r1[*,*,1,t] = m.EL_P[t*nx1:((t+1)*nx1)-1,*,indextm1] &$
	ELP_r1[*,*,2,t] = m.EL_P[t*nx1:((t+1)*nx1)-1,*,indextm2] &$
	ELP_r1[*,*,3,t] = m.EL_P[t*nx1:((t+1)*nx1)-1,*,indextm3] &$

	LOSv_r1[*,*,0,t] = m.V_LOS[t*nx1:((t+1)*nx1)-1,*,indext0] &$
	LOSv_r1[*,*,1,t] = m.V_LOS[t*nx1:((t+1)*nx1)-1,*,indextm1] &$
	LOSv_r1[*,*,2,t] = m.V_LOS[t*nx1:((t+1)*nx1)-1,*,indextm2] &$
	LOSv_r1[*,*,3,t] = m.V_LOS[t*nx1:((t+1)*nx1)-1,*,indextm3] &$

	print,'Running Model batch 1, frame ',t,' of ',nf1 &$
ENDFOR
; > Loop for 2nd half of data...
FOR t = frames1, frames1+(frames2)-1 DO BEGIN &$	
	Bfield_r1[*,*,0,t] = m2.B_LOS_z[(t-frames1)*nx1:(((t-frames1)+1)*nx1)-1,*,indext0] &$
	Bfield_r1[*,*,1,t] = m2.B_LOS_z[(t-frames1)*nx1:(((t-frames1)+1)*nx1)-1,*,indextm1] &$
	Bfield_r1[*,*,2,t] = m2.B_LOS_z[(t-frames1)*nx1:(((t-frames1)+1)*nx1)-1,*,indextm2] &$
	Bfield_r1[*,*,3,t] = m2.B_LOS_z[(t-frames1)*nx1:(((t-frames1)+1)*nx1)-1,*,indextm3] &$

	Temp_r1[*,*,0,t] = m2.T[(t-frames1)*nx1:(((t-frames1)+1)*nx1)-1,*,indext0] &$
	Temp_r1[*,*,1,t] = m2.T[(t-frames1)*nx1:(((t-frames1)+1)*nx1)-1,*,indextm1] &$
	Temp_r1[*,*,2,t] = m2.T[(t-frames1)*nx1:(((t-frames1)+1)*nx1)-1,*,indextm2] &$
	Temp_r1[*,*,3,t] = m2.T[(t-frames1)*nx1:(((t-frames1)+1)*nx1)-1,*,indextm3] &$

	GasP_r1[*,*,0,t] = m2.GAS_P[(t-frames1)*nx1:(((t-frames1)+1)*nx1)-1,*,indext0] &$
	GasP_r1[*,*,1,t] = m2.GAS_P[(t-frames1)*nx1:(((t-frames1)+1)*nx1)-1,*,indextm1] &$
	GasP_r1[*,*,2,t] = m2.GAS_P[(t-frames1)*nx1:(((t-frames1)+1)*nx1)-1,*,indextm2] &$
	GasP_r1[*,*,3,t] = m2.GAS_P[(t-frames1)*nx1:(((t-frames1)+1)*nx1)-1,*,indextm3] &$

	Rho_r1[*,*,0,t] = m2.RHO[(t-frames1)*nx1:(((t-frames1)+1)*nx1)-1,*,indext0] &$
	Rho_r1[*,*,1,t] = m2.RHO[(t-frames1)*nx1:(((t-frames1)+1)*nx1)-1,*,indextm1] &$
	Rho_r1[*,*,2,t] = m2.RHO[(t-frames1)*nx1:(((t-frames1)+1)*nx1)-1,*,indextm2] &$
	Rho_r1[*,*,3,t] = m2.RHO[(t-frames1)*nx1:(((t-frames1)+1)*nx1)-1,*,indextm3] &$

	ELP_r1[*,*,0,t] = m2.EL_P[(t-frames1)*nx1:(((t-frames1)+1)*nx1)-1,*,indext0] &$
	ELP_r1[*,*,1,t] = m2.EL_P[(t-frames1)*nx1:(((t-frames1)+1)*nx1)-1,*,indextm1] &$
	ELP_r1[*,*,2,t] = m2.EL_P[(t-frames1)*nx1:(((t-frames1)+1)*nx1)-1,*,indextm2] &$
	ELP_r1[*,*,3,t] = m2.EL_P[(t-frames1)*nx1:(((t-frames1)+1)*nx1)-1,*,indextm3] &$

	LOSv_r1[*,*,0,t] = m2.V_LOS[(t-frames1)*nx1:(((t-frames1)+1)*nx1)-1,*,indext0] &$
	LOSv_r1[*,*,1,t] = m2.V_LOS[(t-frames1)*nx1:(((t-frames1)+1)*nx1)-1,*,indextm1] &$
	LOSv_r1[*,*,2,t] = m2.V_LOS[(t-frames1)*nx1:(((t-frames1)+1)*nx1)-1,*,indextm2] &$
	LOSv_r1[*,*,3,t] = m2.V_LOS[(t-frames1)*nx1:(((t-frames1)+1)*nx1)-1,*,indextm3] &$
	
	print,'Running Model batch 2, frame ',t,' of ',nf1 &$
ENDFOR

SAVE,FILENAME=dir+'MBP'+STRCOMPRESS(STRING(FIX(MBP_num)),/REMOVE_ALL)+'_model_r1_timeseries.sav',Bfield_r1,LOSv_r1,Temp_r1,GasP_r1,Rho_r1,ELP_r1,info_mod


;--------------------------------------
; OUTPUT A SUBFIELD OF THE TRACKING 
;  FILE TO OVERPLOT ETC..
;--------------------------------------
restore,'/data/solarstore3/phk/data/27Jul2014/final_tracked_mbps.sav',/ver

; >Reload the xcoords etc to get the more accurate x0,y0 vals...

readcol,'/data/solarstore3/phk/data/27Jul2014/final_tracking_file',f='D,D,D,D,D,x', btpt_num, frame, xloc, yloc, area

;MBP_num = 35.

xcoords = xloc[WHERE (btpt_num eq MBP_num)]  
ycoords = yloc[WHERE (btpt_num eq MBP_num)]   				
frame_num = frame[WHERE (btpt_num eq MBP_num)]  
initial_frame = frame_num[0] 
final_frame = frame_num[N_ELEMENTS(frame_num) - 1]   
n_frames = N_ELEMENTS(frame_num) 
dims = SIZE(MBPs_new)

	; > Calculate the biggest/lowest x and y values to get dimensiosn of sub field to invert
	; > Use IFs to make sure the low and high values aren't outside image dimensions
x0 = MIN(xcoords)-50 
x1 = MAX(xcoords)+50 
y0 = MIN(ycoords)-50 
y1 = MAX(ycoords)+50 	
IF (x0 LT 0) THEN x0 = 0. 
IF (x1 GT dims[1]) THEN x1 = dims[1] - 1. 
IF (y0 LT 0) THEN y0 = 0. 
IF (y1 GT dims[2]) THEN y1 = dims[2] - 1. 	

mbp_t = MBPs_new[x0:x1,y0:y1,initial_frame:final_frame]

MBP_track = FLTARR(nx1,ny1,nf1)
MBP_track[*,*,2:nf1-3] = mbp_t

SAVE,FILENAME=dir+'MBP'+STRCOMPRESS(STRING(FIX(MBP_num)),/REMOVE_ALL)+'_tracks.sav',MBP_track

;----------------------------------------
; TRIMMIING THE TRACKS
;
; Maybe it's a good idea to trim the tracks 
; down to only account for the MBP by itself.
; Can then use the trimmed tracks to work out 
; some info from the models.
;
; This is outside th escope of automation 
; as it requires a human eye to remove things 
;
;----------------------------------------
;
;mbp_tracksolo = MBP_track
;
;tvim,MBP_track[*,*,2+0]
tvim,MBP_tracksolo[*,*,2+6]
verline,xcoords[33]-x0
horline,ycoords[33]-y0

;
; MBP 35:
; (Did 2 here, one that removes the split after frame 18 and one that leaves it as is)
;mbp_tracksolo[*,0:34,*] = 0.
;mbp_tracksolo[*,100:*,*] =0.
;mbp_tracksolo[*,61:*,2] = 0.
;mbp_tracksolo[0:58,0:49,7] = 0.
;mbp_tracksolo[0:58,0:50,8] = 0.
;mbp_tracksolo[*,100:*,15] = 0.
;mbp_tracksolo[57:*,60:*,16] = 0.
;mbp_tracksolo[57:*,60:*,17] = 0.
;SAVE,FILENAME=dir+'MBP'+STRCOMPRESS(STRING(FIX(MBP_num)),/REMOVE_ALL)+'_tracks_with_split.sav',MBP_tracksolo
;mbp_tracksolo[52:*,55:*,18] = 0.
;mbp_tracksolo[52:*,55:*,19] = 0.
;mbp_tracksolo[52:*,55:*,20] = 0.
;mbp_tracksolo[52:*,55:*,18:*] = 0.
;SAVE,FILENAME=dir+'MBP'+STRCOMPRESS(STRING(FIX(MBP_num)),/REMOVE_ALL)+'_tracks_solo.sav',MBP_tracksolo

; MBP 174:
; (Did 2 here, one that removes the split after frame 18 and one that leaves it as is)
;mbp_tracksolo[*,0:45,*] = 0.
;mbp_tracksolo[*,80:*,*] = 0.
;mbp_tracksolo[0:40,*,*] = 0.
;mbp_tracksolo[52:57,49:53,3] = 0.
;mbp_tracksolo[52:57,47:52,4] = 0.
;mbp_tracksolo[51:63,46:52,5] = 0.
;mbp_tracksolo[51:56,45:51,6] = 0.
;mbp_tracksolo[*,0:48,11:20] = 0.
;mbp_tracksolo[66:*,*,27:32] = 0.
;mbp_tracksolo[*,0:50,33:*] = 0.

;----------------------------------------

END

; Quick test of MBP Bfield / LOSv properties over whole MBP

;new_field = FLTARR(nx1,ny1,nf1)
;FOR t = 0, nf1-1 DO new_field[*,*,t] = Bfield[*,*,0,t] * mbp_tracksolo[*,*,t]

;B= DBLARR(nf1) 
;FOR t = 0, nf1-1 DO B[t] = MEAN(WHERE(ABS(new_field[*,*,t]) GT 0)) 

;END
