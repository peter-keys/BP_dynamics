;-----------------------------------------------------------------------------
; 			LCT FULL FoV SST IMAGES
;-----------------------------------------------------------------------------
; 	- Want to run LCT on the full FoV as it might give a better 
;	  indication of the flow directions for the granules
;	- Have to think how I want to deal with it 
;	
;	Do I just run it on all the frames then take time averages 
;	near the MBPs existence to see flow directions?
;
;	Do I look at diff images near the MBPs peaks?
;
;	Do I look for a general flow about that time?
;
;	- Need LCT results first before I can answer these questions
;	- Can maybe add a bit at the end to extract the relevant info 
;	  for a particular MBP.

;-------------------------------------------------------
; 	LOAD IN THE DATA...
;-------------------------------------------------------
; 27Jul2014 QS 6302 Data:
;file = '/data/solarstore4/vh/SSTreduc/2014.07.27_6302dc_remod_2018/crispex/14:18:38/crispex.stokes.6302.14:18:38.time_corrected.fcube'
file = '/data/solarstore4/vh/SSTreduc/2014.07.27_6302dc_remod_2018/crispex/14:18:38/crispex.stokes.6302.14:18:38.time_corrected_cmapcorr.fcube'
file1 = '/data/solarstore4/vh/SSTreduc/2014.07.27_6302dc_remod_2018/crispex/14:18:38/crispex.stokes.6302.14:18:38.time_corrected_sp.fcube'
restore,'/data/solarstore4/vh/SSTreduc/2014.07.27_6302dc/crispex/14:18:38/wav.idl',/ver

;;
;; Get dimensions from header
;;
lp_header, file, nx = nx, ny = ny
lp_header, file1,nx = nw, ny = nt

FeImages = FLTARR(nx,ny,nt)

FOR t = 0, nt-1 DO BEGIN &$
	openr, lun, file, /get_lun  &$
		dat = assoc(lun, fltarr(nx, ny, nw, 4, /nozero), 512)  &$
		d = float(dat[t]) &$ ;; for example snapshot 10 		; USE 30 here as it is the 1st frame of this MBP (174) 
	free_lun, lun  &$
	FeImages[*,*,t] =  d[*,*,16,0] &$
ENDFOR
FeImages = FeImages[60:912,65:940,3:*]
restore,'/data/solarstore3/phk/data/27Jul2014/wide_band_data.sav',/ver

CaImages = dataset
;delvar,dataset

; > Normalise them to the Median so that I don't have as many issues with massive differences 
;   between frames for the difference imaging later on....

FOR t = 0, N_ELEMENTS(FeImages[0,0,*])-1 DO FeImages[*,*,t] = FeImages[*,*,t]/MEDIAN(FeImages[*,*,t])
;FOR t = 0, N_ELEMENTS(CaImages[0,0,*])-1 DO CaImages[*,*,t] = CaImages[*,*,t]/MEDIAN(CaImages[*,*,t])

; > Possible solution to the intensity threshold issue related to tracking in Ca/Fe....
;t = (caImages[*,*,0]-FeImages[*,*,0])/MEDIAN((caImages[*,*,0]-FeImages[*,*,0]))

;-------------------------------------------------------
;   RUN LCT ON BOTH AND SAVE OUT THE RESULTS
;-------------------------------------------------------
; > Run as Full temporal res, then apodised temporal res.

;-------------------------
; 	FE IMAGES 1st
;-------------------------
; > Image res and the size of a granule (for spatial window)
scale = 50.
granule = 1500.

; > Time scale of scan and the half scan (for temporal windowing)
time = 33.
time_2 = time/2.

; > Lifetime of a granule (to work out the number of images needed)
granule_time = (60*5.)
half_granule_t = granule_time/2.

sizing = SIZE(FeImages)

; > Move onto the LCT of this nonsense....
gauss_apod = granule/scale

; > Work out how many of these fit in both x and y
division_x = FIX(sizing[1]/(gauss_apod))
division_y = FIX(sizing[2]/(gauss_apod))

half_divisionX = 0.5 * division_x
half_divisionY = 0.5 * division_y

Fe_vel_array = FLTARR(sizing[1],sizing[2],(sizing[3]-1))
Fe_shift_x = FLTARR(sizing[1],sizing[2],(sizing[3]-1))
Fe_shift_y = FLTARR(sizing[1],sizing[2],(sizing[3]-1))
Fe_angles = FLTARR(sizing[1],sizing[2],(sizing[3]-1))
Fe_array = FLTARR(4,sizing[1],sizing[2],(sizing[3]-1))

counter = 0
!p.background=255.
!p.color=0.
WINDOW,0,xsize=600,ysize=600
FOR t = 0 , (sizing[3] - 2) DO BEGIN &$
 print,'Shifts between frame '+arr2str(t,/trim)+' and '+arr2str((t+1),/trim) &$
 LOADCT,0,/SILENT &$
 tvim,FeImages[*,*,(t+1)],title='Fe Images - Time '+arr2str((t*time),/trim)+' s' &$
 FOR x = 0, (((sizing[1])/gauss_apod)-1) DO BEGIN &$
  FOR y = 0, (((sizing[2])/gauss_apod)-1) DO BEGIN &$
  	f = FeImages[(x*gauss_apod):((x*gauss_apod)+gauss_apod)-1,(y*gauss_apod):((y*gauss_apod)+gauss_apod)-1,t] &$   
  	g = FeImages[(x*gauss_apod):((x*gauss_apod)+gauss_apod)-1,(y*gauss_apod):((y*gauss_apod)+gauss_apod)-1,(t+1)] &$   

	arrinfo=SIZE(f) &$
	nx=arrinfo(1) &$
	ny=arrinfo(2) &$

	finv=fft(f,-1) &$
	ginv=fft(g,-1) &$

	spec=ginv*conj(finv) &$
	ccor=fft(spec,1) &$
	ccor=shift(ccor,nx/2,ny/2) &$
	absccor=abs(ccor) &$
	cmax=max(absccor,indmax) &$

	ix = indmax mod nx &$
	iy = indmax/nx &$

	shiftx0=float(ix) &$
	shifty0=float(iy) &$

	n_interp = 10 &$
	interp_pts = findgen(10*n_interp+1)/n_interp-5. &$
	xinterp = interp_pts+shiftx0 &$
	yinterp = interp_pts+shifty0 &$
	peakarea=interpolate(absccor,xinterp,yinterp,cubic=-0.5,/grid,missing=0.) &$
	cmaxmax=max(peakarea,indmaxmax) &$
	ixx=indmaxmax mod (n_interp*10+1) &$
	iyy=indmaxmax / (n_interp*10+1) &$
	shiftxx = xinterp(ixx) &$
	shiftyy = yinterp(iyy) &$
	shiftx=shiftxx-float(nx/2) &$
	shifty=shiftyy-float(ny/2) &$
	
	theta = ATAN(shifty / shiftx) &$ 
	theta_cube = FLTARR((gauss_apod),(gauss_apod)) &$
	theta_cube[*] = theta &$
	
	d = ((SQRT(((shiftx * scale)^2) + ((shifty * scale)^2)))) &$
	;IF (d LT (scale * 2)) THEN vel = 0. &$
	;IF (d GE (scale * 2)) THEN vel =  (d - (scale * 2)) /(time) &$
	vel = d/time &$
	;IF (vel GT 1.5) THEN vel = 0. &$
	
	cube = FLTARR((gauss_apod),(gauss_apod)) &$
	cube[*] = vel &$
	
	shiftx_cube = FLTARR((gauss_apod),(gauss_apod)) &$
	shiftx_cube[*] = shiftx &$
	shifty_cube = FLTARR((gauss_apod),(gauss_apod)) &$
	shifty_cube[*] = shifty &$

	Fe_vel_array[(x*(gauss_apod)):(((x+1)*(gauss_apod))-1),(y*(gauss_apod)):(((y+1)*(gauss_apod))-1),t] = cube &$
	Fe_angles[(x*(gauss_apod)):(((x+1)*(gauss_apod))-1),(y*(gauss_apod)):(((y+1)*(gauss_apod))-1),t] = theta_cube &$
	Fe_shift_x[(x*(gauss_apod)):(((x+1)*(gauss_apod))-1),(y*(gauss_apod)):(((y+1)*(gauss_apod))-1),t] = shiftx_cube &$
	Fe_shift_y[(x*(gauss_apod)):(((x+1)*(gauss_apod))-1),(y*(gauss_apod)):(((y+1)*(gauss_apod))-1),t] = shifty_cube &$

	print,'(x,y) shifts:',shiftx,shifty, vel
	IF (vel EQ 0.) THEN CONTINUE &$
	;IF (vel LT 0.1) THEN CONTINUE &$
	;IF (vel GT 1.5) THEN CONTINUE &$
  	LOADCT,3,/SILENT &$
	tek_color &$
ARROW,((x*gauss_apod)+half_divisionx),((y*gauss_apod)+half_divisiony),(((x*gauss_apod)+half_divisionx)+shiftx),(((y*gauss_apod)+half_divisiony)+shifty),/DATA,HSIZE=(!D.X_SIZE/192.),thick=2,color=7 &$ 
  	;Scale of 8.3 means that 1km/s will be represented by a shift of 10pixels
  ENDFOR &$
 ENDFOR &$
 max_vel = MAX(Fe_vel_array[*,*,t]) &$
 min_vel = MIN(Fe_vel_array[*,*,t]) &$
 av_vel = MEAN(Fe_vel_array[*,*,t]) &$
 tab = DBLARR(4,1) &$
 tab[0,*] = t &$
 tab[1,*] = min_vel &$
 tab[2,*] = max_vel &$
 tab[3,*] = av_vel &$
 ;OPENW,unit,''+file+'',/APPEND,WIDTH=250  &$	
 ;printf,unit,tab &$
 ;CLOSE,unit &$
 Fe_array[0,*,*,t] = Fe_vel_array[*,*,t] &$
 Fe_array[1,*,*,t] = Fe_shift_x[*,*,t] &$
 Fe_array[2,*,*,t] = Fe_shift_y[*,*,t] &$
 Fe_array[3,*,*,t] = Fe_angles[*,*,t] &$
 IF (counter le 9) THEN name = '0000' + arr2str(counter,/trim) &$
 IF (counter gt 9) AND (counter le 99) THEN name = '000' + arr2str(counter,/trim) &$
 IF (counter gt 99) AND (counter le 999) THEN name = '00' + arr2str(counter,/trim) &$
 IF (counter gt 999) AND (counter le 9999) THEN name = '0' + arr2str(counter,/trim) &$
 IF (counter gt 9999) AND (counter le 99999) THEN name = '' + arr2str(counter,/trim) &$
        
 X2JPEG,'/data/solarstore3/phk/data/27Jul2014/LCT_images/Fe_noapodisation_LCT_'+name+'.jpg' &$
 counter = counter + 1  &$
ENDFOR

Fediffimages = FLTARR(sizing[1],sizing[2],(sizing[3]-1))
FOR t = 0, sizing[3]-2 DO Fediffimages[*,*,t] = Feimages[*,*,t+1]-Feimages[*,*,t]

; >>>> TEMPORAL APODISATION NOW...

newFrame = FIX(sizing[3]/FIX(granule_time/time))
apodFeImages = FLTARR(sizing[1],sizing[2],newFrame)
FOR t = 0, newFrame-1 DO apodFeImages[*,*,t] = FeImages[*,*,t*FIX(granule_time/time)]

apodFevel_array = FLTARR(sizing[1],sizing[2],(newFrame-1))
apodFeshift_x = FLTARR(sizing[1],sizing[2],(newFrame-1))
apodFeshift_y = FLTARR(sizing[1],sizing[2],(newFrame-1))
apodFeangles = FLTARR(sizing[1],sizing[2],(newFrame-1))
apodFearray = FLTARR(4,sizing[1],sizing[2],(newFrame-1))

counter = 0
!p.background=255.
!p.color=0.
WINDOW,0,xsize=600,ysize=600
FOR t = 0 , (newframe - 2) DO BEGIN &$
 print,'Shifts between frame '+arr2str(t,/trim)+' and '+arr2str((t+1),/trim) &$
 LOADCT,0,/SILENT &$
 tvim,apodFeImages[*,*,(t+1)],title=' Fe - Temporal Apodisation - Time '+arr2str((t*time),/trim)+' s' &$
 FOR x = 0, (((sizing[1])/gauss_apod)-1) DO BEGIN &$
  FOR y = 0, (((sizing[2])/gauss_apod)-1) DO BEGIN &$
	;f = intensity_images[(x*gauss_apod):(((x+1)*gauss_apod)-1),(y*gauss_apod):(((y+1)*gauss_apod)-1),t] &$   
	;g = intensity_images[(x*gauss_apod):(((x+1)*gauss_apod)-1),(y*gauss_apod):(((y+1)*gauss_apod)-1),(t+1)] &$
  	f = apodFeImages[(x*gauss_apod):((x*gauss_apod)+gauss_apod)-1,(y*gauss_apod):((y*gauss_apod)+gauss_apod)-1,t] &$   
  	g = apodFeImages[(x*gauss_apod):((x*gauss_apod)+gauss_apod)-1,(y*gauss_apod):((y*gauss_apod)+gauss_apod)-1,(t+1)] &$   

	arrinfo=SIZE(f) &$
	nx=arrinfo(1) &$
	ny=arrinfo(2) &$

	finv=fft(f,-1) &$
	ginv=fft(g,-1) &$

	spec=ginv*conj(finv) &$
	ccor=fft(spec,1) &$
	ccor=shift(ccor,nx/2,ny/2) &$
	absccor=abs(ccor) &$
	cmax=max(absccor,indmax) &$

	ix = indmax mod nx &$
	iy = indmax/nx &$

	shiftx0=float(ix) &$
	shifty0=float(iy) &$

	n_interp = 10 &$
	interp_pts = findgen(10*n_interp+1)/n_interp-5. &$
	xinterp = interp_pts+shiftx0 &$
	yinterp = interp_pts+shifty0 &$
	peakarea=interpolate(absccor,xinterp,yinterp,cubic=-0.5,/grid,missing=0.) &$
	cmaxmax=max(peakarea,indmaxmax) &$
	ixx=indmaxmax mod (n_interp*10+1) &$
	iyy=indmaxmax / (n_interp*10+1) &$
	shiftxx = xinterp(ixx) &$
	shiftyy = yinterp(iyy) &$
	shiftx=shiftxx-float(nx/2) &$
	shifty=shiftyy-float(ny/2) &$
	
	theta = ATAN(shifty / shiftx) &$ 
	theta_cube = FLTARR((gauss_apod),(gauss_apod)) &$
	theta_cube[*] = theta &$
	
	d = ((SQRT(((shiftx * scale)^2) + ((shifty * scale)^2)))) &$
	;IF (d LT (scale * 2)) THEN vel = 0. &$
	;IF (d GE (scale * 2)) THEN vel =  (d - (scale * 2)) /(time) &$
	vel = d/time &$
	;IF (vel GT 1.5) THEN vel = 0. &$
	
	cube = FLTARR((gauss_apod),(gauss_apod)) &$
	cube[*] = vel &$
	
	shiftx_cube = FLTARR((gauss_apod),(gauss_apod)) &$
	shiftx_cube[*] = shiftx &$
	shifty_cube = FLTARR((gauss_apod),(gauss_apod)) &$
	shifty_cube[*] = shifty &$

	apodFevel_array[(x*(gauss_apod)):(((x+1)*(gauss_apod))-1),(y*(gauss_apod)):(((y+1)*(gauss_apod))-1),t] = cube &$
	apodFeangles[(x*(gauss_apod)):(((x+1)*(gauss_apod))-1),(y*(gauss_apod)):(((y+1)*(gauss_apod))-1),t] = theta_cube &$
	apodFeshift_x[(x*(gauss_apod)):(((x+1)*(gauss_apod))-1),(y*(gauss_apod)):(((y+1)*(gauss_apod))-1),t] = shiftx_cube &$
	apodFeshift_y[(x*(gauss_apod)):(((x+1)*(gauss_apod))-1),(y*(gauss_apod)):(((y+1)*(gauss_apod))-1),t] = shifty_cube &$

	print,'(x,y) shifts:',shiftx,shifty, vel
	IF (vel EQ 0.) THEN CONTINUE &$
	;IF (vel LT 0.1) THEN CONTINUE &$
	;IF (vel GT 1.5) THEN CONTINUE &$
  	LOADCT,3,/SILENT &$
	tek_color &$
ARROW,((x*gauss_apod)+half_divisionx),((y*gauss_apod)+half_divisiony),(((x*gauss_apod)+half_divisionx)+shiftx),(((y*gauss_apod)+half_divisiony)+shifty),/DATA,HSIZE=(!D.X_SIZE/192.),thick=2,color=7 &$ 
  	;Scale of 8.3 means that 1km/s will be represented by a shift of 10pixels
  ENDFOR &$
 ENDFOR &$
 max_vel = MAX(apodFevel_array[*,*,t]) &$
 min_vel = MIN(apodFevel_array[*,*,t]) &$
 av_vel = MEAN(apodFevel_array[*,*,t]) &$
 tab = DBLARR(4,1) &$
 tab[0,*] = t &$
 tab[1,*] = min_vel &$
 tab[2,*] = max_vel &$
 tab[3,*] = av_vel &$
 ;OPENW,unit,''+file+'',/APPEND,WIDTH=250  &$	
 ;printf,unit,tab &$
 ;CLOSE,unit &$
 apodFearray[0,*,*,t] = apodFevel_array[*,*,t] &$
 apodFearray[1,*,*,t] = apodFeshift_x[*,*,t] &$
 apodFearray[2,*,*,t] = apodFeshift_y[*,*,t] &$
 apodFearray[3,*,*,t] = apodFeangles[*,*,t] &$
 IF (counter le 9) THEN name = '0000' + arr2str(counter,/trim) &$
 IF (counter gt 9) AND (counter le 99) THEN name = '000' + arr2str(counter,/trim) &$
 IF (counter gt 99) AND (counter le 999) THEN name = '00' + arr2str(counter,/trim) &$
 IF (counter gt 999) AND (counter le 9999) THEN name = '0' + arr2str(counter,/trim) &$
 IF (counter gt 9999) AND (counter le 99999) THEN name = '' + arr2str(counter,/trim) &$
        
 X2JPEG,'/data/solarstore3/phk/data/27Jul2014/LCT_images/Fe_apodisation_LCT_'+name+'.jpg' &$
 counter = counter + 1  &$
ENDFOR

print,'MBP is about 1/2 way through data'

apodFediffimages = FLTARR(sizing[1],sizing[2],newFrame-1)
FOR t = 0, newFrame-2 DO apodFediffimages[*,*,t] = apodFeImages[*,*,t+1]-apodFeImages[*,*,t]

; > Save Out the Fe Data For Future Use

note1 = 'No temporal apodisation applied'
note2 = 'Temporal apodiation applied - 5min'
SAVE,FILENAME='/data/solarstore3/phk/data/27Jul2014/FeLCT_results_NoApodisation.sav',FeImages,Fediffimages,Fe_array,Fe_vel_array,Fe_shift_x,Fe_shift_y,Fe_angles,note1
SAVE,FILENAME='/data/solarstore3/phk/data/27Jul2014/FeLCT_results_Apodisation.sav',ApodFeImages,ApodFediffimages,ApodFearray,ApodFevel_array,ApodFeshift_x,ApodFeshift_y,ApodFeangles,note2


;-------------------------
; 	CA IMAGES 2nd
;-------------------------

; > Image res and the size of a granule (for spatial window)
scale = 50.
granule = 1500.

; > Time scale of scan and the half scan (for temporal windowing)
time = 33.
time_2 = time/2.

; > Lifetime of a granule (to work out the number of images needed)
granule_time = (60*5.)
half_granule_t = granule_time/2.

sizing = SIZE(CaImages)

; > Move onto the LCT of this nonsense....
gauss_apod = granule/scale

; > Work out how many of these fit in both x and y
division_x = FIX(sizing[1]/(gauss_apod))
division_y = FIX(sizing[2]/(gauss_apod))

half_divisionX = 0.5 * division_x
half_divisionY = 0.5 * division_y

Ca_vel_array = FLTARR(sizing[1],sizing[2],(sizing[3]-1))
Ca_shift_x = FLTARR(sizing[1],sizing[2],(sizing[3]-1))
Ca_shift_y = FLTARR(sizing[1],sizing[2],(sizing[3]-1))
Ca_angles = FLTARR(sizing[1],sizing[2],(sizing[3]-1))
Ca_array = FLTARR(4,sizing[1],sizing[2],(sizing[3]-1))

counter = 0
!p.background=255.
!p.color=0.
WINDOW,0,xsize=600,ysize=600
FOR t = 0 , (sizing[3] - 2) DO BEGIN &$
 print,'Shifts between frame '+arr2str(t,/trim)+' and '+arr2str((t+1),/trim) &$
 LOADCT,0,/SILENT &$
 tvim,CaImages[*,*,(t+1)],title='Ca Images - Time '+arr2str((t*time),/trim)+' s' &$
 FOR x = 0, (((sizing[1])/gauss_apod)-1) DO BEGIN &$
  FOR y = 0, (((sizing[2])/gauss_apod)-1) DO BEGIN &$
  	f = CaImages[(x*gauss_apod):((x*gauss_apod)+gauss_apod)-1,(y*gauss_apod):((y*gauss_apod)+gauss_apod)-1,t] &$   
  	g = CaImages[(x*gauss_apod):((x*gauss_apod)+gauss_apod)-1,(y*gauss_apod):((y*gauss_apod)+gauss_apod)-1,(t+1)] &$   

	arrinfo=SIZE(f) &$
	nx=arrinfo(1) &$
	ny=arrinfo(2) &$

	finv=fft(f,-1) &$
	ginv=fft(g,-1) &$

	spec=ginv*conj(finv) &$
	ccor=fft(spec,1) &$
	ccor=shift(ccor,nx/2,ny/2) &$
	absccor=abs(ccor) &$
	cmax=max(absccor,indmax) &$

	ix = indmax mod nx &$
	iy = indmax/nx &$

	shiftx0=float(ix) &$
	shifty0=float(iy) &$

	n_interp = 10 &$
	interp_pts = findgen(10*n_interp+1)/n_interp-5. &$
	xinterp = interp_pts+shiftx0 &$
	yinterp = interp_pts+shifty0 &$
	peakarea=interpolate(absccor,xinterp,yinterp,cubic=-0.5,/grid,missing=0.) &$
	cmaxmax=max(peakarea,indmaxmax) &$
	ixx=indmaxmax mod (n_interp*10+1) &$
	iyy=indmaxmax / (n_interp*10+1) &$
	shiftxx = xinterp(ixx) &$
	shiftyy = yinterp(iyy) &$
	shiftx=shiftxx-float(nx/2) &$
	shifty=shiftyy-float(ny/2) &$
	
	theta = ATAN(shifty / shiftx) &$ 
	theta_cube = FLTARR((gauss_apod),(gauss_apod)) &$
	theta_cube[*] = theta &$
	
	d = ((SQRT(((shiftx * scale)^2) + ((shifty * scale)^2)))) &$
	;IF (d LT (scale * 2)) THEN vel = 0. &$
	;IF (d GE (scale * 2)) THEN vel =  (d - (scale * 2)) /(time) &$
	vel = d/time &$
	;IF (vel GT 1.5) THEN vel = 0. &$
	
	cube = FLTARR((gauss_apod),(gauss_apod)) &$
	cube[*] = vel &$
	
	shiftx_cube = FLTARR((gauss_apod),(gauss_apod)) &$
	shiftx_cube[*] = shiftx &$
	shifty_cube = FLTARR((gauss_apod),(gauss_apod)) &$
	shifty_cube[*] = shifty &$

	Ca_vel_array[(x*(gauss_apod)):(((x+1)*(gauss_apod))-1),(y*(gauss_apod)):(((y+1)*(gauss_apod))-1),t] = cube &$
	Ca_angles[(x*(gauss_apod)):(((x+1)*(gauss_apod))-1),(y*(gauss_apod)):(((y+1)*(gauss_apod))-1),t] = theta_cube &$
	Ca_shift_x[(x*(gauss_apod)):(((x+1)*(gauss_apod))-1),(y*(gauss_apod)):(((y+1)*(gauss_apod))-1),t] = shiftx_cube &$
	Ca_shift_y[(x*(gauss_apod)):(((x+1)*(gauss_apod))-1),(y*(gauss_apod)):(((y+1)*(gauss_apod))-1),t] = shifty_cube &$

	print,'(x,y) shifts:',shiftx,shifty, vel
	IF (vel EQ 0.) THEN CONTINUE &$
	;IF (vel LT 0.1) THEN CONTINUE &$
	;IF (vel GT 1.5) THEN CONTINUE &$
  	LOADCT,3,/SILENT &$
	tek_color &$
ARROW,((x*gauss_apod)+half_divisionx),((y*gauss_apod)+half_divisiony),(((x*gauss_apod)+half_divisionx)+shiftx),(((y*gauss_apod)+half_divisiony)+shifty),/DATA,HSIZE=(!D.X_SIZE/192.),thick=2,color=7 &$ 
  	;Scale of 8.3 means that 1km/s will be represented by a shift of 10pixels
  ENDFOR &$
 ENDFOR &$
 max_vel = MAX(Ca_vel_array[*,*,t]) &$
 min_vel = MIN(Ca_vel_array[*,*,t]) &$
 av_vel = MEAN(Ca_vel_array[*,*,t]) &$
 tab = DBLARR(4,1) &$
 tab[0,*] = t &$
 tab[1,*] = min_vel &$
 tab[2,*] = max_vel &$
 tab[3,*] = av_vel &$
 ;OPENW,unit,''+file+'',/APPEND,WIDTH=250  &$	
 ;printf,unit,tab &$
 ;CLOSE,unit &$
 Ca_array[0,*,*,t] = Ca_vel_array[*,*,t] &$
 Ca_array[1,*,*,t] = Ca_shift_x[*,*,t] &$
 Ca_array[2,*,*,t] = Ca_shift_y[*,*,t] &$
 Ca_array[3,*,*,t] = Ca_angles[*,*,t] &$
 IF (counter le 9) THEN name = '0000' + arr2str(counter,/trim) &$
 IF (counter gt 9) AND (counter le 99) THEN name = '000' + arr2str(counter,/trim) &$
 IF (counter gt 99) AND (counter le 999) THEN name = '00' + arr2str(counter,/trim) &$
 IF (counter gt 999) AND (counter le 9999) THEN name = '0' + arr2str(counter,/trim) &$
 IF (counter gt 9999) AND (counter le 99999) THEN name = '' + arr2str(counter,/trim) &$
        
 X2JPEG,'/data/solarstore3/phk/data/27Jul2014/LCT_images/Ca_noapodisation_LCT_'+name+'.jpg' &$
 counter = counter + 1  &$
ENDFOR

Cadiffimages = FLTARR(sizing[1],sizing[2],(sizing[3]-1))
FOR t = 0, sizing[3]-2 DO Cadiffimages[*,*,t] = Caimages[*,*,t+1]-Caimages[*,*,t]

; >>>> TEMPORAL APODISATION NOW...

newFrame = FIX(sizing[3]/FIX(granule_time/time))
apodCaImages = FLTARR(sizing[1],sizing[2],newFrame)
FOR t = 0, newFrame-1 DO apodCaImages[*,*,t] = CaImages[*,*,t*FIX(granule_time/time)]

apodCavel_array = FLTARR(sizing[1],sizing[2],(newFrame-1))
apodCashift_x = FLTARR(sizing[1],sizing[2],(newFrame-1))
apodCashift_y = FLTARR(sizing[1],sizing[2],(newFrame-1))
apodCaangles = FLTARR(sizing[1],sizing[2],(newFrame-1))
apodCaarray = FLTARR(4,sizing[1],sizing[2],(newFrame-1))

counter = 0
!p.background=255.
!p.color=0.
WINDOW,0,xsize=600,ysize=600
FOR t = 0 , (newframe - 2) DO BEGIN &$
 print,'Shifts between frame '+arr2str(t,/trim)+' and '+arr2str((t+1),/trim) &$
 LOADCT,0,/SILENT &$
 tvim,apodCaImages[*,*,(t+1)],title=' Ca - Temporal Apodisation - Time '+arr2str((t*time),/trim)+' s' &$
 FOR x = 0, (((sizing[1])/gauss_apod)-1) DO BEGIN &$
  FOR y = 0, (((sizing[2])/gauss_apod)-1) DO BEGIN &$
	;f = intensity_images[(x*gauss_apod):(((x+1)*gauss_apod)-1),(y*gauss_apod):(((y+1)*gauss_apod)-1),t] &$   
	;g = intensity_images[(x*gauss_apod):(((x+1)*gauss_apod)-1),(y*gauss_apod):(((y+1)*gauss_apod)-1),(t+1)] &$
  	f = apodCaImages[(x*gauss_apod):((x*gauss_apod)+gauss_apod)-1,(y*gauss_apod):((y*gauss_apod)+gauss_apod)-1,t] &$   
  	g = apodCaImages[(x*gauss_apod):((x*gauss_apod)+gauss_apod)-1,(y*gauss_apod):((y*gauss_apod)+gauss_apod)-1,(t+1)] &$   

	arrinfo=SIZE(f) &$
	nx=arrinfo(1) &$
	ny=arrinfo(2) &$

	finv=fft(f,-1) &$
	ginv=fft(g,-1) &$

	spec=ginv*conj(finv) &$
	ccor=fft(spec,1) &$
	ccor=shift(ccor,nx/2,ny/2) &$
	absccor=abs(ccor) &$
	cmax=max(absccor,indmax) &$

	ix = indmax mod nx &$
	iy = indmax/nx &$

	shiftx0=float(ix) &$
	shifty0=float(iy) &$

	n_interp = 10 &$
	interp_pts = findgen(10*n_interp+1)/n_interp-5. &$
	xinterp = interp_pts+shiftx0 &$
	yinterp = interp_pts+shifty0 &$
	peakarea=interpolate(absccor,xinterp,yinterp,cubic=-0.5,/grid,missing=0.) &$
	cmaxmax=max(peakarea,indmaxmax) &$
	ixx=indmaxmax mod (n_interp*10+1) &$
	iyy=indmaxmax / (n_interp*10+1) &$
	shiftxx = xinterp(ixx) &$
	shiftyy = yinterp(iyy) &$
	shiftx=shiftxx-float(nx/2) &$
	shifty=shiftyy-float(ny/2) &$
	
	theta = ATAN(shifty / shiftx) &$ 
	theta_cube = FLTARR((gauss_apod),(gauss_apod)) &$
	theta_cube[*] = theta &$
	
	d = ((SQRT(((shiftx * scale)^2) + ((shifty * scale)^2)))) &$
	;IF (d LT (scale * 2)) THEN vel = 0. &$
	;IF (d GE (scale * 2)) THEN vel =  (d - (scale * 2)) /(time) &$
	vel = d/time &$
	;IF (vel GT 1.5) THEN vel = 0. &$
	
	cube = FLTARR((gauss_apod),(gauss_apod)) &$
	cube[*] = vel &$
	
	shiftx_cube = FLTARR((gauss_apod),(gauss_apod)) &$
	shiftx_cube[*] = shiftx &$
	shifty_cube = FLTARR((gauss_apod),(gauss_apod)) &$
	shifty_cube[*] = shifty &$

	apodCavel_array[(x*(gauss_apod)):(((x+1)*(gauss_apod))-1),(y*(gauss_apod)):(((y+1)*(gauss_apod))-1),t] = cube &$
	apodCaangles[(x*(gauss_apod)):(((x+1)*(gauss_apod))-1),(y*(gauss_apod)):(((y+1)*(gauss_apod))-1),t] = theta_cube &$
	apodCashift_x[(x*(gauss_apod)):(((x+1)*(gauss_apod))-1),(y*(gauss_apod)):(((y+1)*(gauss_apod))-1),t] = shiftx_cube &$
	apodCashift_y[(x*(gauss_apod)):(((x+1)*(gauss_apod))-1),(y*(gauss_apod)):(((y+1)*(gauss_apod))-1),t] = shifty_cube &$

	print,'(x,y) shifts:',shiftx,shifty, vel
	IF (vel EQ 0.) THEN CONTINUE &$
	;IF (vel LT 0.1) THEN CONTINUE &$
	;IF (vel GT 1.5) THEN CONTINUE &$
  	LOADCT,3,/SILENT &$
	tek_color &$
ARROW,((x*gauss_apod)+half_divisionx),((y*gauss_apod)+half_divisiony),(((x*gauss_apod)+half_divisionx)+shiftx),(((y*gauss_apod)+half_divisiony)+shifty),/DATA,HSIZE=(!D.X_SIZE/192.),thick=2,color=7 &$ 
  	;Scale of 8.3 means that 1km/s will be represented by a shift of 10pixels
  ENDFOR &$
 ENDFOR &$
 max_vel = MAX(apodCavel_array[*,*,t]) &$
 min_vel = MIN(apodCavel_array[*,*,t]) &$
 av_vel = MEAN(apodCavel_array[*,*,t]) &$
 tab = DBLARR(4,1) &$
 tab[0,*] = t &$
 tab[1,*] = min_vel &$
 tab[2,*] = max_vel &$
 tab[3,*] = av_vel &$
 ;OPENW,unit,''+file+'',/APPEND,WIDTH=250  &$	
 ;printf,unit,tab &$
 ;CLOSE,unit &$
 apodCaarray[0,*,*,t] = apodCavel_array[*,*,t] &$
 apodCaarray[1,*,*,t] = apodCashift_x[*,*,t] &$
 apodCaarray[2,*,*,t] = apodCashift_y[*,*,t] &$
 apodCaarray[3,*,*,t] = apodCaangles[*,*,t] &$
 IF (counter le 9) THEN name = '0000' + arr2str(counter,/trim) &$
 IF (counter gt 9) AND (counter le 99) THEN name = '000' + arr2str(counter,/trim) &$
 IF (counter gt 99) AND (counter le 999) THEN name = '00' + arr2str(counter,/trim) &$
 IF (counter gt 999) AND (counter le 9999) THEN name = '0' + arr2str(counter,/trim) &$
 IF (counter gt 9999) AND (counter le 99999) THEN name = '' + arr2str(counter,/trim) &$
        
 X2JPEG,'/data/solarstore3/phk/data/27Jul2014/LCT_images/Ca_apodisation_LCT_'+name+'.jpg' &$
 counter = counter + 1  &$
ENDFOR

apodCadiffimages = FLTARR(sizing[1],sizing[2],newFrame-1)
FOR t = 0, newFrame-2 DO apodCadiffimages[*,*,t] = apodCaImages[*,*,t+1]-apodCaImages[*,*,t]
LOADCT,0,/SILENT
; > Save Out the Ca Data For Future Use

note1 = 'No temporal apodisation applied'
note2 = 'Temporal apodiation applied - 5min'
SAVE,FILENAME='/data/solarstore3/phk/data/27Jul2014/CaLCT_results_NoApodisation.sav',CaImages,Cadiffimages,Ca_array,Ca_vel_array,Ca_shift_x,Ca_shift_y,Ca_angles,note1
SAVE,FILENAME='/data/solarstore3/phk/data/27Jul2014/CaLCT_results_Apodisation.sav',ApodCaImages,ApodCadiffimages,ApodCaarray,ApodCavel_array,ApodCashift_x,ApodCashift_y,ApodCaangles,note2


END
