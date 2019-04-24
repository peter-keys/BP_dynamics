;-----------------------------------------------------------------------------
;	LCT TO WORK OUT THE MOTION OF GRANULES ABOUT MBPS TO SEE IF FORCED
;-----------------------------------------------------------------------------
; 	- Want to work out the motion of granules about MBPs
;	- Use LCT
; 	- Will give an indication of direction of granular motion
;	- Pair this with the MBP evolutions to see if forcing occurs on MBPs

;-----------------------------------------------------------------------------
; 		INDIVIDUAL MBP APPROACH (1) [LCT about MBP]
;-----------------------------------------------------------------------------

;-------------------------------------------------------
; 	LOAD IN MBP INFO...
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

readcol,'/data/solarstore3/phk/data/27Jul2014/final_tracking_file',f='D,D,D,D,D,x', btpt_num, frame, xloc, yloc, area

MBP_num = 35.
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

openr, lun, file, /get_lun 
dat = assoc(lun, fltarr(nx, ny, nw, 4, /nozero), 512) 
d = float(dat[initial_frame]) &$ ;; for example snapshot 10 		; USE 30 here as it is the 1st frame of this MBP (174) 
free_lun, lun 
im_dims = SIZE(d)

	; > Ensure the frames don't appear before/after the data
IF (initial_frame LT 0) THEN initial_frame = 0. 
IF (final_frame GT 91) THEN final_frame = 91. 	

	; > Calculate the biggest/lowest x and y values to get dimensiosn of sub field to invert
	; > Use IFs to make sure the low and high values aren't outside image dimensions
;..........
; See below, but htis might not be the best approach if you 
; are trying to create an array as a square.
; Alternatively, take the MEAN location of the MBP and 
; get values for the box from that. Will have issues on dimensions 
; if the box falls outside the dimensions of the images
; I'll add in a print statement to make me aware of that scenario
; Not sure how I would automate a fix for that....
;x0 = MIN(xcoords)-200 
;x1 = MAX(xcoords)+200 
;y0 = MIN(ycoords)-200 
;y1 = MAX(ycoords)+200 	
;..........

x0 = FIX(MEAN(xcoords)-200)
x1 = FIX(MEAN(xcoords)+200)
y0 = FIX(MEAN(ycoords)-200)
y1 = FIX(MEAN(ycoords)+200)

IF (x0 LT 0) THEN x0 = 0. 
IF (x1 GT im_dims[1]) THEN x1 = im_dims[1] - 1. 
IF (y0 LT 0) THEN y0 = 0. 
IF (y1 GT im_dims[2]) THEN y1 = im_dims[2] - 1. 	

xscale = (x1 - x0) + 1
yscale = (y1 - y0) + 1

IF (xscale NE yscale) THEN print,'x an y values dont match... you wont have a square box mate...'

;..............
;Trying doing this to force them to be equal...
IF (xscale NE yscale) THEN BEGIN &$
	print,' ' &$
	print,'x an y values dont match... you wont have a square box mate...' &$
	; - Get the smaller number
	dims = MIN([xscale,yscale]) &$
	; - Get the 1/2 value
	half_dim = FIX(dims/2.) &$
	; - Add in a bit for working out if it is odd or even & what to do
	IF (half_dim*2 NE dims) THEN BEGIN &$
		dimposlow = half_dim &$
		dimposhigh = (dims - half_dim)-1 &$
	ENDIF ELSE BEGIN &$
		dimposlow = half_dim &$
		dimposhigh = half_dim &$
	ENDELSE &$
	; - dimposlow/high are the low/high values for the box being resized
	IF (xscale GT yscale) THEN BEGIN &$ 
		x0 = FIX(MEAN(xcoords))-dimposlow &$ 
		x1 = FIX(MEAN(xcoords))+dimposhigh &$
	ENDIF &$
	IF (yscale GT xscale) THEN BEGIN &$
		y0 = FIX(MEAN(ycoords))-dimposlow &$ 
		y1 = FIX(MEAN(ycoords))+dimposhigh &$
	ENDIF &$
	print, 'Dimensions now being changed to: ',x1-x0+1,y1-y0+1 &$
	print,' ' &$
ENDIF
;..............

nx1 = x1 - x0 + 1
ny1 = y1 - y0 + 1
nf1 = final_frame - initial_frame + 1

fs = INDGEN(nf1)+initial_frame

;-------------------------------------------------------
; 	LOAD IN IMAGES FOR THE LCT TO WORK
;-------------------------------------------------------
; - Want to load in some images to do LCT
; - Would like it to be centred on the MBP in time/space to make things easier
; - To do that need to do some LCT related values to get the windowding right

; > Image res and the size of a granule (for spatial window)
scale = 50.
granule = 1500.

; > Time scale of scan and the half scan (for temporal windowing)
time = 33.
time_2 = time/2.

; > Lifetime of a granule (to work out the number of images needed)
granule_time = (60*5.)
half_granule_t = granule_time/2.

; > Want 3 granule lifecycles either side of the MBP (this value ensures that)
start_end_frames = (half_granule_t * 3.)/time

start_frame = initial_frame - ROUND(start_end_frames)
End_frame = final_frame + ROUND(start_end_frames)
number_frames = (end_frame - start_frame) + 1

; > Add in 2 IF statements to get rid of any issues whereby the frame is outside the data...
IF (start_frame LT 0) THEN start_frame = 0.
IF (End_frame GT nt) THEN End_frame = (nt-1)

; > Get the dimensions of the images you want. (You may want to force xy to be a square)

;..........
; note doing it this way of taking the largest x/y won't work cause you'd have to 
; try and work out some way of 'fixing' the smaller dimensions to that of the other.
; Also have to be careful with those near the origin/end of image.
;I'll keep this here for prosertity, but I need to change the x0,y0 above...
;
;Think I fixed it for the case where I fix the values to the smaller.
;xscale = MAX([(x1-x0+1),(y1-y0+1)])
;yscale = xscale
;..........

fscale = (End_frame-start_frame) +1
intensity_images = FLTARR(nx1,ny1,fscale)

IF (FIX(xscale/(granule/scale)) LT 10) THEN print,'Youre below 10 granules in image dimensions... Will be dodgy LCT....'

; > Start a loop to get the frames for the LCT to be run on
FOR t= start_frame, End_frame DO BEGIN &$
	openr, lun, file, /get_lun  &$
		dat = assoc(lun, fltarr(nx, ny, nw, 4, /nozero), 512)  &$
		d = float(dat[t]) &$ ;; for example snapshot 10 		; USE 30 here as it is the 1st frame of this MBP (174) 
	free_lun, lun  &$
	intensity_images[*,*,t] = d[x0:x1,y0:y1,16,0] &$
ENDFOR

; > Need to check the actual images as there may be some drift in the image edges

tvim,intensity_images[*,*,0],title='Sample image from 1st image in the data'
wait,5
tvim,intensity_images[*,*,FIX(fscale*(1/4.))],title='Sample image from 1/4 in to the data'
wait,5
tvim,intensity_images[*,*,FIX(fscale*(2/4.))],title='Sample image from 1/2 in to the data'
wait,5
tvim,intensity_images[*,*,FIX(fscale*(3/4.))],title='Sample image from 3/4 in to the data'
wait,5
tvim,intensity_images[*,*,fscale-1],title='Sample image from last image in the data'
wait,5

;..........
; > This should be commented out when you are using this as a standalone program
;xstepper,intensity_images,xsize=600,ysize=600
;
; > Work out the part where you should trim the data
;tvim,intensity_images[*,*,3]
;.r picks
;pick
;..........

; > Resize to remove crap from the boundaries
intensity_images = intensity_images[0:299,*,*]

; > Get new sizes of the resized images and check for the 

resizedims = SIZE(intensity_images)
rexsize = resizedims[1]
reysize = resizedims[2]

; > Check that the new images are still useful for the granule scaling etc...
print,rexsize*scale/granule,reysize*scale/granule

; > Make a square again
IF (rexsize GT reysize) THEN intensity_images = intensity_images[0:reysize-1,*,*]
IF (reysize GT rexsize) THEN intensity_images = intensity_images[*,0:rexsize-1,*]
sizing = SIZE(intensity_images)

; > Move onto the LCT of this nonsense....
gauss_apod = granule/scale

; > Work out how many of these fit in both x and y
division_x = FIX(sizing[1]/(gauss_apod))
division_y = FIX(sizing[2]/(gauss_apod))

half_division = 0.5 * division_x
vel_array = FLTARR(sizing[1],sizing[2],(fscale-1))
shift_x = FLTARR(sizing[1],sizing[2],(fscale-1))
shift_y = FLTARR(sizing[1],sizing[2],(fscale-1))
angles = FLTARR(sizing[1],sizing[2],(fscale-1))
array = FLTARR(4,sizing[1],sizing[2],(fscale-1))

counter = 0
!p.background=255.
!p.color=0.
WINDOW,0,xsize=600,ysize=600
FOR t = 0 , (fscale - 2) DO BEGIN &$
 print,'Shifts between frame '+arr2str(t,/trim)+' and '+arr2str((t+1),/trim) &$
 LOADCT,0,/SILENT &$
 tvim,intensity_images[*,*,(t+1)],title='Time '+arr2str((t*time),/trim)+' s' &$
 FOR x = 0, (((sizing[1])/gauss_apod)-1) DO BEGIN &$
  FOR y = 0, (((sizing[2])/gauss_apod)-1) DO BEGIN &$
	;f = intensity_images[(x*gauss_apod):(((x+1)*gauss_apod)-1),(y*gauss_apod):(((y+1)*gauss_apod)-1),t] &$   
	;g = intensity_images[(x*gauss_apod):(((x+1)*gauss_apod)-1),(y*gauss_apod):(((y+1)*gauss_apod)-1),(t+1)] &$
  	f = intensity_images[(x*gauss_apod):((x*gauss_apod)+gauss_apod)-1,(y*gauss_apod):((y*gauss_apod)+gauss_apod)-1,t] &$   
  	g = intensity_images[(x*gauss_apod):((x*gauss_apod)+gauss_apod)-1,(y*gauss_apod):((y*gauss_apod)+gauss_apod)-1,(t+1)] &$   

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

	vel_array[(x*(gauss_apod)):(((x+1)*(gauss_apod))-1),(y*(gauss_apod)):(((y+1)*(gauss_apod))-1),t] = cube &$
	angles[(x*(gauss_apod)):(((x+1)*(gauss_apod))-1),(y*(gauss_apod)):(((y+1)*(gauss_apod))-1),t] = theta_cube &$
	shift_x[(x*(gauss_apod)):(((x+1)*(gauss_apod))-1),(y*(gauss_apod)):(((y+1)*(gauss_apod))-1),t] = shiftx_cube &$
	shift_y[(x*(gauss_apod)):(((x+1)*(gauss_apod))-1),(y*(gauss_apod)):(((y+1)*(gauss_apod))-1),t] = shifty_cube &$

	print,'(x,y) shifts:',shiftx,shifty, vel
	IF (vel EQ 0.) THEN CONTINUE &$
	;IF (vel LT 0.1) THEN CONTINUE &$
	;IF (vel GT 1.5) THEN CONTINUE &$
  	LOADCT,3,/SILENT &$
	tek_color &$
ARROW,((x*gauss_apod)+half_division),((y*gauss_apod)+half_division),(((x*gauss_apod)+half_division)+shiftx),(((y*gauss_apod)+half_division)+shifty),/DATA,HSIZE=(!D.X_SIZE/192.),thick=2,color=7 &$ 
  	;Scale of 8.3 means that 1km/s will be represented by a shift of 10pixels
  ENDFOR &$
 ENDFOR &$
 max_vel = MAX(vel_array[*,*,t]) &$
 min_vel = MIN(vel_array[*,*,t]) &$
 av_vel = MEAN(vel_array[*,*,t]) &$
 tab = DBLARR(4,1) &$
 tab[0,*] = t &$
 tab[1,*] = min_vel &$
 tab[2,*] = max_vel &$
 tab[3,*] = av_vel &$
 ;OPENW,unit,''+file+'',/APPEND,WIDTH=250  &$	
 ;printf,unit,tab &$
 ;CLOSE,unit &$
 array[0,*,*,t] = vel_array[*,*,t] &$
 array[1,*,*,t] = shift_x[*,*,t] &$
 array[2,*,*,t] = shift_y[*,*,t] &$
 array[3,*,*,t] = angles[*,*,t] &$
 IF (counter le 9) THEN name = '0000' + arr2str(counter,/trim) &$
 IF (counter gt 9) AND (counter le 99) THEN name = '000' + arr2str(counter,/trim) &$
 IF (counter gt 99) AND (counter le 999) THEN name = '00' + arr2str(counter,/trim) &$
 IF (counter gt 999) AND (counter le 9999) THEN name = '0' + arr2str(counter,/trim) &$
 IF (counter gt 9999) AND (counter le 99999) THEN name = '' + arr2str(counter,/trim) &$
        
 ;X2JPEG,'Global_LCT_'+name+'.jpg' &$
 counter = counter + 1  &$
ENDFOR

; > Produce some difference images too of the subfield. All could be useful...

diffimages = FLTARR(sizing[1],sizing[2],(fscale-1))
FOR t = 0, fscale-2 DO diffimages[*,*,t] = intensity_images[*,*,t+1]-intensity_images[*,*,t]

; XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
; YOU NEED TO CHANGE THIS PATH 
; XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
; > Save them out So I can overwrite with more useful stuff....
print,'Saving out this info'
SAVE,FILENAME='/data/solarstore3/phk/data/27Jul2014/StrongMBPs/interesting_plots/MBP035/LCT_results_MBP35_notime_apod.sav',intensity_images,diffimages,array,vel_array,shift_x,shift_y,angles

; > Remember white (+VE) values are those representing where the feature has moved to while 
; > Darker hues (-VE) values are where the thing moved from. Null (Zero) values are for things 
; > that remain where they are.
; 
; > You have to be careful though, because in bad frames it won't work. You end up getting 
; > something that looks like a granule or an inverse granule depending on whether the poor 
; > frame was before/after the good one....

;-------------------------------------------------------
;  RUN LCT WITH THE PROPER TEMPORAL APODISATION NOW TOO
;-------------------------------------------------------

print,' '
print,'Given the dimensions in time, the temporal apodisation will be...',FIX(granule_time/time),' frames'
print,'This cuts your image length from ',fscale,' to ', FIX(fscale/FIX(granule_time/time))
print,' '

; > Lifetime of a granule (to work out the number of images needed)
granule_time = (60*5.)
half_granule_t = granule_time/2.

newFrame = FIX(fscale/FIX(granule_time/time))
nintensity_images=intensity_images[*,*,0:newFrame-1]
FOR t = 0, newFrame-1 DO nintensity_images[*,*,t] = intensity_images[*,*,t*FIX(granule_time/time)]

nvel_array = FLTARR(sizing[1],sizing[2],(newFrame-1))
nshift_x = FLTARR(sizing[1],sizing[2],(newFrame-1))
nshift_y = FLTARR(sizing[1],sizing[2],(newFrame-1))
nangles = FLTARR(sizing[1],sizing[2],(newFrame-1))
narray = FLTARR(4,sizing[1],sizing[2],(newFrame-1))

counter = 0
!p.background=255.
!p.color=0.
WINDOW,0,xsize=600,ysize=600
FOR t = 0 , (newframe - 2) DO BEGIN &$
 print,'Shifts between frame '+arr2str(t,/trim)+' and '+arr2str((t+1),/trim) &$
 LOADCT,0,/SILENT &$
 tvim,nintensity_images[*,*,(t+1)],title='Time '+arr2str((t*time),/trim)+' s' &$
 FOR x = 0, (((sizing[1])/gauss_apod)-1) DO BEGIN &$
  FOR y = 0, (((sizing[2])/gauss_apod)-1) DO BEGIN &$
	;f = intensity_images[(x*gauss_apod):(((x+1)*gauss_apod)-1),(y*gauss_apod):(((y+1)*gauss_apod)-1),t] &$   
	;g = intensity_images[(x*gauss_apod):(((x+1)*gauss_apod)-1),(y*gauss_apod):(((y+1)*gauss_apod)-1),(t+1)] &$
  	f = nintensity_images[(x*gauss_apod):((x*gauss_apod)+gauss_apod)-1,(y*gauss_apod):((y*gauss_apod)+gauss_apod)-1,t] &$   
  	g = nintensity_images[(x*gauss_apod):((x*gauss_apod)+gauss_apod)-1,(y*gauss_apod):((y*gauss_apod)+gauss_apod)-1,(t+1)] &$   

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

	nvel_array[(x*(gauss_apod)):(((x+1)*(gauss_apod))-1),(y*(gauss_apod)):(((y+1)*(gauss_apod))-1),t] = cube &$
	nangles[(x*(gauss_apod)):(((x+1)*(gauss_apod))-1),(y*(gauss_apod)):(((y+1)*(gauss_apod))-1),t] = theta_cube &$
	nshift_x[(x*(gauss_apod)):(((x+1)*(gauss_apod))-1),(y*(gauss_apod)):(((y+1)*(gauss_apod))-1),t] = shiftx_cube &$
	nshift_y[(x*(gauss_apod)):(((x+1)*(gauss_apod))-1),(y*(gauss_apod)):(((y+1)*(gauss_apod))-1),t] = shifty_cube &$

	print,'(x,y) shifts:',shiftx,shifty, vel
	IF (vel EQ 0.) THEN CONTINUE &$
	;IF (vel LT 0.1) THEN CONTINUE &$
	;IF (vel GT 1.5) THEN CONTINUE &$
  	LOADCT,3,/SILENT &$
	tek_color &$
ARROW,((x*gauss_apod)+half_division),((y*gauss_apod)+half_division),(((x*gauss_apod)+half_division)+shiftx),(((y*gauss_apod)+half_division)+shifty),/DATA,HSIZE=(!D.X_SIZE/192.),thick=2,color=7 &$ 
  	;Scale of 8.3 means that 1km/s will be represented by a shift of 10pixels
  ENDFOR &$
 ENDFOR &$
 max_vel = MAX(vel_array[*,*,t]) &$
 min_vel = MIN(vel_array[*,*,t]) &$
 av_vel = MEAN(vel_array[*,*,t]) &$
 tab = DBLARR(4,1) &$
 tab[0,*] = t &$
 tab[1,*] = min_vel &$
 tab[2,*] = max_vel &$
 tab[3,*] = av_vel &$
 ;OPENW,unit,''+file+'',/APPEND,WIDTH=250  &$	
 ;printf,unit,tab &$
 ;CLOSE,unit &$
 narray[0,*,*,t] = vel_array[*,*,t] &$
 narray[1,*,*,t] = shift_x[*,*,t] &$
 narray[2,*,*,t] = shift_y[*,*,t] &$
 narray[3,*,*,t] = angles[*,*,t] &$
 IF (counter le 9) THEN name = '0000' + arr2str(counter,/trim) &$
 IF (counter gt 9) AND (counter le 99) THEN name = '000' + arr2str(counter,/trim) &$
 IF (counter gt 99) AND (counter le 999) THEN name = '00' + arr2str(counter,/trim) &$
 IF (counter gt 999) AND (counter le 9999) THEN name = '0' + arr2str(counter,/trim) &$
 IF (counter gt 9999) AND (counter le 99999) THEN name = '' + arr2str(counter,/trim) &$
        
 ;X2JPEG,'Global_LCT_'+name+'.jpg' &$
 counter = counter + 1  &$
ENDFOR

print,'MBP is about 1/2 way through data'

ndiffimages = FLTARR(sizing[1],sizing[2],newFrame-1)
FOR t = 0, newFrame-2 DO ndiffimages[*,*,t] = nintensity_images[*,*,t+1]-nintensity_images[*,*,t]


print,'Saving out this info'
SAVE,FILENAME='/data/solarstore3/phk/data/27Jul2014/StrongMBPs/interesting_plots/MBP035/LCT_results_MBP35_time_apod.sav',nintensity_images,ndiffimages,narray,nvel_array,nshift_x,nshift_y,nangles
LOADCT,0,/SILENT

END
