;+
;ROUTINE:	LCT_data
;
;PURPOSE:	LCT squares of a data cube from frame to frame to establish a velocity field
;
;USEAGE:	LCT_data,datacube,distance,time
;
;INPUT:		The datacube under consideration 
;		The name of a file to output values to 
;		The distance of each pixel in the data set
;		The time for each frame in the data set
;		The number of cubes that the data cube is to be divided into in x and y for LCT
;
;OUTPUT:	Velocities of the data cube at each of the squares that the data is divided into 
;		The Max, Min and Average velocities in each frame
;		The inclination angle to be used in the future for plotting arrows to show flows
;
;AUTHOR:	Peter H. Keys, January 2012
;
;-

FUNCTION LCT_data,datacube,file=file,distance,time,division

sizing = SIZE(datacube,/dimensions)
IF N_ELEMENTS(sizing) eq 2. THEN frames = 1.
IF N_ELEMENTS(sizing) eq 3. THEN frames = sizing[2]

IF NOT KEYWORD_SET(distance) THEN distance = 50.
IF NOT KEYWORD_SET(time) THEN time = 2.112
IF NOT KEYWORD_SET(division) THEN time = 60.

half_division = 0.5 * division

GET_LUN, unit

OPENW,unit, ''+file+''
printf,unit,'       Frame          Min Velocity   Max Velocity     Average Velocity'
CLOSE,unit


vel_array = FLTARR(sizing[0],sizing[1],(frames-1))
angles = FLTARR(sizing[0],sizing[1],(frames-1))
array = FLTARR(2,sizing[0],sizing[1],(frames-1))

counter = 0
!p.background=255.
!p.color=0.
WINDOW,0,xsize=600,ysize=600
FOR t = 0 , (frames - 2) DO BEGIN &$
 print,'Shifts between frame '+arr2str(t,/trim)+' and '+arr2str((t+1),/trim) &$
 LOADCT,0,/SILENT &$
 tvframe,datacube[*,*,(t+1)]>0.6<1.8,ticklen=-0.015;,title='Time '+arr2str((t*time),/trim)+' s' &$
 FOR x = 0, (((sizing[0])/division)-1) DO BEGIN &$
  FOR y = 0, (((sizing[1])/division)-1) DO BEGIN &$
	f = datacube[(x*division):(((x+1)*division)-1),(y*division):(((y+1)*division)-1),t] &$   
	g = datacube[(x*division):(((x+1)*division)-1),(y*division):(((y+1)*division)-1),(t+1)] &$
  
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
	theta_cube = FLTARR((division),(division)) &$
	theta_cube[*] = theta &$
	
	d = ((SQRT(((shiftx * distance)^2) + ((shifty * distance)^2)))) &$
	;IF (d LT (distance * 2)) THEN vel = 0. &$
	;IF (d GE (distance * 2)) THEN vel =  (d - (distance * 2)) /(time) &$
	vel = d/time &$
	IF (vel GT 1.5) THEN vel = 0. &$
	
	cube = FLTARR((division),(division)) &$
	cube[*] = vel &$
	
	vel_array[(x*(division)):(((x+1)*(division))-1),(y*(division)):(((y+1)*(division))-1),t] = cube &$
	angles[(x*(division)):(((x+1)*(division))-1),(y*(division)):(((y+1)*(division))-1),t] = theta_cube &$
	print,'(x,y) shifts:',shiftx,shifty, vel
	IF (vel LT 0.1) THEN CONTINUE &$
	IF (vel GT 1.5) THEN CONTINUE &$
  	LOADCT,3,/SILENT &$
	tek_color &$
ARROW,((x*division)+half_division),((y*division)+half_division),(((x*division)+half_division)+shiftx*20.8333),(((y*division)+half_division)+shifty*20.8333),/DATA,HSIZE=(!D.X_SIZE/192.),thick=2,color=7 &$ 
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
 OPENW,unit,''+file+'',/APPEND,WIDTH=250  &$	
 printf,unit,tab &$
 CLOSE,unit &$
 array[0,*,*,t] = vel_array[*,*,t] &$
 array[1,*,*,t] = angles[*,*,t] &$
 IF (counter le 9) THEN name = '0000' + arr2str(counter,/trim) &$
 IF (counter gt 9) AND (counter le 99) THEN name = '000' + arr2str(counter,/trim) &$
 IF (counter gt 99) AND (counter le 999) THEN name = '00' + arr2str(counter,/trim) &$
 IF (counter gt 999) AND (counter le 9999) THEN name = '0' + arr2str(counter,/trim) &$
 IF (counter gt 9999) AND (counter le 99999) THEN name = '' + arr2str(counter,/trim) &$
        
 X2JPEG,'Global_LCT_'+name+'.jpg' &$
 counter = counter + 1  &$
ENDFOR
LOADCT,0,/SILENT
FREE_LUN,unit
readcol,''+file+'',f='x,x,D,D',maxes,avs,SKIPLINE=1
overall_max = MAX(maxes)
overall_av = MEAN(avs)
print,'  '
print,'Overall Max. Velocity = '+arr2str(overall_max,/trim)
print,'Overall Average Velocity = '+arr2str(overall_av,/trim)
print,'  '
RETURN,array
END
