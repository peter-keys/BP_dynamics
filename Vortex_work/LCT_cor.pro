;+
;ROUTINE:	LCT_cor
;
;PURPOSE:	Uses 2D arrays to find the shifts in features which correlate and plots the shifts with an arrow
;
;USEAGE:	cortrack,f,g,shiftx,shifty
;
;INPUT:		f,g => 2D arrays used for correlation tracking
;		shiftx,shifty => the shifts in x and y if known
;
;OUTPUT:	Shift of feature
;		Velocity of features
;		Plots with arrows showing shifts
;
;AUTHOR:	
;
;-

pro cortrack,f,g,frame1,frame2,shiftx,shifty

;============================================================================================
;		USES CORRELATION TRACKING TO CALC. THE SHIFT IN X AND Y
;============================================================================================

arrinfo=size(f)
nx=arrinfo(1)
ny=arrinfo(2)

finv=fft(f,-1)
ginv=fft(g,-1)

spec=ginv*conj(finv)
ccor=fft(spec,1)
ccor=shift(ccor,nx/2,ny/2)
absccor=abs(ccor)
cmax=max(absccor,indmax)

ix = indmax mod nx
iy = indmax/nx

shiftx0=float(ix)
shifty0=float(iy)

n_interp = 10
interp_pts = findgen(10*n_interp+1)/n_interp-5.
xinterp = interp_pts+shiftx0
yinterp = interp_pts+shifty0
peakarea=interpolate(absccor,xinterp,yinterp,cubic=-0.5,/grid,missing=0.)
cmaxmax=max(peakarea,indmaxmax)
ixx=indmaxmax mod (n_interp*10+1)
iyy=indmaxmax / (n_interp*10+1)
shiftxx = xinterp(ixx)
shiftyy = yinterp(iyy)
shiftx=shiftxx-float(nx/2)
shifty=shiftyy-float(ny/2)

print,'lct done...'
print,'(x,y) shifts:',shiftx,shifty

;============================================================================================
;		VELOCITY IS CALCULATED USING VARIABLES FROM ROSA ETC...
;============================================================================================

distance = 50.
;time = 0.528
time=2.112
;distance = 25.
;time = 16.
;time = 8.7
;time = 4.224

;distance=60.
;time=39.703

d = ((SQRT(((shiftx * distance)^2) + ((shifty * distance)^2)))) 
IF (d LT (distance * 2)) THEN vel = 0.
IF (d GE (distance * 2)) THEN vel =  (d - (distance * 2)) /(time*((frame2-frame1)+1))

print,'Velocity = ',vel

;stop

;============================================================================================
;		PLOTS IMAGES OF THE 2 ARRAYS WITH THE SHIFTS MARKED BY AN ARROW
;============================================================================================

mult,2,1

tvim,f
tvim,g
arrow,ix,iy,ix+shiftx,iy+shifty,/data,color=255

mult,1,1

end
