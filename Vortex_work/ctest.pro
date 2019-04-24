pro ctrack,f,g,shiftx,shifty

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

distance = 75.
time = 0.528
vel = ((SQRT(((shiftx * distance)^2) + ((shifty * distance)^2))) - (distance * 2))/ time 
IF (vel lt 0) THEN vel = 0
print,'Velocity = ',vel

;stop

end


npt=100

seed=15.0

rlev=1.0d0


a1=dblarr(npt,npt)
a2=dblarr(npt,npt)
a3=dblarr(npt,npt)

for i=0,npt-1 do begin &$
 for j=0,npt-1 do begin &$
   a1(i,j)=randomu(seed)*rlev   &$
 endfor &$
endfor

for i=0,npt-1 do begin &$
 for j=0,npt-1 do begin &$
  a1(i,j)=a1(i,j)+3.d0*exp(-(i-40.d0)^2.d0/(2.d0^2.d0)-(j-20.d0)^2.d0/(4.d0^2.d0)) &$
 endfor &$
endfor


for i=0,npt-1 do begin &$
 for j=0,npt-1 do begin &$
  a1(i,j)=a1(i,j)+1.5d0*exp(-(i-80.d0)^2.d0/(4.d0^2.d0)-(j-60.d0)^2.d0/(6.d0^2.d0)) &$
 endfor &$
endfor


a1(randomu(seed)*npt,randomu(seed)*npt)=5.d0


for i=0,npt-1 do begin &$
 for j=0,npt-1 do begin &$
   a2(i,j)=randomu(seed)*rlev  &$
 endfor &$
endfor

for i=0,npt-1 do begin &$
 for j=0,npt-1 do begin &$
  a2(i,j)=a2(i,j)+3.d0*exp(-(i-35.d0)^2.d0/(3.d0^2.d0)-(j-40.d0)^2.d0/(4.d0^2.d0)) &$
 endfor &$
endfor


for i=0,npt-1 do begin &$
 for j=0,npt-1 do begin &$
  a2(i,j)=a2(i,j)+1.5d0*exp(-(i-70.d0)^2.d0/(4.d0^2.d0)-(j-50.d0)^2.d0/(6.d0^2.d0)) &$
 endfor &$
endfor


a2(randomu(seed)*npt,randomu(seed)*npt)=5.d0


for i=0,npt-1 do begin &$
 for j=0,npt-1 do begin &$
   a3(i,j)=randomu(seed)*rlev &$ 
 endfor &$
endfor

for i=0,npt-1 do begin &$
 for j=0,npt-1 do begin &$
  a3(i,j)=a3(i,j)+3.d0*exp(-(i-45.d0)^2.d0/(3.d0^2.d0)-(j-60.d0)^2.d0/(4.d0^2.d0)) &$
 endfor &$
endfor


for i=0,npt-1 do begin &$
 for j=0,npt-1 do begin &$
  a3(i,j)=a3(i,j)+1.5d0*exp(-(i-60.d0)^2.d0/(4.d0^2.d0)-(j-40.d0)^2.d0/(6.d0^2.d0)) &$
 endfor &$
endfor


a3(randomu(seed)*npt,randomu(seed)*npt)=5.d0








shiftx=0.d0
shifty=0.d0

x1=40
y1=20

!p.multi=[0,3,1]

tvframe,a1,/bar,charsize=2

ctrack,a1,a2,shiftx,shifty

print,shiftx,shifty

tvframe,a2,/bar,charsize=2

arrow,x1,y1,x1+shiftx,y1+shifty,/data,thick=2,color=255

xo=x1+shiftx
yo=y1+shifty

ctrack,a2,a3,shiftx,shifty

print,shiftx,shifty

tvframe,a3,/bar,charsize=2

arrow,x1,y1,xo,yo,/data,thick=2,color=255
arrow,xo,yo,xo+shiftx,yo+shifty,/data,thick=2,color=255

mult,1,1

end

