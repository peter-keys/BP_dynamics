; ------------------------------------------------------------------------------------------------------------
; 	GRAB MURAM MODEL AND CONVERT TO A NICOLE MODEL FOR SYNTHESIS TO STOKES I, Q, U, V
; ------------------------------------------------------------------------------------------------------------

; -----------------------------------------
;	GET THE MURAM MODEL
; -----------------------------------------
; >> This could be automated at a 
;    later time to get all the models
;    (Simple FOR-LOOP using a modlist)

;filesr=FILE_SEARCH('/data/espstore1/data/200G_mhd/result_*.595000')
;filese=FILE_SEARCH('/data/espstore1/data/200G_mhd/eos.595000')

ptk='/data/espstore1/data/200G_mhd/'
nmod='595000'

npt=480

coef=sqrt(4*3.14159265)

openr,1,ptk+'eos.'+nmod
a=assoc(1,fltarr(npt,200,npt))
temp=a(0)
pre=a(1)
close,1
print,tem(0,0,0)

openr,1,ptk+'result_0.'+nmod
a=assoc(1,fltarr(npt,200,npt))
rho=a(0)
close,1

openr,1,ptk+'result_1.'+nmod
a=assoc(1,fltarr(npt,200,npt))
vxx=a(0)
close,1

openr,1,ptk+'result_2.'+nmod
a=assoc(1,fltarr(npt,200,npt))
vyy=a(0)
close,1

openr,1,ptk+'result_3.'+nmod
a=assoc(1,fltarr(npt,200,npt))
vzz=a(0)
close,1

openr,1,ptk+'result_4.'+nmod
a=assoc(1,fltarr(npt,200,npt))
ene=a(0)
close,1

openr,1,ptk+'result_5.'+nmod
a=assoc(1,fltarr(npt,200,npt))
bxx=a(0)
close,1

openr,1,ptk+'result_6.'+nmod
a=assoc(1,fltarr(npt,200,npt))
byy=a(0)
close,1

openr,1,ptk+'result_7.'+nmod
a=assoc(1,fltarr(npt,200,npt))
bzz=a(0)
close,1


;openr,1,ptk+'eos.'+nmod
;a=assoc(1,fltarr(npt,200,npt))
;t2=a(0)
;pre=a(1)
;close,1
;npt=480
;openr,1,ptk+'eos.'+nmod
;a=assoc(1,fltarr(npt,100,npt))
;tem=a(0)
;close,1
;print,tem(0,0,0)


print,'data restored...'

vxx=vxx/rho
vyy=vyy/rho
vzz=vzz/rho
ene=ene/rho
bxx=bxx*coef
byy=byy*coef
bzz=bzz*coef
vhh=sqrt(vxx^2.d0+vzz^2.d0)

; >> Need the Z component for the synthesis.... This is 14Mm therefore create
;    an array that is 200 long and goes from -1400 to 1386 in steps of 14
z = (FINDGEN(200)*14)-1400

; -----------------------------------------
; VORTICITY (for future ref)
; -----------------------------------------
dvzdy=fltarr(npt,nptz,npt)
dvydz=fltarr(npt,nptz,npt)
dvxdz=fltarr(npt,nptz,npt)
dvzdx=fltarr(npt,nptz,npt)
dvydx=fltarr(npt,nptz,npt)
dvxdy=fltarr(npt,nptz,npt)

for i=0,npt-1 do begin
 for j=0,npt-1 do begin
   dvzdy(i,*,j)=deriv(vzz(i,*,j))
   dvxdy(i,*,j)=deriv(vxx(i,*,j))
 endfor
endfor

for i=0,npt-1 do begin
 for j=0,nptz-1 do begin
  dvydz(i,j,*)=deriv(vyy(i,j,*))
  dvxdz(i,j,*)=deriv(vxx(i,j,*))

  dvydx(*,j,i)=deriv(vyy(*,j,i))
  dvzdx(*,j,i)=deriv(vzz(*,j,i))
 endfor
endfor



wxx=dvzdy-dvydz
wyy=dvxdz-dvzdx
wzz=dvydx-dvxdy
  
  
endif


aswy=smooth(abs(wyy(*,nptz-1,*)),20,/edge_truncate)



; -----------------------------------------
; Generate the NICOLE Model
; -----------------------------------------
; > Put the MURAM sims into the format 
;   expected of a NICOLE model then 
;   save it out to be synthesised

dum = size(obs, /dim)
nx1 = npt
ny1 = npt
nz = N_ELEMENTS(temp[0,*,0])  

v_los = FLTARR(480,480,200)
FOR t = 0, 199 DO BEGIN &$
	FOR x = 0, 479 DO BEGIN &$
		FOR y = 0, 479 DO BEGIN &$
			;v_los[x,y,t] = SQRT(vxx[x,t,y]^2+vyy[x,t,y]^2+vzz[x,t,y]^2) &$
			v_los[x,y,t] = SQRT(vxx[x,t,y]^2+vyy[x,t,y]^2+vzz[x,t,y]^2) &$
		ENDFOR &$
	ENDFOR &$
ENDFOR

m1={z:fltarr(nx1,ny1,nz), tau:fltarr(nx1,ny1,nz),  t:fltarr(nx1,ny1,nz), $
   gas_p:fltarr(nx1,ny1,nz), rho:fltarr(nx1,ny1,nz), el_p:fltarr(nx1,ny1,nz), $
   v_los:fltarr(nx1,ny1,nz), v_mic:fltarr(nx1,ny1,nz), $
   b_los_z:fltarr(nx1,ny1,nz), b_los_x:fltarr(nx1,ny1,nz), $
   b_los_y:fltarr(nx1,ny1,nz), b_x:fltarr(nx1,ny1,nz), $
   b_y:fltarr(nx1,ny1,nz), b_z:fltarr(nx1,ny1,nz), $
   v_x:fltarr(nx1,ny1,nz), v_y:fltarr(nx1,ny1,nz), v_z:fltarr(nx1,ny1,nz), $
   nH: fltarr(nx1,ny1,nz), nHminus: fltarr(nx1,ny1,nz), nHplus: fltarr(nx1,ny1,nz),$
   nH2: fltarr(nx1,ny1,nz), nh2plus: fltarr(nx1,ny1,nz), $
   v_mac: fltarr(nx1,ny1), stray_frac: fltarr(nx1,ny1),$
   keep_el_p: fltarr(nx1,ny1), keep_gas_p: fltarr(nx1,ny1), $
   keep_rho: fltarr(nx1,ny1), keep_nH: fltarr(nx1,ny1), $
   keep_nHminus: fltarr(nx1,ny1), keep_nHplus: fltarr(nx1,ny1), $
   keep_nH2: fltarr(nx1,ny1), keep_nh2plus: fltarr(nx1,ny1), $
   ffactor: fltarr(nx1,ny1), abundance: fltarr(nx1,ny1,92)}

FOR ii = 0, nz-1 DO BEGIN &$
	m1.z[*,*,ii] = z[ii] &$
	m1.T[*,*,ii] = temp[*,ii,*] &$
	m1.gasp_p[*,*,ii] = pre[*,ii,*] &$
	m1.rho[*,*,ii] = rho[*,ii,*] &$
	m1.b_los_x[*,*,ii] = bxx[*,ii,*] &$
	m1.b_los_y[*,*,ii] = bzz[*,ii,*] &$
	m1.b_los_z[*,*,ii] = byy[*,ii,*] &$
	m1.v_los[*,*,ii] = vyy[*,ii,*] &$
	m1.v_mic[*,*,ii] = vhh[*,ii,*] &$

	ggg.b_los_x[xx,yy,*] = bxx[xx,nptz:nptz_1-1,yy]
	ggg.b_los_y[xx,yy,*] = bzz[xx,nptz:nptz_1-1,yy]
	ggg.b_los_z[xx,yy,*] = byy[xx,nptz:nptz_1-1,yy]
	ggg.rho[xx,yy,*] = rho[xx,nptz:nptz_1-1,yy]
	;ggg.el_p[xx,yy,*] = gas_p[xx,nptz:nptz_1-1,yy]
	ggg.v_los[xx,yy,*] = vyy[xx,nptz:nptz_1-1,yy]*(-1)
	ggg.v_mic[xx,yy,*] = vhh[xx,nptz:nptz_1-1,yy]
ENDFOR
vxx=vxx/rho
vyy=vyy/rho
vzz=vzz/rho
ene=ene/rho
bxx=bxx*coef
byy=byy*coef
bzz=bzz*coef

;; We smooth the model by skipping gridpoints vertically. This way we
;; remove a lot of the structure of the temperature profile
if(0) then begin &$ ; if(0) gives the buggy one 
   sk= 22 &$
   tau = reform(m.tau[0,0,*]) &$
   t = red_bezier3( reform(m.tau[0,0,0:*:sk]), reform(m.t[0,0,0:*:sk]), reform(m.tau[0,0,*]),/lin) &$
   el = exp(red_bezier3( reform(m.tau[0,0,0:*:sk]), alog(reform(m.el_p[0,0,0:*:sk])), reform(m.tau[0,0,*]),/lin)) &$

   for ii = 0,90 do begin &$
      m1.tau[*,*,ii] = tau[ii] &$
      m1.t[*,*,ii] = t[ii] &$
      m1.el_p[*,*,ii] = el[ii] &$
   endfor &$

endif else begin &$
   for ii = 0, nz-1 do begin &$   ; this one uses the previous 2d one that was already smoothed in height
      m1.tau[*,*,ii] = m.tau[0,0,ii] &$
      m1.t[*,*,ii] = m.t[0,0,ii] &$
      m1.el_p[*,*,ii] = m.el_p[0,0,ii] &$
      m1.v_mic[*,*,ii] = 0.0d5 &$
   endfor &$
endelse
   
;write_model2,folder+'modelin.mod',m1
idl_to_nicole, file='modelin.nic', m = m1 ; do m=read_model('modelin.nic'), do help,m,/st  and plot: m.T[0,0,*],  m.el_p[0,0,*],  like this plot,m.tau[0,0,*],m.T[0,0,*] 
Print,'Producing input model called: modelin.nic..'

