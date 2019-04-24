;;	Prep the data to run inversions on both the 6301/6302 lines together
;;
;; Define input file and read one snapshot to invert

; 27Jul2014 QS 6302 Data:
;file = '/data/solarstore4/vh/SSTreduc/2014.07.27_6302dc_remod_2018/crispex/14:18:38/crispex.stokes.6302.14:18:38.time_corrected.fcube'
file = '/data/solarstore4/vh/SSTreduc/2014.07.27_6302dc_remod_2018/crispex/14:18:38/crispex.stokes.6302.14:18:38.time_corrected_cmapcorr.fcube'
file1 = '/data/solarstore4/vh/SSTreduc/2014.07.27_6302dc_remod_2018/crispex/14:18:38/crispex.stokes.6302.14:18:38.time_corrected_sp.fcube'
restore,'/data/solarstore4/vh/SSTreduc/2014.07.27_6302dc/crispex/14:18:38/wav.idl',/ver
;wav = f0('/data/solarstore4/vh/SSTreduc/2014.07.27_6302dc/crispex/14:18:38/wav.6302.f0')	; f0 requires just IDL (not ssw) and ANA enabled

; SET OUTPUT HERE
folder='.'

;;
;; Convert the wavelength array to mA and store as integer
;; Much cleaner to compare integer numbers (further down)
;;
iwav = fix(wav * 1000.d0)


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
   d = float(dat[5]) &$ ;; for example snapshot 10 		; 30 is 174 USE 5 here as it is the 1st frame of this MBP (35) ; frame 21 (-2) is the frame with the issues...
   free_lun, lun &$	
endif

; for 6301 
   ;; Remove cross-talk from I - > Q,U,V
;   offsets = dblarr(4)
;   pos = [0,14]
;   for ii = 1, 3 do begin &$
;      offsets[ii] = median(d[*,*,pos,ii]/d[*,*,pos,0]) &$
;      d[*,*,*,ii] -= d[*,*,*,0] * offsets[ii] &$
;      print, 'Offsets Stk ',ii,' -> ', offsets[ii] &$
;   endfor
   
;for 6302
;   ;; Remove cross-talk from I - > Q,U,V
;   offsets = dblarr(4)
;   pos = [15,29]
;   for ii = 1, 3 do begin &$
;      offsets[ii] = median(d[*,*,pos,ii]/d[*,*,pos,0]) &$
;      d[*,*,*,ii] -= d[*,*,*,0] * offsets[ii] &$
;      print, 'Offsets Stk ',ii,' -> ', offsets[ii] &$
;   endfor

;for 6301 & 6302 Together
;   ;; Remove cross-talk from I - > Q,U,V
   offsets = dblarr(4)
   pos = [0,14,15,29]
   for ii = 1, 3 do begin &$
      offsets[ii] = median(d[*,*,pos,ii]/d[*,*,pos,0]) &$
      d[*,*,*,ii] -= d[*,*,*,0] * offsets[ii] &$
      print, 'Offsets Stk ',ii,' -> ', offsets[ii] &$
   endfor


;;
;; The spectra are not critically sampled. For convolutions with the
;; CRISP profile this is important. We must select a REGULAR grid of
;; wavelength that contain our observed points and that is, at least
;; 0.5 * FWHM of CRISP at a given wavelength. In this case the FWHM of
;; CRISP is ~55 mA. Since these observations are taken in multiples of
;; 35 mA, I choose that wavelength separation.

;;
;; To avoid border effect in the convolutions at the outermost points
;; of the spectra, it is good to extend the grid a bit more than our
;; observed range. In the end, all these extra points will have
;; Weight = 0 except for the observed ones, but we need then to
;; perform the convolution.
;;
;; Grid for 6301 assuming half sampling: 19.5
;; Grid for 6302 assuming half sampling: 18.5 

dlam1=19.5
dlam2=18.5
;iwav=iwav[16:31]
range1 = max(iwav[0:15]) - min(iwav[0:15]) + dlam1*4
range2 = max(iwav[16:31]) - min(iwav[16:31]) + dlam2*4

;range = max(iwav) - min(iwav) + dlam*4 ;; extend the range including two extra points on each side
;np = fix(range / float(dlam) + 1)   ;; number of points in the new grid
;if(np/2 * 2 eq np) then np += 1     ;; Odd number makes more sense for convolutions

np1 = fix(range1 / float(dlam1) + 1)
np2 = fix(range2 / float(dlam2) + 1)
if(np1/2 * 2 eq np1) then np1 += 1     ;; Odd number makes more sense for convolutions
if(np2/2 * 2 eq np2) then np2 += 1     ;; Odd number makes more sense for convolutions

; New Grid. Make each individually then add them into 1 grid
nwav1 = indgen(np1) * dlam1 + (min(iwav[0:15]) - dlam1*2)
nwav2 = indgen(np2) * dlam2 + (min(iwav[16:31]) - dlam2*2)
nwav = DBLARR(N_ELEMENTS(nwav1)+N_ELEMENTS(nwav2))
nwav[0:N_ELEMENTS(nwav1)-1] = nwav1
nwav[N_ELEMENTS(nwav1):N_ELEMENTS(nwav)-1] = nwav2

np = np1 + np2
;; new grid
;nwav = indgen(np) * dlam + (min(iwav) - dlam*2)

;; get the positions of the observed points in the new grid
; -> Vasco: idx determination fixed  
; As we are treating 6301 and 6302 separately (at first) 
; nw needs to be reset to the values of the wavelengths 
; for that particular scan position. This makes things 
; easier later in the file (where divisions by zero will pop up..)
iwav1 = iwav[0:15]
iwav2 = iwav[16:*]
nw = N_ELEMENTS(iwav)
idx1 = intarr(N_ELEMENTS(iwav1))
for ww = 0, (N_ELEMENTS(iwav1))-1 do begin &$
      comparor=abs(nwav1 - iwav1[ww]) &$
      idx1[ww] = where(comparor eq min(comparor)) &$
endfor
print,'idx_s1',idx1

idx2 = intarr(N_ELEMENTS(iwav1))
for ww = 0, (N_ELEMENTS(iwav2))-1 do begin &$
      comparor=abs(nwav2 - iwav2[ww]) &$
      idx2[ww] = where(comparor eq min(comparor)) &$
endfor
print,'idx_s2',idx2

idx = intarr(nw)
for ww = 0, (N_ELEMENTS(iwav))-1 do begin &$
      comparor=abs(nwav - iwav[ww]) &$
      idx[ww] = where(comparor eq min(comparor)) &$
endfor

Print,' The following compare the interpolated grid to the real wavelengths'
Print,' They should be close (+/- 2 or 3), but to not have to 1-to-1'
Print,' They are important for the Weights.pro file that is generated'
FOR i = 0, N_ELEMENTS(iwav)-1 DO print,'Nwav[idx[i]] = ', nwav[idx[i]],'  iWav[i] = ', iwav[i]


; NB idx is the idx of 6301 and 6302 combined, therefore, later you need the individual idx1 and idx2 arrays to get the right value
WRITEFITS,folder+'/idx_6301_6302.fits',idx
WRITEFITS,folder+'/iwav_6301_6302.fits',iwav

;;
;; Select a portion of quiet-Sun (as quiet as possible)
;; to compute a quiet-Sun average profile. We will use it
;; to estimate the continuum intensity of our dataset
;;
;x0 = 50
;x1 = 380
;y0 = 40
;y1 = 255
x0 = 140
x1 = 390
y0 = 100
y1 = 390
if(n_elements(imean) eq 0) then red_show, histo_opt(d[x0:x1, y0:y1, 0, 0]) ;; show for inspection the first time

;;
;; Compute spatial average. Use the median to remove outliers
;; from plage and other stuff
;;

imean1 = FLTARR(N_ELEMENTS(wav[0:15]))
imean2 = FLTARR(N_ELEMENTS(wav[16:*]))
for ww = 0, 15 do imean1[ww] = median(d[x0:x1, y0:y1, ww, 0])
for ww = 16, 31 do imean2[ww-16] = median(d[x0:x1, y0:y1, ww, 0])

imean = FLTARR(nw)
for ww = 0, nw-1 do imean[ww] = median(d[x0:x1, y0:y1, ww, 0])

print,'Plotting the mean profiles....'
window,0,xsize=600,ysize=600
plot,iwav1,imean1,xtitle='Wavelength from lamda0 (mA)',Title='Median Observed Profile 6301',thick=2 
wait,5
plot,iwav2,imean2,xtitle='Wavelength from lamda0 (mA)',Title='Median Observed Profile 6302',thick=2 
wait,5
plot,iwav,imean,xtitle='Wavelength from lamda0 (mA)',Title='Median Observed Profile 6301/6302',thick=2 
wait,5

;;
;; Fit parabola to line center to check if there is an offset
;;
wav1 = wav[0:15]
wav2 = wav[16:*]
dum1 = min(imean1, p)
cc1 = parab_fit(wav1[p-1:p+1], imean1[p-1:p+1])
offset1 = cc1[1] / cc1[2] * 0.5d0

dum2 = min(imean2, p)
cc2 = parab_fit(wav2[p-1:p+1], imean2[p-1:p+1])
offset2 = cc2[1] / cc2[2] * 0.5d0


;;
;; Use the same technique as de la Cruz Rodriguez et al. (2013) to
;; derive the center to limb variations: use a simulation and then
;; multiply that ratio with the intensity of the solar atlas at the
;; outermost point of the observations. The profiles from the
;; simulation are included here as .sav files for many lines.
;;
mu_obs = 1 ;; the heliocentric angle of the observations
int_ratio = 1 ; because mu=1 

;;
;; Get the solar atlas and get the intensity at the outermost
;; observed point. Then apply the center to limb variation
;;

;; FIRST LINE
nw1 = N_ELEMENTS(iwav1)
cw1 = 6301.5d0
red_satlas, min(wav1)+cw1 - 0.5, max(wav1) + cw1 + 0.2, x, y ;for 6301
dum = min(y, p) ;; detect line center
x -= x[p]
dx = x[1] - x[0]

;; get CRISP profile in the grid of the SATLAS
fpi =  cfpi(cw1)
n = 115 < n_elements(x)	; Original value -> had to change to get y and x the same size arrays
tw = (dindgen(n)-n/2) * dx
tr = fpi->dual_fpi(tw+cw1)
tr /= total(tr) ;; normalize by the area
y = fftconvol(y, tr)
normalized_intensity1 = interpol(y,x,(wav1[nw1-1] - offset1)) * int_ratio

;; Normalize observations
scale1 = normalized_intensity1 / imean1[nw1-1]
;d *= scale
d[*,*,0:15,*] = d[*,*,0:15,*] * scale1
imean1 *= scale1

;; SECOND LINE
nw2 = N_ELEMENTS(iwav2)
cw2 = 6302.5d0 ;for the prep 6302 
red_satlas, min(wav)+cw2 - 0.2, max(wav) + cw2 + 0.5, x, y ;for 6302
dum = min(y, p) ;; detect line center
x -= x[p]
dx = x[1] - x[0]

;; get CRISP profile in the grid of the SATLAS
fpi =  cfpi(cw2)
n = 115 < n_elements(x)	; Original value -> had to change to get y and x the same size arrays
tw = (dindgen(n)-n/2) * dx
tr = fpi->dual_fpi(tw+cw2)
tr /= total(tr) ;; normalize by the area
y = fftconvol(y, tr)
normalized_intensity2 = interpol(y,x,(wav2[nw2-1] - offset2)) * int_ratio

;; Normalize observations
scale2 = normalized_intensity2 / imean2[nw2-1]
;d *= scale
d[*,*,16:*,*] = d[*,*,16:*,*] * scale2
imean2 *= scale2

print,'Plotting the normalised mean profiles (check to make sure this step is correct)....'
window,0,xsize=600,ysize=600
plot,iwav1,imean1,xtitle='Wavelength from lamda0 (mA)',Title='Normalised Observed Profile 6301',thick=2
wait,5
plot,iwav2,imean2,xtitle='Wavelength from lamda0 (mA)',Title='Normalised Observed Profile 6302',thick=2
wait,5
plot,iwav,imean,xtitle='Wavelength from lamda0 (mA)',Title='Normalised Observed Profile 6301/6302',thick=2
wait,5

;;
 ;; Select FOV to invert
;;

; MBP35 (get locations using 'MBP_barycentre_timeseries.bat')

;x0 = 807
;x1 = 846
;y0 = 721
;y1 = 759

;x0 = 772
;x1 = 881
;y0 = 686
;y1 = 794

;Full FoV
;174:
;x0 = 35
;x1 = 950
;y0 = 35
;y1 = 950
;x0 = 72
;x1 = 890
;y0 = 70
;y1 = 935
x0 = 410
x1 = 927
y0 = 434
y1 = 935


nx1 = x1 - x0 + 1
ny1 = y1 - y0 + 1


;; Interpolate new line positions. It does not matter which
;; interpolant is used as long as the curve does go over the
;; original data. Just pick something fast, the weight will be zero
;; anyhow
;; Like n (i.e., with x and y not being same dimensions) ss seems to 
;; be an issue here due to the change of not taking the entire spectrum
;; (i.e., treating 6301 & 6302 separately initially). 
obs  = fltarr(nx1, ny1, np, 4)
for yy = y0, y1 do for xx = x0, x1 do for ss = 0, 3 do begin &$
    obs[xx-x0,yy-y0,*,ss] = interpol(reform(d[xx,yy,*,ss]), float(iwav), float(nwav)) &$
endfor

print,'Potting the FoV to be inverted...'
tvim,obs[*,*,0,0],title='FoV that is to be Inverted (Intensity)',/sc
wait,5
;; Show the Stokes V profile to make sure that there is signal there
;; Use the centre of the 6302 scan (i.e., number of points interpolated points in 
;; 6301 (np1) + 1/2 number of points in 6302 (np2)
tvim,obs[*,*,(FIX(np1 + (np2/2.))),3],title='FoV that is to be Inverted (Stokes-V)',/sc
wait,5

;;
;; write interpolated profiles to disk		
;;
;;;with rebinning
;;;obs = (rebin(obs[*,1:*,*,*], [nx1/2, ny1/2, np, 4]))[*,*,*,*]
idl_to_nicole, file=folder+'/observed_FOVq2.nic', i=obs[*,*,*,0], q = obs[*,*,*,1], u = obs[*,*,*,2], v = obs[*,*,*,3]

print,'Saving observed profiles to '+folder+'/observed_FOV.nic'

;;
;; Create instrumental profile, asumming Jaime's modifications to NICOLE
;;
; >> npt requires value for dlam, which I don't have here (as it's 2 lines - dlam1 & dlam2)

; >> cw needs to change above as it is the central wavelength and we have 2 here...
fpi= cfpi(cw1)
fwhm1 = fpi->getwidth() * 6 * 1000;; generate the profile for a range of x3 the FWHM
npt1 = round(fwhm1 / float(dlam1))
if(npt1/2 * 2 eq npt1) then npt1 -= 1 ;; odd number
tw1 = (dindgen(npt1) - npt1/2) * dlam1 * 1.d-3
tr1 = fpi->dual_fpi(tw1 + cw1)
tr1 /= total(tr1, /double)

fpi= cfpi(cw2)
fwhm2 = fpi->getwidth() * 6 * 1000;; generate the profile for a range of x3 the FWHM
npt2 = round(fwhm2 / float(dlam2))
if(npt2/2 * 2 eq npt2) then npt2 -= 1 ;; odd number
tw2 = (dindgen(npt2) - npt2/2) * dlam2 * 1.d-3
tr2 = fpi->dual_fpi(tw2 + cw2)
tr2 /= total(tr2, /double)


; Save it to file, first entry is the header, same size as the
; transmission profile 

; >> Need to get the two profiles in the one Instrumental profile file...
;; Becky's method states.... (for reference from getwav_sept.pro)
;;
;; Transmission profile, using a theoretical CRISP profile. Since each
;; region has a different wavelength step, we need to create one
;; instrumental profile per region, even though the profile shape is
;; the same.
;; 

openw,1,'Instrumental_profile.dat'
kk = dblarr(np)
kk[0:1]  = [npt1,npt2]
writeu,1,kk
kk[0:npt1-1] = shift(tr1, -npt1/2)
kk[npt1:npt1+npt2-1] = shift(tr2, -npt2/2)
for yy = 0, ny1-1 do for xx=0, nx1-1 do writeu, 1, kk
close,1

;; Write the instrumental profile of both regions as many times as
;; pixel in the FOV.
;for yy = 0, ny-1 do for xx = 0, nx-1 do writeu, lun, rec
;free_lun, lun

;;
;; Create weights file. In reality we should give an estimate of the
;; noise, so large values means low weight.
;; there are 5 columns due to legacy input files (lambda, weight_I,
;; weight_Q, Weight_U, Weight_V). The lambda column is dummy, not used
;; but it is good practice to store the wavelength array in case we
;; need to check the values later.
;; Note that giving large weights to Q,U,V can affect the temperature
;; maps with a lot of noise
;;
;;weights[1:4,idx[16+9:16+14]] = 1.d13 ;; mask the points affecte by telluric
;; I'd personally get rid of the points 11:14.. for the telluric


; >> Might need to separate these into 2 separate weights files then add them together after
; >> Use the idx1 and idx2 values to set the weights for both lines
; >> Not needed idx (i.e., the combined one) is fine for what we need.
weights = dblarr(5, np) + 1.d13 ;; set all values by default to weight zero and then modify the observed ones
weights[0,*] = nwav * 1.0d-3
weights[1,idx] = 2.e-3 ;; assume values relative to the continuum intensity.
weights[2,idx] = 1.d13			; Aaron thinks that I should comlpetely ignore Q & u as I don't want Bx By
weights[3,idx] = 1.d13
weights[4,idx] = 4.e-3
;weights[1:4,idx[16+9:16+14]] = 1.d13
Weights[1:4,idx[28:30]] = 1.d13		;Aaron thinks the mask of the telluric should be minimised to get red wing info

;; Save to file
openw, 1, 'Weights.pro', width = 300 &$
printf, 1, weights &$
close, 1
WRITEFITS,folder+'/nwav_6301_6302.fits',nwav

;;
;; Init guessed model, typically something like a smoothed FALC will work (Older not correct version...)
;;

;; read falc from file (1 column)
;m = read_model('/data/solarstore3/phk/nicole_2018/NICOLE/run/modelin_falc.nic')
;dim = size(m.t,/dim)
;m1 = create_model2(nx1,ny1,91)

;; We smooth the model by skipping gridpoints vertically. This way we
;; remove a lot of the structure of the temperature profile
;sk= 2
;tau = m.tau[*]
;t = bezier( reform(m.tau[0,0,0:*:sk]), reform(m.t[0,0,0:*:sk]), reform(m.tau[0,0,*]),/lin)
;el = exp(bezier( reform(m.tau[0,0,0:*:sk]), alog(reform(m.el_p[0,0,0:*:sk])), reform(m.tau[0,0,*]),/lin))

; >> Aaron says that my Pgas initial value is too low
; >> Also, the 
;for ii = 0,101 do begin &$		;Original (may not match the size of the other arrays)
;for ii = 0,90 do begin &$		;new
;   m1.tau[*,*,ii] = tau[ii] &$
;   m1.t[*,*,ii] = t[ii] &$
;   m1.el_p[*,*,ii] = el[ii] &$
;endfor
;write_model2,'/data/solarstore3/phk/nicole_2018/NICOLE/run/modelin.mod',m1

;;
;; Init guessed model, typically something like a smoothed FALC will work
;;

;; read falc from file (1 column)
m = read_model('/data/solarstore3/phk/nicole_2018/NICOLE/run/modelin_falc.nic') ; doesn't get used unless the if down below turns to if(1)
dim = size(m.t,/dim)


;m1 = create_model2(nx1,ny1,91)
dum = size(obs, /dim)
nx1 = dum[0]
ny1 = dum[1]
nz = dim[2]  

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
idl_to_nicole, file='modelin_fovq2.nic', m = m1 ; do m=read_model('modelin.nic'), do help,m,/st  and plot: m.T[0,0,*],  m.el_p[0,0,*],  like this plot,m.tau[0,0,*],m.T[0,0,*] 
Print,'Producing input model called: modelin.nic..'

;;
;; Print details for the NICOLE.input
;;
Print,' '
Print,'Deatails needed for NICOLE.input file...'
print,' '
print, 'Region 1'
print, 'Wavelength step = ', dlam1, ' mA'
print, 'First wavelength = ', nwav1[0]*1.d-3+6301.5012d0-offset1
print, 'number of wavelengths = ', np1
print, ' '
print, 'Region 2'
print, 'Wavelength step = ', dlam2, ' mA'
print, 'First wavelength = ', nwav2[0]*1.d-3+6302.4936d0-offset2
print, 'number of wavelengths = ', np2
print,' '


end

