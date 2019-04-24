;;
;; Define input file and read one snapshot to invert

; 27Jul2014 QS 6302 Data:
;file = '/data/solarstore4/vh/SSTreduc/2014.07.27_6302dc_remod_2018/crispex/14:18:38/crispex.stokes.6302.14:18:38.time_corrected.fcube'
file = '/data/solarstore4/vh/SSTreduc/2014.07.27_6302dc_remod_2018/crispex/14:18:38/crispex.stokes.6302.14:18:38.time_corrected_cmapcorr.fcube'
file1 = '/data/solarstore4/vh/SSTreduc/2014.07.27_6302dc_remod_2018/crispex/14:18:38/crispex.stokes.6302.14:18:38.time_corrected_sp.fcube'
restore,'/data/solarstore4/vh/SSTreduc/2014.07.27_6302dc/crispex/14:18:38/wav.idl',/ver
;wav = f0('/data/solarstore4/vh/SSTreduc/2014.07.27_6302dc/crispex/14:18:38/wav.6302.f0')
;;;;<<<<< WHERE IS THIS FUNCTION?>>>>>>>>


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
   d = float(dat[5]) &$ ;; for example snapshot 10 		; USE 5 here as it is the 1st frame of this MBP (35) 
   free_lun, lun &$	
endif

; for 6301 
;   ;; Remove cross-talk from I - > Q,U,V
;   offsets = dblarr(4)
;   pos = [0,14]
;   for ii = 1, 3 do begin &$
;      offsets[ii] = median(d[*,*,pos,ii]/d[*,*,pos,0]) &$
;      d[*,*,*,ii] -= d[*,*,*,0] * offsets[ii] &$
;      print, 'Offsets Stk ',ii,' -> ', offsets[ii] &$
;   endfor
   
;for 6302
;   ;; Remove cross-talk from I - > Q,U,V
   offsets = dblarr(4)
   pos = [15,29]
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


dlam=18.5
iwav=iwav[16:31]

range = max(iwav) - min(iwav) + dlam*4 ;; extend the range including two extra points on each side
np = fix(range / float(dlam) + 1)   ;; number of points in the new grid
if(np/2 * 2 eq np) then np += 1     ;; Odd number makes more sense for convolutions

;; new grid
nwav = indgen(np) * dlam + (min(iwav) - dlam*2)

;; get the positions of the observed points in the new grid
; -> Vasco: idx determination fixed  
; As we are treating 6301 and 6302 separately (at first) 
; nw needs to be reset to the values of the wavelengths 
; for that particular scan position. This makes things 
; easier later in the file (where divisions by zero will pop up..)
nw = N_ELEMENTS(iwav)
idx = intarr(nw)
for ww = 0, (N_ELEMENTS(iwav))-1 do begin &$
      comparor=abs(nwav - iwav[ww]) &$
      idx[ww] = where(comparor eq min(comparor)) &$
endfor
print,'idx',idx
WRITEFITS,folder+'/idx_6302.fits',idx
WRITEFITS,folder+'/iwav_6302.fits',iwav

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
imean = fltarr(nw)
for ww = 0, nw -1 do imean[ww] = median(d[x0:x1, y0:y1, ww+16, 0])
;;
;; Fit parabola to line center to check if there is an offset
;;
dum = min(imean, p)
cc = parab_fit(wav[p-1:p+1], imean[p-1:p+1])

offset = cc[1] / cc[2] * 0.5d0


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
wav = wav[16:31] ;For the single line calc, 6301, otherwiese [16:31]!
;cw = 6301.5d0
cw = 6302.5d0 ;for the prep 6302 
;red_satlas, min(wav)+cw - 0.5, max(wav) + cw + 0.2, x, y ;for 6301
red_satlas, min(wav)+cw - 0.2, max(wav) + cw + 0.5, x, y ;for 6302
dum = min(y, p) ;; detect line center
x -= x[p]
dx = x[1] - x[0]

;; get CRISP profile in the grid of the SATLAS		; THIS NEEDS TO BE IN IDL NOT SSWIDL (SOMETHING WRONG WITH THE PROGRAMS IT LOADS)
fpi =  cfpi(cw)
n = 115 < n_elements(x)	; Original value -> had to change to get y and x the same size arrays
tw = (dindgen(n)-n/2) * dx
tr = fpi->dual_fpi(tw+cw)
tr /= total(tr) ;; normalize by the area
y = fftconvol(y, tr)
normalized_intensity = interpol(y,x,(wav[nw-1] - offset)) * int_ratio

;; Normalize observations
scale = normalized_intensity / imean[nw-1]
d *= scale
imean *= scale

;;
 ;; Select FOV to invert
;;

; MBP35 (get locations using 'MBP_barycentre_timeseries.bat')

x0 = 807
x1 = 846
y0 = 721
y1 = 759


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
   ;obs[xx-x0,yy-y0,*,ss] = interpol(reform(d[xx,yy,*,ss]), float(iwav), float(nwav))&$	; Original Change * in d[] to [0:15] for scan
    obs[xx-x0,yy-y0,*,ss] = interpol(reform(d[xx,yy,16:31,ss]), float(iwav), float(nwav)) &$
endfor



;;
;; write interpolated profiles to disk		
;;
;;;with rebinning
;;;obs = (rebin(obs[*,1:*,*,*], [nx1/2, ny1/2, np, 4]))[*,*,*,*]
idl_to_nicole, file=folder+'/observed2.nic', i=obs[*,*,*,0], q = obs[*,*,*,1], u = obs[*,*,*,2], v = obs[*,*,*,3]



;;
;; Create instrumental profile, asumming Jaime's modifications to NICOLE
;;
fpi= cfpi(cw)
fwhm = fpi->getwidth() * 6 * 1000;; generate the profile for a range of x3 the FWHM
npt = round(fwhm / float(dlam))
if(npt/2 * 2 eq npt) then npt -= 1 ;; odd number
tw = (dindgen(npt) - npt/2) * dlam * 1.d-3
tr = fpi->dual_fpi(tw + cw)
tr /= total(tr, /double)

; Save it to file, first entry is the header, same size as the
; transmission profile 
openw,1,'Instrumental_profile.dat'
kk = dblarr(np)
kk[0]  = [npt]
writeu,1,kk
kk[0:npt-1] = shift(tr, -npt/2)
for yy = 0, ny1-1 do for xx=0, nx1-1 do writeu, 1, kk
close,1


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
weights = dblarr(5, np) + 1.d13 ;; set all values by default to weight zero and then modify the observed ones
weights[0,*] = nwav * 1.0d-3
weights[1,idx] = 2.e-3 ;; assume values relative to the continuum intensity.
weights[2,idx] = 1.e-3
weights[3,idx] = 1.e-3
weights[4,idx] = 3.e-3
; >> Mask the Telluric....
weights[1:4,idx[9:14]] = 1.d13
;; Save to file
openw, 1, 'Weights.pro', width = 300 &$
printf, 1, weights &$
close, 1



;;
;; Init guessed model, typically something like a smoothed FALC will work
;;

;; read falc from file (1 column)
m = read_model('/data/solarstore3/phk/nicole_2018/NICOLE/run/modelin_falc.nic')
dim = size(m.t,/dim)
m1 = create_model2(nx1,ny1,91)

;; We smooth the model by skipping gridpoints vertically. This way we
;; remove a lot of the structure of the temperature profile
sk= 2
tau = m.tau[*]
t = bezier( reform(m.tau[0,0,0:*:sk]), reform(m.t[0,0,0:*:sk]), reform(m.tau[0,0,*]),/lin)
el = exp(bezier( reform(m.tau[0,0,0:*:sk]), alog(reform(m.el_p[0,0,0:*:sk])), reform(m.tau[0,0,*]),/lin))


;for ii = 0,101 do begin &$		;Original (may not match the size of the other arrays)
for ii = 0,90 do begin &$		;new
   m1.tau[*,*,ii] = tau[ii] &$
   m1.t[*,*,ii] = t[ii] &$
   m1.el_p[*,*,ii] = el[ii] &$
endfor
; >> Just copy the output model from 6301 run to 6302 i.e, cp modelout.mod modelin.mod
;write_model2,'/data/solarstore3/phk/nicole_2018/NICOLE/run/modelin.mod',m1


;;
;; Print details for the NICOLE.input
;;
;; >> Peter you twat you didn't change the 
print, 'Wavelength Start =', nwav[0]*1.d-3 - offset + 6302.4936d0 ;; Assume that line center is at rest (can't do much better in the chromosphere). Use the lab wavelength in LINES
print, 'nwavelength=', np
print, 'dlambda=', string(dlam)+'mA'
end

