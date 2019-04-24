;; ----------------------------------------------------------------------------
	;;		PREP PROFILES FOR SIR INVERSION (ALL MBP PIXELS)
;; ----------------------------------------------------------------------------

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
;; Read in the MBP file to get relevant info
;;

;readcol,'/data/solarstore3/phk/data/27Jul2014/final_tracking_file',f='D,D,D,D,D,x', btpt_num, frame, sxloc, yloc, area
;xcoords = xloc[WHERE (btpt_num eq MBP)]  
;ycoords = yloc[WHERE (btpt_num eq MBP)]   				
;frame_num = frame[WHERE(btpt_num eq MBP)]  
;initial_frame = frame_num[0] 
;final_frame = frame_num[N_ELEMENTS(frame_num) - 1]   
;n_frames = N_ELEMENTS(frame_num) 

restore,'/data/solarstore3/phk/data/27Jul2014/SST_6302_scans_trimmed_to_MBP_detections.sav',/ver
;restore,'/data/solarstore3/phk/data/27Jul2014/final_tracked_mbps.sav',/ver
dims = SIZE(sst_scans,/dim)

	; > Correct values for trimming	
;xcoords = xcoords + 60. 
;ycoords = ycoords + 65. 
;frame_num = frame_num + 3. 
;initial_frame = initial_frame + 3. 
;final_frame = final_frame +3. 

;; >> Start a FOR LOOP to go through the frames of the MBP
;;	- Work out Cross Talk
;;	- Work out the scaling for each frame to nromalise
;;	- Work out the offsets for each profile
;;	- Output the straylight profiles
;;	- Output the Observed corrected profiles
;;	- Output the Observed profiles (with weighting apploed to unobserved interpolations)

dnew = FLTARR(dims[0],dims[1],92,dims[3])
d = FLTARR(dims[0],dims[1],dims[4],dims[3])

FOR fr = 0, dims[2]-1 DO BEGIN &$
	; Load in the image
	FOR s = 0, dims[4]-1 DO d[*,*,s,*] = sst_scans[*,*,fr,*,s] &$	
	
	;for 6301 & 6302 Together
	;   ;; Remove cross-talk from I - > Q,U,V
	   offsets = dblarr(4) &$
	   pos = [0,14,15,29] &$
	   FOR ii = 1, 3 do begin &$
	      offsets[ii] = median(d[*,*,pos,ii]/d[*,*,pos,0]) &$
	      d[*,*,*,ii] -= d[*,*,*,0] * offsets[ii] &$
	      print, 'Offsets Stk ',ii,' -> ', offsets[ii] &$
	   ENDFOR &$

	;; >>>	CHANGE THE SAMPLING FOR THE GRID HERE <<< ;;
	;; Grid for 6301 assuming half sampling: 19.25
	;; Grid for 6302 assuming half sampling: 18.5 

	dlam1=19.25 &$
	dlam2=18.5 &$
	range1 = max(iwav[0:15]) - min(iwav[0:15]) + dlam1*4 &$
	range2 = max(iwav[16:31]) - min(iwav[16:31]) + dlam2*4 &$
	np1 = fix(range1 / float(dlam1) + 1) &$
	np2 = fix(range2 / float(dlam2) + 1)&$
	if(np1/2 * 2 eq np1) then np1 += 1 &$    ;; Odd number makes more sense for convolutions
	if(np2/2 * 2 eq np2) then np2 += 1 &$    ;; Odd number makes more sense for convolutions

	; New Grid. Make each individually then add them into 1 grid
	nwav1 = indgen(np1) * dlam1 + (min(iwav[0:15]) - dlam1*2) &$
	nwav2 = indgen(np2) * dlam2 + (min(iwav[16:31]) - dlam2*2) &$
	nwav = DBLARR(N_ELEMENTS(nwav1)+N_ELEMENTS(nwav2))
	nwav[0:N_ELEMENTS(nwav1)-1] = nwav1 &$
	nwav[N_ELEMENTS(nwav1):N_ELEMENTS(nwav)-1] = nwav2 &$

	; Find the positions in the new grid which are observed / interpolated
	np = np1 + np2 &$
	iwav1 = iwav[0:15] &$
	iwav2 = iwav[16:*] &$
	nw = N_ELEMENTS(iwav) &$
	idx1 = intarr(N_ELEMENTS(iwav1)) &$
	FOR ww = 0, (N_ELEMENTS(iwav1))-1 do begin &$
	      comparor=abs(nwav1 - iwav1[ww]) &$
	      idx1[ww] = where(comparor eq min(comparor)) &$
	ENDFOR &$
	idx2 = intarr(N_ELEMENTS(iwav1)) &$
	FOR ww = 0, (N_ELEMENTS(iwav2))-1 do begin &$
	      comparor=abs(nwav2 - iwav2[ww]) &$
	      idx2[ww] = where(comparor eq min(comparor)) &$
	ENDFOR &$
	idx = intarr(nw) &$
	FOR ww = 0, (N_ELEMENTS(iwav))-1 do begin &$
	      comparor=abs(nwav - iwav[ww]) &$
	      idx[ww] = where(comparor eq min(comparor)) &$
	ENDFOR &$

	; For checking the difference between the observed wavelength and the closest interpolated grid posiiton
	difflam = DBLARR(N_ELEMENTS(iwav)) &$				
	FOR i = 0, N_ELEMENTS(iwav)-1 DO difflam[i] = nwav[idx[i]]-iwav[i] &$
	print,' ' &$
	Print,'The difference between the observed position and the closest interpolated grid position is:' &$
	print,difflam &$

	; NB idx is the idx of 6301 and 6302 combined, therefore, later you need the individual idx1 and idx2 arrays to get the right value
	SAVE,FILENAME=folder+'/wave_grids.sav',iwav,iwav1,iwav2,nwav,nwav1,nwav2,idx,idx1,idx2

	;;
	;; Select a portion of quiet-Sun (as quiet as possible)
	;; to compute a quiet-Sun average profile. We will use it
	;; to estimate the continuum intensity of our dataset
	;;
	x0 = 140 &$
	x1 = 390 &$
	y0 = 100 &$
	y1 = 390 &$

	;;
	;; Compute spatial average. Use the median to remove outliers
	;; from plage and other stuff
	;;

	imean1 = FLTARR(N_ELEMENTS(wav[0:15])) &$
	imean2 = FLTARR(N_ELEMENTS(wav[16:*]))  &$
	for ww = 0, 15 do imean1[ww] = median(d[x0:x1, y0:y1, ww, 0]) &$
	for ww = 16, 31 do imean2[ww-16] = median(d[x0:x1, y0:y1, ww, 0]) &$

	imean = FLTARR(nw) &$
	for ww = 0, nw-1 do imean[ww] = median(d[x0:x1, y0:y1, ww, 0]) &$

	;;
	;; Fit parabola to line center to check if there is an offset
	;;
	wav1 = wav[0:15] &$
	wav2 = wav[16:*] &$
	dum1 = min(imean1, p) &$
	cc1 = parab_fit(wav1[p-1:p+1], imean1[p-1:p+1]) &$
	offset1 = cc1[1] / cc1[2] * 0.5d0 &$

	dum2 = min(imean2, p) &$
	cc2 = parab_fit(wav2[p-1:p+1], imean2[p-1:p+1]) &$
	offset2 = cc2[1] / cc2[2] * 0.5d0 &$


	mu_obs = 1 &$ ;; the heliocentric angle of the observations
	int_ratio = 1 &$ ; because mu=1  

	;;
	;; Get the solar atlas and get the intensity at the outermost
	;; observed point. Then apply the center to limb variation
	;;

	;; FIRST LINE
	nw1 = N_ELEMENTS(iwav1) &$
	cw1 = 6301.5d0 &$
	red_satlas, min(wav1)+cw1 - 0.5, max(wav1) + cw1 + 0.2, x, y  &$;for 6301
	dum = min(y, p)  &$;; detect line center
	x -= x[p] &$
	dx = x[1] - x[0] &$

	;; get CRISP profile in the grid of the SATLAS
	fpi =  cfpi(cw1) &$
	n = 115 < n_elements(x) &$	; Original value -> had to change to get y and x the same size arrays
	tw = (dindgen(n)-n/2) * dx &$
	tr = fpi->dual_fpi(tw+cw1) &$
	tr /= total(tr) &$ ;; normalize by the area
	y = fftconvol(y, tr) &$
	normalized_intensity1 = interpol(y,x,(wav1[nw1-1] - offset1)) * int_ratio &$

	;; Normalize observations
	scale1 = normalized_intensity1 / imean1[nw1-1] &$
	d[*,*,0:15,*] = d[*,*,0:15,*] * scale1 &$
	imean1 *= scale1 &$

	;; SECOND LINE
	nw2 = N_ELEMENTS(iwav2) &$
	cw2 = 6302.5d0  &$ ;for the prep 6302 
	red_satlas, min(wav)+cw2 - 0.2, max(wav) + cw2 + 0.5, x, y  &$ ;for 6302
	dum = min(y, p)  &$ ;; detect line center
	x -= x[p] &$
	dx = x[1] - x[0]  &$

	;; get CRISP profile in the grid of the SATLAS
	fpi =  cfpi(cw2) &$
	n = 115 < n_elements(x)	 &$ ; Original value -> had to change to get y and x the same size arrays
	tw = (dindgen(n)-n/2) * dx &$
	tr = fpi->dual_fpi(tw+cw2) &$
	tr /= total(tr)  &$ ;; normalize by the area
	y = fftconvol(y, tr) &$ 
	normalized_intensity2 = interpol(y,x,(wav2[nw2-1] - offset2)) * int_ratio &$

	;; Normalize observations
	scale2 = normalized_intensity2 / imean2[nw2-1] &$
	;d *= scale
	d[*,*,16:*,*] = d[*,*,16:*,*] * scale2 &$
	imean2 *= scale2 &$

	;; Compare selected positions to the Observed lines from the atlas.
	;Full line scan
	nw1 = N_ELEMENTS(iwav1) &$
	cw1 = 6301.5d0 &$
	nw2 = N_ELEMENTS(iwav2) &$
	cw2 = 6302.5d0  &$

	red_satlas,min(wav1)+cw1 - 0.5,max(wav) + cw2 + 0.5,xfull,yfull &$
	;dnew = FLTARR(scan_dims[0],scan_dims[1],scan_dims[4],scan_dims[3],scan_dims[2])
	;d = FLTARR(scan_dims[0],scan_dims[1],scan_dims[4],scan_dims[3])
	; d = x, y, s, stk
	;dnew = x, y, s, stk, fr
	
	FOR xx = 0, dims[0]-1 DO BEGIN &$
	FOR yy = 0, dims[0]-1 DO BEGIN &$
	FOR ss = 0, 3 DO BEGIN &$
        	dnew[xx,yy,*,ss] = interpol(reform(d[xx,yy,*,ss]), float(iwav), float(nwav)) &$
    	ENDFOR &$
	ENDFOR &$
    	ENDFOR &$

	SAVE,FILENAME='/data/solarstore3/phk/data/27Jul2014/Inversions/SIR/Interpolated_profile'+STRING(fr,FORMAT='(I3.3)')+'.sav',dnew &$
	print,'Finished interpolating frame ',fr &$
ENDFOR

END
