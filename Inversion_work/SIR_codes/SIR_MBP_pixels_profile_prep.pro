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

restore,'/data/solarstore3/phk/data/27Jul2014/Inversions/SIR/Interpolated_profile000.sav',/ver
restore,'/data/solarstore3/phk/data/27Jul2014/final_tracked_mbps.sav',/ver
restore,'/data/solarstore3/phk/data/27Jul2014/Inversions/SIR/wave_grids.sav',/ver
dims = SIZE(dnew,/dim)
frames = 89
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
	
	;; >>> START PREPPING THE PROFILES TO SAVE OUT
counter = 0. 
bounter = 0.
founter = 0.
np = dims[2]
profile = FLTARR(np,dims[3])
FOR fr = 0, frames-1 DO BEGIN &$
restore,'/data/solarstore3/phk/data/27Jul2014/Inversions/SIR/Interpolated_profile'+STRING(fr,FORMAT='(I3.3)')+'.sav',/ver &$	
FOR x = 0, dims[0]-1 DO BEGIN &$
	FOR y = 0, dims[1] - 1 DO BEGIN &$
	
		; > Root out non-MBP pixels first...
		IF (MBPS_new[x,y,fr] EQ 0) THEN CONTINUE &$
	
		; > Set up 2 IF statements to determine whether the average of the 6301 left lobe is +ve 
		;	or -ve. This then will be used to name the file if it is +VE (p) or -VE (n)
		;	The two separte counters will be used to name the end of the file so they are 
		;	sequential and can be run as a group. x & y is kept in the name so that 
		; 	this info can be gathered later on again to reconstruct the data.
		;	
		; > Might need to consider the scenario that the lobe is zero. Does inclination 
		;	really matter then? I assume it doesn't ....
		
		; > For the Observed Profiles want to create 2 varieties
		;       1: _adjusted_ is those where the interpolated positions are ignored (with val -2)
		;       2: other is those that are interpolated but the values the same as the interpolant to overplot/test
	
		FOR s = 0,np-1 DO profile[s,*] = dnew[x,y,s,*] &$
	
	;>>>	FIRST BATCH - First lobe = -VE & is: gamma (a+b*gam+c*logtau) a,b,c?:180,0,0 !give a,b,c (default 0,1,0), in guess.mtrol
	
		IF (MEAN(profile[21:30,3]) LT 0) THEN BEGIN &$
		;   1:
		inten = DBLARR(np, 4) &$
		inten[*,0] = -2. &$
		inten[idx,0] = profile[idx,0] &$
		inten[83:89,0] = -2 &$ 
		inten[*,3] = profile[*,3] &$
		    fileobsad='Observed_profile_adjusted_x'+STRING(x,FORMAT='(I3.3)')+'_y'+STRING(y,FORMAT='(I3.3)')+'_fr'+STRING(fr,FORMAT='(I3.3)')+'_n'+STRING(counter,FORMAT='(I6.6)')+'.per' &$
		    GET_LUN,unit &$
		    FOR lam= 0, N_ELEMENTS(nwav)-1 DO BEGIN &$
        		tab = DBLARR(6,1) &$
        		IF (lam LT N_ELEMENTS(nwav1)) THEN index = 1. &$
        		IF (lam GE N_ELEMENTS(nwav1)) THEN index = 2. &$
        		IF (lam LT N_ELEMENTS(nwav1)) THEN wavelength = nwav1[lam] &$
        		IF (lam GE N_ELEMENTS(nwav1)) THEN wavelength = nwav2[lam-N_ELEMENTS(nwav1)] &$
        		tab[0,*] = index &$
        		tab[1,*] = wavelength &$
        		tab[2,*] = inten[lam,0] &$
        		tab[3,*] = inten[lam,1] &$
        		tab[4,*] = inten[lam,2] &$
        		tab[5,*] = inten[lam,3] &$
        		OPENU,unit,''+fileobsad+'',width=350,/APPEND &$
        		    printf,unit,tab &$
        		CLOSE,unit &$
		    ENDFOR &$
		    FREE_LUN,unit &$

		; COMMENT THIS OUT AS IT'S NOT REALLY NECESSARY - ADDS TO MANY FILES....
		;   2:
		;inten = DBLARR(np, 4) &$
		;inten[*,0] = profile[*,0] &$
		;; Mask out the telluric in Stokes I
		;inten[83:89,0] = -2 &$ 
		;inten[*,3] = profile[*,3] &$
		;fileobs='Observed_profile_x'+STRING(x,FORMAT='(I3.3)')+'_y'+STRING(y,FORMAT='(I3.3)')+'_fr'+STRING(fr,FORMAT='(I3.3)')+'_n'+STRING(counter,FORMAT='(I6.6)')+'.per' &$
		;GET_LUN,unit &$
		;FOR lam= 0, N_ELEMENTS(nwav)-1 DO BEGIN &$
		;	tab = DBLARR(6,1) &$
		;	IF (lam LT N_ELEMENTS(nwav1)) THEN index = 1. &$
		;	IF (lam GE N_ELEMENTS(nwav1)) THEN index = 2. &$
		;	IF (lam LT N_ELEMENTS(nwav1)) THEN wavelength = nwav1[lam] &$
		;	IF (lam GE N_ELEMENTS(nwav1)) THEN wavelength = nwav2[lam-N_ELEMENTS(nwav1)] &$
		;	tab[0,*] = index &$
		;	tab[1,*] = wavelength &$
		;	tab[2,*] = inten[lam,0] &$
		;	tab[3,*] = inten[lam,1] &$
		;	tab[4,*] = inten[lam,2] &$
		;	tab[5,*] = inten[lam,3] &$
		;	OPENU,3,''+fileobs+'',width=350,/APPEND &$
		;		printf,3,tab &$
		;	CLOSE,3 &$
		;ENDFOR &$
		;FREE_LUN,unit &$

		;   Straylight Profile: Already producedf21 in this case...
		counter = counter + 1. &$
		ENDIF &$
				
	;>>>	SECOND BATCH - First lobe = +VE & is: gamma (a+b*gam+c*logtau) a,b,c?:0,0,0 !give a,b,c (default 0,1,0), in guess.mtrol
	
		IF (MEAN(profile[21:30,3]) GT 0) THEN BEGIN &$
		;   1:
		inten = DBLARR(np, 4) &$
		inten[*,0] = -2. &$
		inten[idx,0] = profile[idx,0] &$
		inten[83:89,0] = -2 &$ 
		inten[*,3] = profile[*,3] &$
		    fileobsad='Observed_profile_adjusted_x'+STRING(x,FORMAT='(I3.3)')+'_y'+STRING(y,FORMAT='(I3.3)')+'_fr'+STRING(fr,FORMAT='(I3.3)')+'_p'+STRING(bounter,FORMAT='(I6.6)')+'.per' &$
		    GET_LUN,unit &$
		    FOR lam= 0, N_ELEMENTS(nwav)-1 DO BEGIN &$
        		tab = DBLARR(6,1) &$
        		IF (lam LT N_ELEMENTS(nwav1)) THEN index = 1. &$
        		IF (lam GE N_ELEMENTS(nwav1)) THEN index = 2. &$
        		IF (lam LT N_ELEMENTS(nwav1)) THEN wavelength = nwav1[lam] &$
        		IF (lam GE N_ELEMENTS(nwav1)) THEN wavelength = nwav2[lam-N_ELEMENTS(nwav1)] &$
        		tab[0,*] = index &$
        		tab[1,*] = wavelength &$
        		tab[2,*] = inten[lam,0] &$
        		tab[3,*] = inten[lam,1] &$
        		tab[4,*] = inten[lam,2] &$
        		tab[5,*] = inten[lam,3] &$
        		OPENU,unit,''+fileobsad+'',width=350,/APPEND &$
        		    printf,unit,tab &$
        		CLOSE,unit &$
		    ENDFOR &$
		    FREE_LUN,unit &$

		;   2:
		;inten = DBLARR(dims[2], 4) &$
		;inten[*,0] = oi[x, y,*] &$
		;; Mask out the telluric in Stokes I
		;inten[83:89,0] = -2 &$ 
		;inten[*,3] = ov[x, y,*] &$
		;fileobs='Observed_profile_x'+STRING(x,FORMAT='(I3.3)')+'_y'+STRING(y,FORMAT='(I3.3)')+'_fr'+STRING(fr,FORMAT='(I3.3)')+'_p'+STRING(bounter,FORMAT='(I6.6)')+'.per' &$
		;GET_LUN,unit &$
		;FOR lam= 0, N_ELEMENTS(nwav)-1 DO BEGIN &$
		;	tab = DBLARR(6,1) &$
		;	IF (lam LT N_ELEMENTS(nwav1)) THEN index = 1. &$
		;	IF (lam GE N_ELEMENTS(nwav1)) THEN index = 2. &$
		;	IF (lam LT N_ELEMENTS(nwav1)) THEN wavelength = nwav1[lam] &$
		;	IF (lam GE N_ELEMENTS(nwav1)) THEN wavelength = nwav2[lam-N_ELEMENTS(nwav1)] &$
		;	tab[0,*] = index &$
		;	tab[1,*] = wavelength &$
		;	tab[2,*] = inten[lam,0] &$
		;	tab[3,*] = inten[lam,1] &$
		;	tab[4,*] = inten[lam,2] &$
		;	tab[5,*] = inten[lam,3] &$
		;	OPENU,3,''+fileobs+'',width=350,/APPEND &$
		;		printf,3,tab &$
		;	CLOSE,3 &$
		;ENDFOR &$
		;FREE_LUN,unit &$

		;   Straylight Profile: Already producedf21 in this case...
		bounter = bounter + 1. &$
		ENDIF &$
	ENDFOR &$
	ENDFOR &$
	;tvim,oi[0:x,0:y-1,0] &$
		medprofile = FLTARR(nw) &$
	
	; >> Straylight Profiles... per frame
	FOR i = 0, nw-1 DO medprofile[i] = MEDIAN(dnew[120:790,170:809,i,0]) &$
        interstray = FLTARR(np) &$		
	interstray = interpol(medprofile, float(iwav), float(nwav)) &$
	filestray='Straylight_profile_f'+STRING(fr,FORMAT='(I3.3)')+'.per' &$
	;GET_LUN,unit &$
	FOR lam= 0, N_ELEMENTS(nwav)-1 DO BEGIN &$
		tab = DBLARR(6,1) &$
		IF (lam LT N_ELEMENTS(nwav1)) THEN index = 1. &$
		IF (lam GE N_ELEMENTS(nwav1)) THEN index = 2. &$
		IF (lam LT N_ELEMENTS(nwav1)) THEN wavelength = nwav1[lam] &$
		IF (lam GE N_ELEMENTS(nwav1)) THEN wavelength = nwav2[lam-N_ELEMENTS(nwav1)] &$
		tab[0,*] = index &$
		tab[1,*] = wavelength &$
		tab[2,*] = interstray[lam] &$
		tab[3,*] = 0. &$
		tab[4,*] = 0. &$
		tab[5,*] = 0. &$
		OPENU,unit,''+filestray+'',width=350,/APPEND &$
		    printf,unit,tab &$
		CLOSE,unit &$
	ENDFOR &$
	FREE_LUN,unit &$

	
	print,'Finished Frame ',fr, ' of ',frames &$
	
ENDFOR

END
