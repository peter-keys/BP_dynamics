;==================================================================================================
; 				SURFACE/BODY WORK	Apr2014 DEFINITIVE RESULTS	
;==================================================================================================
;						DATA
;==================================================================================================
;15Apr2014 	R. Erdelyi Data
restore,'/data/rosa3/oldrosa1/phk/data/wavelet_datasets/Apr2014/Gband_data.sav',/ver
mardata = data
mardata = mardata[50:965,95:965,*]	;Gband
restore,'/data/rosa3/oldrosa1/phk/data/wavelet_datasets/Apr2014/15Apr2014_Gband_osc_power.sav',/ver
marpower = osc_power
marpower = marpower[50:965,95:965,*]	;Gband
restore,'/data/rosa3/oldrosa1/phk/data/wavelet_datasets/Apr2014/wavelet_periods_4s_1533s_gband.sav',/ver
marperiod = save_period

delvar,data,osc_power,save_period

;   -	KEY PARAMETERS FILE
restore,'/data/solarstore3/phk/data/wavelet_datasets/Apr2014/key_parameters_gband_aprdata.sav',/ver
restore,'/data/rosa3/oldrosa1/phk/data/wavelet_datasets/Apr2014/filtered_wavelet_power_gband_Apr2014.sav',/ver
;==================================================================================================
;SLICE_ARRAYS    		DOUBLE    = Array[535, 1452, 36]
;POWER_ARRAYS    		DOUBLE    = Array[535, 136, 36]
;MARPERIOD       		FLOAT     = Array[136]
;MARPOREAREAS    		DOUBLE    = Array[1452]
;MARRMS          		DOUBLE    = Array[1452]
;MARPOREVALS     		DOUBLE    = Array[1452]
;AV_XC           		LONG      =          390
;AV_YC           		LONG      =          454
;MARPOREBINARY   		FLOAT     = Array[535, 1452, 36]
;LEFT_MARTPS     		DOUBLE    = Array[1452, 36]
;RIGHT_MARTPS    		DOUBLE    = Array[1452, 36]
;CENTROID        		DOUBLE    = Array[2, 1452]
;INTEN_LEFTTPS   		DOUBLE    = Array[1452]
;INTEN_RIGHTTPS  		DOUBLE    = Array[1452]
;PERIOD_NORMED_POWER		DOUBLE    = Array[535, 136, 36]
;INTEGRATED_POWER_OVER_SLICE	DOUBLE    = Array[535, 36]
;TOTAL_BINARY    		FLOAT     = Array[191, 171]
;==================================================================================================
;					PROGRAMS
;==================================================================================================
;some programs to be used
.r picks.pro
.r dbjslice.pro
.r /home/phk/idl/wavelet/wave_recon.pro
.r /home/phk/idl/wavelet/wavelet.pro
.r /home/phk/idl/wavelet/wave_signif.pro
.r /home/phk/idl/wavelet/sdbwave.pro

;==================================================================================================
;==================================================================================================
;					PORE PARAMETERS
;==================================================================================================
;	- Establish some key parameters (area & inten related)

x1 = 535
x2 = 860
y1 = 640
y2 = 640
sizing = SIZE(mardata)
xsize = sizing[1]
ysize = sizing[2]
tsize = sizing[3]

;   -	BACKGROUND - choose a 'quiet' region of the data
; here we choose - 75:855,45:240

;   -	VALUE TO CONTOUR PORES WITH
marporevals = DBLARR(tsize)	
FOR i = 0, tsize-1 DO marporevals[i] = (MEDIAN(mardata[75:855,45:240,i]))-(2.5*(STDDEV(mardata[75:855,45:240,i])))

;   -	PORE AREAS:
marporeareas = DBLARR(tsize)	
FOR i = 0, tsize-1 DO marporeareas[i] = TOTAL(mardata[300:490,390:550,i] LE marporevals[i])

;   -	RMS of images:
marrms = DBLARR(tsize)		
FOR i = 0, tsize-1 DO marrms[i] = SQRT(TOTAL(mardata[*,*,i]^2)/N_ELEMENTS(mardata[*,*,i]))

;   -	INTEGRATED INTENSITY OVER ENTIRE PORE:
marporeinten = DBLARR(tsize)
FOR i = 0, tsize-1 DO BEGIN &$ 
	test_bin = mardata[300:490,390:550,i] LE marporevals[i] &$
	test_bin = test_bin * mardata[300:490,390:550,i] &$
	area_num = N_ELEMENTS(WHERE(test_bin GT 0)) &$
	marporeinten[i] = TOTAL(test_bin,/DOUBLE)/area_num &$
ENDFOR
norminten = marporeinten/MAX(marporeinten)
normarea = marporeareas/MAX(marporeareas)

;   -	WORK OUT THE MAIN PERIODS IN THE DATA (AREA & INTEN):
temparea = WAVE_RECON(marporeareas,2.1)
sdbwave,temparea,delt=2.1,/CONE,/FAST,/YLOG
;   -	115.5, 198.9, 357.2 (1.9mins, 3.3mins,5.95mins)

tempinten = WAVE_RECON(marporeinten,2.1)
sdbwave,tempinten,delt=2.1,/CONE,/FAST,/YLOG
;   -	60, 190.8, 349.86  (1min, 3.2min, 5.8mins)

;==================================================================================================
; 				CENTRE OF GRAVITY OF THE PORE
;==================================================================================================
;   - 	Find the CofG of the pore in each frame


tvim,mardata[x1:x2,y1-30:y2+30,0]                        
contour,mardata[x1:x2,y1-30:y2+30,0],level=marporevals[0],/over
x1 = 300
x2 = 490

;   -	TOTAL AREA UNDER THE THRESHOLD
total_area = TOTAL(mardata[x1:x2,y1-50:y2+90,0] LE marporevals[0])  

;   -	The CENTRE OF GRAVITY across the pore
binary = FLTARR(x2-x1+1,141)
binary[*] = 0.
centroid = DBLARR(2,N_ELEMENTS(mardata[0,0,*]))
FOR t = 0, N_ELEMENTS(mardata[0,0,*])-1 DO BEGIN &$
	FOR x = 0, (x2-x1) DO BEGIN &$
		FOR y = 0, 140 DO BEGIN &$
			IF (mardata[x1+x,(y1-50)+y,t] LE marporevals[t]) THEN binary[x,y] = 1. &$
		ENDFOR &$
	ENDFOR &$
	dim = SIZE(binary,/DIMENSIONS) &$
	totalMass = TOTAL(binary) &$
	xcm = TOTAL(TOTAL(binary,2)*INDGEN(dim[0]))/totalMass &$
	ycm = TOTAL(TOTAL(binary,1)*INDGEN(dim[1]))/totalMass &$
	centroid[0,t] = xcm+x1 &$
	centroid[1,t] = ycm+y1-30 &$
	binary[*] = 0. &$
ENDFOR

;   -	Average CENTRE OF GRAVITY
av_xc = MEAN(centroid[0,*])
av_yc = MEAN(centroid[1,*])

;======================================================================================================
;				 MAY NEED TO MAKE SURE DATA IS ALIGNED				
;======================================================================================================
;  -	Check the data to make sure that it is aligned properly
;  -	May have shifts due to the motion of the pore
;  -	Want to correct these
;  -	Look at the centroids and centroid range to see if it is necessary

mult,1,2
plot,centroid[0,*],xst=1,yst=1,title='X drift'
plot,centroid[1,*],xst=1,yst=1,title='Y drift'
mult,1,1

print,MAX(centroid[0,*])-MIN(centroid[0,*])
print,MAX(centroid[1,*])-MIN(centroid[1,*])

;For the Apr2014 Gband data the drift is 27 in x and 18 in y
;Therefore need more aligning:

image = mardata[*,*,0]>0.
x1 = 300
x2 = 490
y1 = 390
y2 = 550
temp = FLTARR(x2-x1+1,y2-y1+1,2)
;temp = FLTARR(911,883,2)
temp[*,*,0] = image[x1:x2,y1:y2]
;temp[*,*,0] = image
ccshifts = FLTARR(2,2)

xshifts = FLTARR(tsize-1)
yshifts = FLTARR(tsize-1)
xshifts[*] = 0.
yshifts[*] = 0.

new_sub = mardata
new_sub[*] = 0.
new_sub[*,*,0] = mardata[*,*,0]>0.
mult,1,2
FOR i = 1,(tsize-1) DO BEGIN &$
    image = mardata[*,*,i]>0. &$
   temp[*,*,1] = image[x1:x2,y1:y2] &$
    ;temp[*,*,1] = image &$
    ccshifts[*] = 0. &$
    FOR j = 0,4 DO ccshifts = ccshifts + TR_GET_DISP(temp,/shift) &$
    ;IF ccshifts[0,1] lt 1. THEN ccshifts[0,1] = 0. &$
    ;IF ccshifts[1,1] lt 1. THEN ccshifts[1,1] = 0. &$
    image = SHIFT(image,ccshifts[0,1],ccshifts[1,1]) &$
    ;image = temp[*,*,1] &$
    xshifts[i-1] = ccshifts[0,1] &$
    yshifts[i-1] = ccshifts[1,1] &$
    ;PUT THE NEXT LINE IN IF YOU WANT CUMULATIVE CO-ALIGNMENT
    ;temp[*,*,0] = image[x1:x2,y1:y2] &$
    ;temp[*,*,0] = image &$
    ;PUT THE NEXT LINE IN IF YOU WANT SEMI-CUMULATIVE CO-ALIGNMENT
    IF i MOD 50 EQ 0 THEN temp[*,*,0] = image[x1:x2,y1:y2] &$
    ;IF i MOD 50 EQ 0 THEN temp[*,*,0] = image &$
    PRINT,'Processing image number ',i,' of ',(tsize-1) &$
    IF (i MOD 10) eq 0 THEN plot,xshifts,thick=2,xtitle='Frame',ytitle='x Pixel shift',xst=1,yst=1,charsize=2 &$
    IF (i MOD 10) eq 0 THEN plot,yshifts,thick=2,xtitle='Frame',ytitle='y Pixel shift',xst=1,yst=1,charsize=2 &$
    new_sub[*,*,i] = image &$
ENDFOR
mult,1,1


;NOW TRY TO DO THE OTHER SHIFT PROCESS TO BRING IT BACK TO THE 
; POINT IT SHOULD BE AT (IE REMOVE RESIDUAL SHIFTS)

;elements = N_ELEMENTS(new_sub[0,0,*])
;temp = new_sub[*,*,0]
;xsize = N_ELEMENTS(temp[*,0])
;ysize = N_ELEMENTS(temp[0,*])
;temp = FLTARR(xsize,ysize,2)
;temp[*,*,0] = new_sub[*,*,0]
;xshifts = FLTARR(elements-1)
;yshifts = FLTARR(elements-1)
;mult,1,2
;FOR i = 1,(elements-1) DO BEGIN &$
;    temp[*,*,1] = new_sub[*,*,i] &$
;    FOR j = 0,4 DO BEGIN &$
;        ccshifts = TR_GET_DISP(temp,/shift) &$
;	xshifts[i-1] = xshifts[i-1] + ccshifts[0,1] &$
;	yshifts[i-1] = yshifts[i-1] + ccshifts[1,1] &$
;    ENDFOR &$	
;    new_sub[*,*,i] = REFORM(temp[*,*,1]) &$
;    temp[*,*,0] = temp[*,*,1] &$
;    IF (i MOD 100) eq 0 THEN print,'Processing image '+STRTRIM(i,2)+' of '+STRTRIM(elements,2) &$
;    IF (i MOD 10) eq 0 THEN plot,xshifts,thick=2,xtitle='Frame',ytitle='x Pixel shift',xst=1,yst=1,charsize=2 &$
;    IF (i MOD 10) eq 0 THEN plot,yshifts,thick=2,xtitle='Frame',ytitle='y Pixel shift',xst=1,yst=1,charsize=2 &$
;ENDFOR
;mult,1,1


;WORK OUT THE CENTROID AND THE BINARY OF THIS TO SEE IF IT IMPROVES IT (RANGE)

;   -	The CENTRE OF GRAVITY across the pore

x1 = 300
x2 = 490
y1 = 390
y2 = 550
binary = FLTARR(xsize,ysize,tsize)
binary[*] = 0.
centroid = DBLARR(2,tsize)
centroid[*] = 0.
FOR t = 0, (tsize-1) DO BEGIN &$
	FOR x = 0, (x2-x1) DO BEGIN &$
		FOR y = 0, (y2-y1) DO BEGIN &$
			IF (new_sub[x+x1,y+y1,t] LE marporevals[t]) THEN binary[x+x1,y+y1,t] = 1. &$
		ENDFOR &$
	ENDFOR &$
	dim = SIZE(binary,/DIMENSIONS) &$
	totalMass = TOTAL(binary[*,*,t]) &$
	xcm = TOTAL(TOTAL(binary[*,*,t],2)*INDGEN(dim[0]))/totalMass &$
	ycm = TOTAL(TOTAL(binary[*,*,t],1)*INDGEN(dim[1]))/totalMass &$
	centroid[0,t] = xcm &$
	centroid[1,t] = ycm &$
	;binary[*] = 0. &$
ENDFOR

plot,centroid[0,*],yst=1

;   -	Still not great alignment (doesn't seem to have changed it at all....) - might look at it again but will continue for now.

mardata = new_sub
delvar,new_sub
;======================================================================================================
;				 RE-ESTABLISH SOME PARAMETERS AFTER ALIGNMENT				
;======================================================================================================
;   -	The CENTRE OF GRAVITY across the pore
;	Other values that may have altered (eg x1,x2,sizing etc.)

x1 = 300
x2 = 490

binary = FLTARR(x2-x1+1,61)
binary[*] = 0.
centroid = DBLARR(2,N_ELEMENTS(mardata[0,0,*]))
FOR t = 0, N_ELEMENTS(mardata[0,0,*])-1 DO BEGIN &$
	FOR x = 0, (x2-x1) DO BEGIN &$
		FOR y = 0, 59 DO BEGIN &$
			IF (mardata[x1+x,(y1-30)+y,t] LE marporevals[t]) THEN binary[x,y] = 1. &$
		ENDFOR &$
	ENDFOR &$
	dim = SIZE(binary,/DIMENSIONS) &$
	totalMass = TOTAL(binary) &$
	xcm = TOTAL(TOTAL(binary,2)*INDGEN(dim[0]))/totalMass &$
	ycm = TOTAL(TOTAL(binary,1)*INDGEN(dim[1]))/totalMass &$
	centroid[0,t] = xcm+x1 &$
	centroid[1,t] = ycm+y1-30 &$
	binary[*] = 0. &$
ENDFOR

;   -	Average CENTRE OF GRAVITY
av_xc = MEAN(centroid[0,*])
av_yc = MEAN(centroid[1,*])

x1 = 140
x2 = 675
y1 = FIX(av_yc)
y2 = FIX(av_yc)
sizing = SIZE(mardata)
xsize = sizing[1]
ysize = sizing[2]
tsize = sizing[3]

;==================================================================================================
;		TIME AVERAGED BINARY OF PORE
;==================================================================================================

;  -	Make a sub field of the pore
x1 = 300
x2 = 490
marsub = mardata[x1:x2,y1-50:y2+120,*]
sizing = SIZE(marsub)
xsub = sizing[1]
ysub = sizing[2]
tsub = sizing[3]

;  -	Get the binaries for each frame of the sub-field
binary = FLTARR(xsub,ysub,tsub)
binary[*] = 0.
FOR t = 0, (tsub-1) DO BEGIN &$
	FOR x = 0, xsub-1 DO BEGIN &$
		FOR y = 0, ysub-1 DO BEGIN &$
			IF (marsub[x,y,t] LE marporevals[t]) THEN binary[x,y,t] = 1. &$
		ENDFOR &$
	ENDFOR &$
ENDFOR


;   -	Add up the binaries to make a complete binary (time-averaged)

total_binary = FLTARR(xsub,ysub)
FOR x = 0, xsub-1 DO BEGIN &$
	FOR y = 0, ysub-1 DO BEGIN &$
		total_binary[x,y] = TOTAL(binary[x,y,*]) &$
	ENDFOR &$
ENDFOR

;	May need to play with the threshold levels before this step

FOR x = 0, xsub-1 DO BEGIN &$
	FOR y = 0, ysub-1 DO BEGIN &$
		IF (total_binary[x,y] LT 900) THEN total_binary[x,y] = 0. &$
		IF (total_binary[x,y] GE 900) THEN total_binary[x,y] = 1. &$
	ENDFOR &$
ENDFOR

;   -	Add up each plot from the bottom of pore to the top (use y vals)
elements = WHERE(total_binary GT 0.)       
xy_elements = xyposition(total_binary,elements)

av_pore_minX = MIN(xy_elements[0,*])	;For boundary locations
av_pore_maxX = MAX(xy_elements[0,*])	;For boundary locations

av_pore_minY = MIN(xy_elements[1,*])	;For vertical totalling 
av_pore_maxY = MAX(xy_elements[1,*])	;For vertical totalling

;   -	Work out the average xc and yc over the pore lifetime (for rotating later)
av_xc1 = ROUND(MEAN(xy_elements[0,*]))
av_yc1 = ROUND(MEAN(xy_elements[1,*]))

;   -	Ammend this for the whole FOV
av_xc = av_xc1 + x1
av_yc = av_yc1 + (y1 - 50)


;==================================================================================================
;					SLICE ARRAYS
;==================================================================================================
;   -	Create arrays at each 5 degree rotation in slice (36 POINTS AS 180/5)
;   -	Done in both intensity and power maps

;  -	Change the x1, x2 and y1/2 values to suit a larger slice (and the y centre)
x1 = 140
x2 = 675
y1 = av_yc	;449
y2 = av_yc

;   -	Slices in Intensity
slice_arrays2 = DBLARR(x2-x1,tsize,36)
rotate_data = FLTARR(xsize,ysize,tsize)
	

;   -	For this data set it has an issue with the av_xc/yc being swayed (due to pore shape) a bit so manually entered the pivot point I want 
;		to ensure that the slice doesn't completely miss the pore.					
FOR i = 0, 35 DO BEGIN &$
	;print,'Working out angle: '+arr2str(i*5,/trim) &$
	FOR j =0, tsize-1 DO BEGIN &$
		;rotate_data[*,*,j] = ROT(mardata[*,*,j],i*5,1,380,445,/PIVOT) &$ ;+VE IS CLOCKWISE!!!
		rotate_data[*,*,j] = ROT(mardata[*,*,j],i*5,1,av_xc,av_yc,/PIVOT) &$ ;+VE IS CLOCKWISE!!!
	ENDFOR &$
	slice =  DBJSLICE(rotate_data,3,x1=x1,x2=x2,y1=y1,y2=y1) &$
	slice = TOTAL(slice,2) &$
	slice_arrays2[*,*,i] = slice &$
	print,'Working out angel: '+arr2str(i*5,/trim) &$
ENDFOR 

;   -	Slices in Power
power_arrays = DBLARR(x2-x1,N_ELEMENTS(marperiod),36)
period_num = N_ELEMENTS(marperiod)

rotate_data = FLTARR(xsize,ysize,N_ELEMENTS(marperiod))
FOR i = 0, 35 DO BEGIN &$
	print,'Working out angle: '+arr2str(i*5,/trim) &$
	FOR j =0, period_num-1 DO BEGIN &$
	rotate_data[*,*,j] = ROT(marpower[*,*,j],i*5,1,380,445,/PIVOT) &$ ;+VE IS CLOCKWISE!!!
	ENDFOR &$
	slice =  DBJSLICE(rotate_data,3,x1=x1,x2=x2,y1=y1,y2=y2) &$
	slice = TOTAL(slice,2) &$
	power_arrays[*,*,i] = slice &$
ENDFOR 


;==================================================================================================
;				PORE BINARY MAPS/TPs
;==================================================================================================
;   -	Work out the pore binary map in slice images
;   -	Work out a binary map of the pore in data images
;   -	Work out the TPs in the slice/data

marporebinary = FLTARR(x2-x1,tsize,36)
marporebinary[*] = 0.

;   -	Remember to change the values below 'y' for the slice to encompass the pore properly

FOR i = 0, 35 DO BEGIN &$
	print,'Working on angle = '+arr2str((i*5),/trim) + ' degrees' &$
	FOR x = 0, tsize-1 DO BEGIN &$
		sigma = MEAN(slice_arrays[*,x,i]) - (1.5*STDDEV(slice_arrays[*,x,i])) &$
		;print,sigma,marporevals[x] &$
		FOR y = 190, 310 DO BEGIN &$		;To ignore points outside pore that fall into the sigma
			IF (slice_arrays[y,x,i] LE sigma) THEN marporebinary[y,x,i] = 1.&$
		ENDFOR &$
	ENDFOR &$
ENDFOR 

;   -	Intensity Turning Points
inten_lefttps = DBLARR(tsize)
inten_righttps = DBLARR(tsize)

FOR t = 0, tsize-1 DO BEGIN &$
	elements = WHERE(mardata[x1:x2,y1,t] LE marporevals[t]) &$       
	xy_elements = xyposition(mardata[x1:x2,y1,t],elements) &$
	inten_lefttps[t] = MIN(xy_elements[0,*]) &$
	inten_righttps[t] = MAX(xy_elements[0,*]) &$
	print,'Frame num processing = '+arr2str(t,/trim) &$
ENDFOR

;   -	Binary Turning Points
left_martps = DBLARR(tsize,36)
right_martps = DBLARR(tsize,36)

FOR i = 0, 35 DO BEGIN &$
	print,'Currently working out angle = ' +arr2str(i*5,/trim) + ' degress' &$
	FOR t = 0, tsize-1 DO BEGIN &$
		elements = WHERE(marporebinary[*,t,i] EQ 1.) &$       
		xy_elements = xyposition(marporebinary[*,t,i],elements) &$
		left_martps[t,i] = MIN(xy_elements[0,*]) &$
		right_martps[t,i] = MAX(xy_elements[0,*]) &$
	ENDFOR &$
ENDFOR


;==================================================================================================
;		INTEGRATED AND NORMALISED POWER OVER SLICE
;==================================================================================================

;THEREFORE WILL WORK OUT THE INTEGRATED POWER FOR THE SLICE OVER ALL PERIODS TO SEE PEAK LOCATION
;Want to mormalise the power_arrays by P^2 first before you attempt it
; period_normed power = the power at each period normalised by period squared
period_normed_power = DBLARR(x2-x1,period_num,36)

FOR i = 0, 35 DO BEGIN &$
	FOR j = 0, period_num-1 DO BEGIN &$
		period_normed_power[*,j,i] = power_arrays[*,j,i] / (marperiod[j]^2) &$
	ENDFOR &$
ENDFOR


integrated_power_over_slice = DBLARR(x2-x1,36)

FOR j = 0, 35 DO BEGIN &$
	FOR i = 0, (x2-x1-1) DO BEGIN &$
		integrated_power_over_slice[i,j] = TOTAL(period_normed_power[i,*,j]) &$
	ENDFOR &$
ENDFOR

;==================================================================================================
; UN INTERLUDE OF STUFF TO DATE:
;==================================================================================================
;KEY VALUES HAVE BEEN SAVED IN THIS FILE:
;SAVE,FILENAME='/data/rosa3/oldrosa1/phk/data/wavelet_datasets/Apr2014/key_parameters_gband_aprdata.sav',slice_arrays,power_arrays,marperiod,marporeareas,marporeinten,marrms,marporevals,av_xc,av_yc,marporebinary,left_martps,right_martps,centroid,inten_lefttps,inten_righttps,period_normed_power,integrated_power_over_slice,total_binary
;KEY VALUES SAVED ARE:
;	+ slice_arrays = array of 5deg rotations of slices about the pore
;	+ power_arrays = corresponding array of 5deg power plot slices
;	+ marperiod = periods from the power info
;	+ marporeareas = Area of pore in each image
;	+ marrms = rms at each frame in mardata
;	+ marporevals = the mean-3sigma level for each frame
;	+ av_xc = average centre of the pore (x)
;	+ av_yc = average centre of the pore (y)
;	+ marporebinary = binary map of pore boundaries of slice for each angle
;	+ left_martps = left edges of pore in each frame at each angle (with the slice arrays)
;	+ right_martps = right edges of pore in each frame at each angle (with the slice arrays)
;
;	+ centroid = (x,y) of cofg of the pore & for phase relations
;	+ inten_left/righttps = the left and right TPs based on the inten images (opposed to slice)
;	+ period_normed_power = the power slices normalised by the period
;	+ integrated_power_over_slice = The normalised power integrated over all periods (various rotations)
;	+ total_binary = The time-averaged binary of the pore. Threshold used to get the best binary 
;==================================================================================================
;   -	Little RE-CAP:

;   -	WORK OUT THE MAIN PERIODS IN THE DATA (AREA & INTEN):
temparea = WAVE_RECON(marporeareas,2.1)
sdbwave,temparea,delt=2.1,/CONE,/FAST,/YLOG
;   -	115.5, 198.9, 357.2 (1.9mins, 3.3mins,5.95mins)

tempinten = WAVE_RECON(marporeinten,2.1)
sdbwave,tempinten,delt=2.1,/CONE,/FAST,/YLOG
;   -	60, 190.8, 349.86  (1min, 3.2min, 5.8mins)
;==================================================================================================
;				 PHASE AND COHERENCE IN DATA				
;======================================================================================================
;  -	Between Intensity and Area

;   -	DEFINE YOUR LIGHTCURVES
lc_1 = marporeareas
lc_2 = marporeinten
    
;   -	REMOVE LOW-FREQUENCY POWER WHICH ARTIFICIALLY SWAMPS SIGNAL
;	   -	HERE 'delt' IS THE CADENCE
delt = 2.112
lc_1 = WAVE_RECON(lc_1, delt)
lc_2 = WAVE_RECON(lc_2, delt)

;   -	TAKE WAVELET TRANSFORMS OF THE LCs
;	   -	HERE I WOULD SET 'deejay' TO 16. OR 32. DEPENDING ON THE FREQUENCY RESOLUTION YOU WANT
deejay = 16.
wave_arr_1 = WAVELET(lc_1 , delt, PERIOD=period, SCALE=scale, COI=coi, MOTHER='morlet', DJ=1./deejay) 
wave_arr_2 = WAVELET(lc_2 , delt, PERIOD=period, SCALE=scale, COI=coi, MOTHER='morlet', DJ=1./deejay) 
    
;   -	ESTABLISH THE POWER OF LC_1 AND LC_2
P_1 = ABS(wave_arr_1)^2.
P_2 = ABS(wave_arr_2)^2.

;   -	NORMALISE THE POWER SPECTRA
P_1 = P_1 / MAX(P_1)
P_2 = P_2 / MAX(P_2)

;   -	CREATE THE CROSS-POWER SPECTRUM
cps = wave_arr_1 * CONJ( wave_arr_2 ) 
     
;   -	ESTABLISH THE PHASE DIFFERENCES BETWEEN THE POWER SPECTRA
phd = REFORM(ATAN(IMAGINARY(cps),REAL_PART(cps)))*180./!pi 

;   -	DETERMINE THE COHERENCE OF THE OSCILLATIONS FROM THE CROSS-POWER SPECTRA
cps = cps / MAX(cps)
coherence = ABS(cps)^2.

;   -	PHASE DIFFERENCE DIAGRAM:

@velocity_colour_table.bat
tvim, phd, range=[-90,90], /sc
contour, coherence, levels=[0.4, 0.6, 0.8], /over

;where the y axis will be the period elements (i.e. use the 'print, period' to display actual values) 
;and the x axis will be time (i.e. the same as FINDGEN(elements)*delt). I sometimes make a mask 
;of the coherence array (i.e. masking out all values below a certain threshold) so I can multiply this 
;binary-map by the phase difference array to give a better indication of where important results lie. 

;   -	MAKE A BINARY MAP OF THE COHERENCE AT A SPECIFIC LEVEL

level = 0.4
coherence_bin = coherence
coherence_bin[*] = 0.

FOR x = 0, (N_ELEMENTS(coherence[*,0])-1) DO BEGIN &$
	FOR y = 0, (N_ELEMENTS(coherence[0,*])-1) DO BEGIN &$
		IF (coherence[x,y] GE level) THEN coherence_bin[x,y] = 1. &$
	ENDFOR &$
ENDFOR
tvim, phd*coherence_bin, /sc,title='Phase Difference for Mar2013 Data',xtitle='Frame Num',ytitle='Period Array Val'            


;======================================================================================================
;			PHASE WITH RESPECT TO THE CENTRE OF THE PORE 
;======================================================================================================
;   -	Work out the phase w.r.t the centre of the pore

;	ADD IN THE POINT OF HIGHEST POWER VERSION!!!! (Run over the weekend cause it takes a while)
tvim,mardata[x1:x2,y1-50:y1+90,0]
contour,mardata[x1:x2,y1-50:y2+90,0],/OVER,level=marporevals[0] 

;   -	Lightcurve of the Centre-point
deejay = 16.
centre_lc = DBLARR(tsize)
FOR t = 0, tsize-1 DO centre_lc[t] = mardata[(ROUND(av_xc)),(ROUND(av_yc)),t]

phase_array = FLTARR(x2-x1,141,(N_ELEMENTS(period)))
phase_array[*] = 0.

FOR x = 0, (x2-x1-1) DO BEGIN &$
	FOR y = 0, 140 DO BEGIN &$
		point = DBLARR(tsize) &$
		IF (mardata[x1+x,(y1-50)+y,0] LE marporevals[0]) THEN BEGIN &$
			FOR t = 0, tsize-1 DO point[t] =  mardata[x1+x,(y1-50)+y,t] &$

			lc_1 = centre_lc &$
			lc_2 = point &$

			; REMOVE LOW-FREQUENCY POWER WHICH ARTIFICIALLY SWAMPS SIGNAL
			; HERE 'delt' IS THE CADENCE
			delt = 2.112 &$
			lc_1 = WAVE_RECON(lc_1, delt)  &$
			lc_2 = WAVE_RECON(lc_2, delt) &$

			; TAKE WAVELET TRANSFORMS OF THE LCs
			; HERE I WOULD SET 'deejay' TO 16. OR 32. DEPENDING ON THE FREQUENCY RESOLUTION YOU WANT
			deejay = 16. &$
			wave_arr_1 = WAVELET(lc_1 , delt, PERIOD=period, SCALE=scale, COI=coi, MOTHER='morlet', DJ=1./deejay)  &$
			wave_arr_2 = WAVELET(lc_2 , delt, PERIOD=period, SCALE=scale, COI=coi, MOTHER='morlet', DJ=1./deejay)  &$

			; ESTABLISH THE POWER OF LC_1 AND LC_2
			P_1 = ABS(wave_arr_1)^2. &$
			P_2 = ABS(wave_arr_2)^2. &$

			; NORMALISE THE POWER SPECTRA
			P_1 = P_1 / MAX(P_1) &$
			P_2 = P_2 / MAX(P_2) &$

			; CREATE THE CROSS-POWER SPECTRUM
			cps = wave_arr_1 * CONJ(wave_arr_2)  &$

			; ESTABLISH THE PHASE DIFFERENCES BETWEEN THE POWER SPECTRA
			phd1 = REFORM(ATAN(IMAGINARY(cps),REAL_PART(cps)))*180./!pi  &$

			; DETERMINE THE COHERENCE OF THE OSCILLATIONS FROM THE CROSS-POWER SPECTRA
			cps = cps / MAX(cps) &$
			coherence = ABS(cps)^2. &$
			
			;   -	FIND THE MEAN PHASE FOR EACH PIXEL AT EACH PERIOD
			FOR q = 0, (N_ELEMENTS(period)-1) DO phase_array[x,y,q] = MEAN(phd1[*,q]) &$
		ENDIF &$

	ENDFOR &$
ENDFOR

; NOT 100% CONVINCED THIS IS RIGHT BUT I'LL HAVE A LOOK AT IT AGAIN, AFTER WORKING OUT SOME OTHER THINGS

;==================================================================================================
; SAVE OUT THE PHASE MILARKY
;==================================================================================================
;KEY VALUES HAVE BEEN SAVED IN THIS FILE:
;SAVE,FILENAME='/data/rosa3/oldrosa1/phk/data/wavelet_datasets/Apr2014/phase_parameters_gband_aprdata.sav',phd,period,coherence_bin,phase_array,imfs_areas,imfs_inten
;	+ phd = Phase of the Area/Intensity Lightcurves
;	+ coherence_bin = coherence of the area/intensity lightcurves in 2D
;	+ phase_array = the average phase at each location in the pore
;==================================================================================================

;======================================================================================================
;					EMD OF AREA/INTEN SIGNALS
;======================================================================================================

; I'VE PUT THIS OFF FOR LIKE FOR EVS SO OFF I GO TO DO EMD ON THE INTENSITY/AREA SIGNALS
; Made some modifications off the other code to try and get it to out put the images to the terminal
; This will make it easier to actually analyse over the PS file version

.r /home/phk/idl/wavelet/nabil_code/emd_test.pro

signal = marporeareas
sda = 0.2
print,sda

;EMD,signal,sda,DT=2.112
imfs_areas = EMD_FUNCTION(signal,sda,DT=2.112)

; Have EMD working alright now and outputting the plots to the terminal.
; The function version out puts the individual IMFs with EMD so that I can look at them

c_x = DBLARR(tsize)		; The individual IMF at Amplitude x (e.g. 2, 3,.... ignore 1 as it is noise)
c_x[*] = imfs_inten[11,*]
c_x = WAVE_RECON(c_x,2.112)

sdbwave,c_x,delt=2.112,/FAST,/YLOG,/CONE
 ;______________________________________________________________________|_______________________________________________________________________|
 ;RESULTS		(G-BAND)					|									|
 ;	marporeinten:							|	marporeareas:							|
 ;		- c2 = 	17.33s						|		- c2 = 	12.9s (22.28s)					|
 ;		- c3 = 	31.79s						|		- c3 = 	25.26s 						|
 ;		- c4 = 	55.90s						|		- c4 = 	31.79s (47.29 blip & 118.65s v. small) 		|
 ;		- c5 = 	98.3s  (67.47 less dominant)			|		- c5 = 	45.36s (83.16s blip & 169.3s v. small) 		|
 ;		- c6 =  129s  	(200s less dominant)			|		- c6 = 	98.3s 						|
 ;		- c7 =  165.8s 						|		- c7 = 	134.51s (217.55s dom but out & 415.91s) 	|
 ;		- c8 =  303.96s 					|		- c8 = 	310.38s 					|
 ;		- c9 =  359.29s	(569.08s less dominant)			|		- c9 = 	390.62s  (557.31s slightly out)  		|
 ;		- c10 = 523.4s (starting to move out) 			|		- c10 = 762.57s (just out) 				|
 ;		- c11 = 658.8s  (mostly out) 				|		- c11 = 605.91s (out)					|
 ;		- c12 = 846.58s (well out)				|		- c12 = 658.76s (well out)				|
 ;______________________________________________________________________|_______________________________________________________________________|

 ; Want to look at the same periods as the Wavelet (i.e. 
 ; Major periods in the Data (with wavelet) are:
 ;	-  Area peaks: 115.5, 198.9, 357.2 (1.9mins, 3.3mins,5.95mins)
 ;	-  Intensity peaks: 60, 190.8, 349.86  (1min, 3.2min, 5.8mins)
 
 ; Therefore look at:
 ;			98s, 130s, 190s, 350s	(1.6mins, 2.2mins, 3.2mins, 5.8mins)		

 ; For IMFs look at:
 ; 			c5/c6, c6/c7, c8
			 
c8_inten = IMFS_inten[7,*]/MAX(ABS(IMFS_inten[7,*]))
c8_areas = IMFS_areas[7,*]/MAX(ABS(IMFS_areas[7,*]))

c7_inten = IMFS_inten[5,*]/MAX(ABS(IMFS_inten[5,*]))
c7_areas = IMFS_areas[6,*]/MAX(ABS(IMFS_areas[6,*]))

c5_inten = IMFS_inten[4,*]/MAX(ABS(IMFS_inten[4,*]))
c5_areas = IMFS_areas[5,*]/MAX(ABS(IMFS_areas[5,*]))


;SAVE,FILENAME='/data/solarstore3/phk/data/wavelet_datasets/Apr2014/15Apr2014_IMFs.sav',imfs_areas,imfs_inten,c5_inten,c5_areas,c7_inten,c7_areas,c8_inten,c8_areas
;print,MEANABSDEV(c7_inten),MEANABSDEV(c8_inten),MEANABSDEV(c10_inten)
print,(MEANABSDEV(c7_inten)/MEAN(marporeinten))*100.,(MEANABSDEV(c8_inten)/MEAN(marporeinten))*100.,(MEANABSDEV(c5_inten)/MEAN(marporeinten))*100.

;   -  DELTA I AND DELTA A ARE THE MEANABSDEV OF THE IMFS
;   -  TO GET THE ERROR IN THESE MEASUREMENTS LOOK AT MAX/MIN
;   -  NB USE A cX_inten/area that is not normalised (as it is above)
deltaI = DBLARR(3)
deltaI[0] = (MEANABSDEV(c5_inten)/MEAN(marporeinten))*100.
deltaI[1] = (MEANABSDEV(c7_inten)/MEAN(marporeinten))*100.
deltaI[2] = (MEANABSDEV(c8_inten)/MEAN(marporeinten))*100.

	; In % ....
deltaA = DBLARR(3)
deltaA[0] = (MEANABSDEV(c5_areas)/MEAN(marporeareas))*100.
deltaA[1] = (MEANABSDEV(c7_areas)/MEAN(marporeareas))*100.
deltaA[2] = (MEANABSDEV(c8_areas)/MEAN(marporeareas))*100.

;   - ERRORS IN THE DELTA COMPONENETS SHOULD BE +/- MAX -> MIN value....

err_dI = (MAX(deltaI) - MIN(deltaI)) /2.
err_dA = (MAX(deltaA) - MIN(deltaA)) /2.

;   - ERRORS IN THE AREA / INTENSITY RAW SIGNALS....

;err_inten = MEANABSDEV(marporeinten)
;err_area = MEANABSDEV(marporeareas)	; IN pixels

err_inten = (MEAN(deltaI)/100.) * MEAN(marporeinten)
err_area = (MEAN(deltaA)/100.) * (SQRT(((MEAN(marporeareas))*(((0.069^2)*(725/1000.)))/!PI))*2)

av_diameter = SQRT(((MEAN(marporeareas))*(((0.069^2)*(725/1000.)))/!PI))*2	; IN Mm
av_inten = MEAN(marporeinten)
;print,SQRT(((MEANABSDEV(marporeareas))*(((0.069^2)*(725/1000.)))/!PI))*2	; IN Mm

;   - PRINT OUT THE ERRORS NOW FOR CONSISTENCY....
print,' ' &$
Print,'Average Diameter: '+arr2str(av_diameter,/trim)+' Mm '&$
Print,'Error in Area: '+arr2str(err_area,/trim)+' Mm (diameter)' &$
Print,'Delta A: '+arr2str(MEAN(deltaA),/trim)+' %' &$
Print,'Error in delta Area: '+arr2str(err_dA,/trim)+' %' &$
Print,'Average Intensity: '+arr2str(av_inten,/trim) &$
Print,'Error in Inten: '+arr2str(err_inten,/trim) &$
Print,'Delta I: '+arr2str(MEAN(deltaI),/trim)+' %' &$
Print,'Error in delta Inten: '+arr2str(err_dI,/trim)+' %' &$
Print,' ' 


;Take Intensity err as MEANABSDEV of inten here...



 ;   -	Overplot the associated IMFs for each period
time_frame = FINDGEN(tsize)*2.112
window,0,xsize=1200,ysize=800			

;C5
cgplot,time_frame,c5_inten,xst=1,yst=7,title = 'c!d'+STRCOMPRESS(strtrim(5),/REMOVE)+'!n ',xtitle='Time (s)',charsize=2.,charthick=1.8,background='white',thick=2
cgoplot,time_frame,c5_areas,color='red',thick=2
cgAXIS,yaxis=0,yrange=[MIN(imfs_inten[4,*]),MAX(imfs_inten[4,*])],/SAVE,Ytitle=' Intensity ',charsize=2.,charthick=1.8,yticklen=0.0055
cgAXIS,yaxis=1,yrange=[MIN(imfs_areas[5,*])*(100^2),MAX(imfs_areas[5,*])*(100^2)],/SAVE,Ytitle=' Area ',charsize=2.,charthick=1.8,yticklen=0.0055
X2JPEG,'c5_15Apr2014_98s.jpg'

;C7
cgplot,time_frame,c7_inten,xst=1,yst=7,title = 'c!d'+STRCOMPRESS(strtrim(7),/REMOVE)+'!n ',xtitle='Time (s)',charsize=2.,charthick=1.8,background='white'
cgoplot,time_frame,c7_areas,color='red'
cgAXIS,yaxis=0,yrange=[MIN(imfs_inten[5,*]),MAX(imfs_inten[5,*])],/SAVE,Ytitle=' Intensity ',charsize=2.,charthick=1.8,yticklen=0.0055
cgAXIS,yaxis=1,yrange=[MIN(imfs_areas[6,*])*(100^2),MAX(imfs_areas[6,*])*(100^2)],/SAVE,Ytitle=' Area ',charsize=2.,charthick=1.8,yticklen=0.0055
X2JPEG,'c7_15Apr2014_130s.jpg'

;C8
cgplot,time_frame,c8_inten,xst=1,yst=7,title = 'c!d'+STRCOMPRESS(strtrim(8),/REMOVE)+'!n ',xtitle='Time (s)',charsize=2.,charthick=1.8,background='white'
cgoplot,time_frame,c8_areas,color='red'
cgAXIS,yaxis=0,yrange=[MIN(imfs_inten[7,*]),MAX(imfs_inten[7,*])],/SAVE,Ytitle=' Intensity ',charsize=2.,charthick=1.8,yticklen=0.0055
cgAXIS,yaxis=1,yrange=[MIN(imfs_areas[7,*])*(100^2),MAX(imfs_areas[7,*])*(100^2)],/SAVE,Ytitle=' Area ',charsize=2.,charthick=1.8,yticklen=0.0055
X2JPEG,'c8_15Apr2014_305s.jpg'

;======================================================================================================
;			FFT FILTERED POWER
;====================================================================================================== 100, 190, 350 
;  -	FFT Power plots across the pore
;  -	Have been filtered at the desired frequencies
;  -	Takes total power at each point along the slice
;  -	This may not be perfect yet
;  -	Freqs that we want to look at are:
;	  +	350s = 2.86mHz
;	  +	190s = 5.26mHz
;	  +	100s = 10mHz
;  -	The widths (from Rick) appear to be roughly double the centre
;  -	Therefore I'll use the following widths:
;	  +	700s = 1.43mHz
;	  +	380s = 2.63mHz
;	  +	200s = 5mHz
;  -	Arrays are named as:
;	  +	filtered_data3mhz 	DOUBLE    = Array[916, 871, 1452]
;	  +	filtered_data5mhz	DOUBLE    = Array[916, 871, 1452]
;	  +	filtered_data10mhz	DOUBLE    = Array[916, 871, 1452]
;	  +	filtered_power3mhz 	DOUBLE    = Array[916, 871, 1452]
;	  +	filtered_power5mhz	DOUBLE    = Array[916, 871, 1452]
;	  +	filtered_power10mhz	DOUBLE    = Array[916, 871, 1452]

; -	Convert to freq in mHz
print,(1/(90.))*1000.

sizing = SIZE(mardata)
xsize = sizing[1]
ysize = sizing[2]
tsize = sizing[3]

freq = 10./1000.
;width = 5./1000.
width = (freq + freq/10.) - (freq - freq/10.) 
dt = 2.112
sz=tsize			;Size of data set
n21=sz/2+1				;Half width
f=(findgen(sz))/sz/dt			;Freq. Array
f[n21]=(n21-sz+findgen(n21-2))/sz/dt	;Freq (-ves)
cent=freq*dt*sz 
width2=width*dt*sz
psf=psf_gaussian(npixel=sz,fwhm=width2,centroid=cent,ndimen=1)
;psf2=psf_gaussian(npixel=sz,fwhm=width2,centroid=zsize-cent,ndimen=1)
psf2=psf_gaussian(npixel=sz,fwhm=width2,centroid=sz-cent,ndimen=1)

filtered_power10mHz = DBLARR(xsize,ysize,(tsize))
filtered_data10mHz = DBLARR(xsize,ysize,(tsize))
FOR x = 0, xsize-1 DO BEGIN &$
	FOR y = 0, ysize-1 DO BEGIN &$
		signal = DBLARR(tsize) &$
		signal[*] = 0. &$
		signal[*] = mardata[x,y,*] &$	;Changed for the properly aligned data
		freqspect = FFT(signal,-1) &$
		filt_freq = freqspect*psf + freqspect*psf2 &$
		filtered_data10mHz[x,y,*] = FFT(filt_freq,1) &$
		power = (ABS(filt_freq))^2 &$
		;filt_power = power[0:(tsize/2)-1] &$
		power = power &$;/tsize &$
		;filtered_power[x,y,*] = FFT(power,1) &$
		filtered_power10mHz[x,y,*] = power &$
	ENDFOR &$
ENDFOR
;SAVE,FILENAME='/data/rosa3/oldrosa1/phk/data/wavelet_datasets/30Sept2012/fourier_filtered_parameters_gband_aprdata.sav',filtered_data3mhz,filtered_data5mhz,filtered_data11mhz,filtered_power3mhz,filtered_power5mhz,filtered_power11mhz,$
;							filt_slice3mHz,filt_slice5mHz,filt_slice11mHz,filt_pslice3mHz,filt_pslice5mHz,filt_pslice11mHz,power_2d_3mhz,$
;							power_2d_5mhz,power_2d_11mhz

;  -	Save them individually:
SAVE,FILENAME='/data/rosa3/oldrosa1/phk/data/wavelet_datasets/Apr2014/fourier_filtered_aprdata_3mHz_narrow.sav',filtered_data3mhz,filtered_power3mhz,filt_slice3mHz,filt_pslice3mHz

;----------------------------------------------
;One of the Figs looks better with filtering
;the slice as opposed to the data. Therefore:
;----------------------------------------------
;  -	Figure looks better this way
;  -	Take the same freq and widths
;  -	Remember to change the freqs above
;  -	The eventual arrays will be:
;	  +	filt_slice3mHz
;	  +	filt_slice5mHz
;	  +	filt_slice10mHz
;	  +	filt_pslice3mHz
;	  +	filt_pslice5mHz
;	  +	filt_pslice10mHz


sizing = SIZE(slice_arrays)
xdim = SIZING[1]			;xsize of slice
ydim = SIZING[2]			;ysize of slice
rotdim = SIZING[3]		;rotation of slice

;  -	FIND FOR EACH ROTATION OF THE SLICE (AND EACH FREQ.) :
filt_slice10mHz = FLTARR(xdim,ydim,rotdim)
filt_pslice10mhz = FLTARR(xdim,ydim,rotdim)
FOR r = 0, rotdim-1 DO BEGIN &$
  FOR x = 0, xdim-1 DO BEGIN &$
    signal = DBLARR(ydim) &$
    signal[*] = 0. &$
    signal[*] = slice_arrays[x,*,r] &$	
    freqspect = FFT(signal,-1) &$
    filt_freq = freqspect*psf + freqspect*psf2 &$
    filt_slice10mHz[x,*,r] = FFT(filt_freq,1) &$		;<--- change for freq.
    power = (ABS(filt_freq))^2 &$
    ;power = power/ydim &$
    filt_pslice10mHz[x,*,r] = power &$			;<--- change for freq.
  ENDFOR &$
ENDFOR

;SAVE,FILENAME='/data/rosa3/oldrosa1/phk/data/wavelet_datasets/Apr2014/fourier_filtered_gband_aprdata_narrow.sav',filtered_data3mhz,filtered_data5mhz,filtered_data10mhz,filtered_power3mhz,filtered_power5mhz,filtered_power10mhz,$
;							filt_slice3mHz,filt_slice5mHz,filt_slice10mHz,filt_pslice3mHz,filt_pslice5mHz,filt_pslice10mHz

;  -	Save them individually:
SAVE,FILENAME='/data/rosa3/oldrosa1/phk/data/wavelet_datasets/Apr2014/fourier_filtered_aprdata_3mHz_narrow.sav',filtered_data3mhz,filtered_power3mhz,filt_slice3mHz,filt_pslice3mHz
SAVE,FILENAME='/data/rosa3/oldrosa1/phk/data/wavelet_datasets/Apr2014/fourier_filtered_aprdata_5mHz_narrow.sav',filtered_data5mhz,filtered_power5mhz,filt_slice5mHz,filt_pslice5mHz
SAVE,FILENAME='/data/rosa3/oldrosa1/phk/data/wavelet_datasets/Apr2014/fourier_filtered_aprdata_10mHz_narrow.sav',filtered_data10mhz,filtered_power10mhz,filt_slice10mHz,filt_pslice10mHz

;----------------------------------------------------------------------
;THE FOURIER POWER IN 2D ACROSS THE PORE (NORMALISED TO THE AVERAGE):
;----------------------------------------------------------------------
;  -	Produce 2D Power Map of the filtered power
;  -	This allows a plot of power across pore
;  -	Gives image to determine if it is surface/body mode
;  -	Rotate the image to get it at multiple angles
;  -	Eventual data sets are:
;	  +	power_2d_3mhz
;	  +	power_2d_5mhz
;	  +	power_2d_11mhz

av = MEAN(filtered_power10mHz)
power_2d_10mhz = FLTARR(xsize,ysize,rotdim)
rotated_images = FLTARR(xsize,ysize,tsize)
FOR r = 0, rotdim-1 DO BEGIN &$
	FOR t = 0, tsize-1 DO BEGIN &$
		rotated_images[*,*,t] = ROT(filtered_power10mhz[*,*,t],r*5,1,av_xc,av_yc,/PIVOT) &$
	ENDFOR &$
	FOR x = 0, xsize-1 DO BEGIN &$
		FOR y = 0, ysize-1 DO BEGIN &$	
			power_2d_10mhz[x,y,r] = (TOTAL(rotated_images[x,y,*]))/av &$
		ENDFOR &$
	ENDFOR &$
	print,'Working out angle: '+arr2str(r*5,/trim) &$
ENDFOR
plot,power_2d_10mHz[x1:x2,y1,0]

;MAYBE PRODUCE A VIDEO OF THE PLOTS TO SEE WHICH ANGLE IS THE BEST FOR EACH

;==================================================================================================
;SAVE THESE FUCKERS OUT
;SAVE,FILENAME='/data/rosa3/oldrosa1/phk/data/wavelet_datasets/Apr2014/fourier_filtered_parameters_gband_aprdata_new.sav',filtered_data3mhz,filtered_data5mhz,filtered_data10mhz,filtered_power3mhz,filtered_power5mhz,filtered_power10mhz,$
;							filt_slice3mHz,filt_slice5mHz,filt_slice10mHz,filt_pslice3mHz,filt_pslice5mHz,filt_pslice10mHz,power_2d_3mhz,$
;							power_2d_5mhz,power_2d_10mhz
;	+ filtered_data3mhz	=  dataset fourier filtered by 2.86mHz	
;	+ filtered_data5mhz	=  dataset fourier filtered by 5.26mHz
;	+ filtered_data10mhz	=  dataset fourier filtered by 10mHz
;	+ filtered_power3mhz	=  filtered power of whole FoV
;	+ filtered_power5mhz	=  filtered power of whole FoV
;	+ filtered_power10mhz	=  filtered power of whole FoV
;	+ filt_slice3mHz	=  slice_arrays fourier filtered by 2.86mHz
;	+ filt_slice5mHz	=  slice_arrays fourier filtered by 5.26mHz
;	+ filt_slice10mHz	=  slice_arrays fourier filtered by 10mHz
;	+ filt_pslice3mHz	=  filtered power of slice_arrays
;	+ filt_pslice5mHz	=  filtered power of slice_arrays
;	+ filt_pslice10mHz	=  filtered power of slice_arrays
;	+ power_2d_3mhz		=  2D power map of total/averaged power (at multiple angles)
;	+ power_2d_5mhz		=  2D power map of total/averaged power (at multiple angles)
;	+ power_2d_10mhz	=  2D power map of total/averaged power (at multiple angles)
;==================================================================================================


;======================================================================================================
;			WAVELET POWER
;======================================================================================================
;  -	Wavelet power of the pore
;  -	Plot this at the suitable periods to see the power
;  -	Do this on the filtered data may haps?
;	  +	350s = 2.86mHz
;	  +	190s = 5.26mHz
;	  +	100s = 10mHz
;  -	Widths:
;	  +	700s = 1.43mHz
;	  +	380s = 2.63mHz
;	  +	200s = 5mHz

osc_power = TOTAL_WAVELET_POWER(filtered_data10mhz,2.112,4.224,1533.31,16.)
osc_power5mHz = osc_power

;  -	Add up power across selected periods to get a 2d power map

wavelet_power2d_5mHz = FLTARR(xsize,ysize,rotdim)
cent_f = 5.26/1000.
;width_f = 2.63/1000.
width_f = (cent_f + cent_f/10.) - (cent_f - cent_f/10.) ;1.79/1000.
lower_range = cent_f - (width_f * 0.5)
upper_range = cent_f + (width_f * 0.5)
period_range = [1/upper_range,1/lower_range]

wavelet_power2d_5mHz[*] = 0.
rotated_osc = osc_power5mHz
rotated_osc[*] = 0.
FOR r = 0, rotdim -1 DO BEGIN &$
	FOR t = 0, (N_ELEMENTS(save_period)-1) DO BEGIN &$
		rotated_osc[*,*,t] = ((ROT(osc_power5mHz[*,*,t],r*5,1,av_xc,av_yc,/PIVOT))/(save_period[t]^2)) &$
	ENDFOR &$
	FOR t = 0, (N_ELEMENTS(save_period)-1) DO BEGIN &$
		IF save_period[t] GE period_range[0] THEN BEGIN &$
		IF save_period[t] LE period_range[1] THEN BEGIN &$
			FOR x = 0, xsize-1 DO BEGIN &$
				FOR y = 0, ysize-1 DO BEGIN &$
					wavelet_power2d_5mHz[x,y,r] = wavelet_power2d_5mHz[x,y,r] + rotated_osc[x,y,t] &$
				ENDFOR &$
			ENDFOR &$
		ENDIF &$
		ENDIF &$
	ENDFOR &$
	print,'Finished rotation '+arr2str(r*5,/trim) &$
ENDFOR

;  -	Save them out individually again:

SAVE,FILENAME='/data/rosa3/oldrosa1/phk/data/wavelet_datasets/Apr2014/filtered_wavelet_mardata_3mHz_narrow.sav',osc_power3mHz,wavelet_power2d_3mHz
SAVE,FILENAME='/data/rosa3/oldrosa1/phk/data/wavelet_datasets/Apr2014/filtered_wavelet_mardata_5mHz_narrow.sav',osc_power5mHz,wavelet_power2d_5mHz
SAVE,FILENAME='/data/rosa3/oldrosa1/phk/data/wavelet_datasets/Apr2014/filtered_wavelet_mardata_10mHz_narrow.sav',osc_power10mHz,wavelet_power2d_10mHz


;==================================================================================================
;SAVE,FILENAME='/data/rosa3/oldrosa1/phk/data/wavelet_datasets/Apr2014/filtered_wavelet_power_gband_Apr2014.sav',osc_power3mHz,osc_power5mHz,osc_power10mHz,wavelet_power2d_3mHz,wavelet_power2d_5mHz,wavelet_power2d_10mHz
;	+ osc_power3mHz	= Wavelet power of filtered data
;	+ osc_power5mHz	= Wavelet power of filtered data
;	+ osc_power10mHz = Wavelet power of filtered data
;	+ wavelet_power2d_3mHz	= Wavelet power of filtered data in 2D
;	+ wavelet_power2d_5mHz	= Wavelet power of filtered data in 2D
;	+ wavelet_power2d_10mHz	= Wavelet power of filtered data in 2D
;==================================================================================================

;======================================================================================================
;			NICE PLOT OF POWER WITH PORE IMAGES FOR EACH FREQ
;======================================================================================================

;  -	Across Pore (0/45/90/135deg) at each freq.
;x1 = 300 ;binary
;x2 = 490 ;binary
x1 = 140
x2 = 675
y1 = av_yc	;449
y2 = av_yc
;print,160+32,160+139 ;As binary is smaller field than the slice.
;marsub = mardata[x1:x2,y1-50:y2+120,*]

!p.background=255.
!p.color=0.

; This works out the x and y positions of a line based on the input angle and the fact that the rotations were about the centre point!
theta = (!PI/180)*0
hyp = (x2-x1)+1
half_hyp = hyp/2.

opp1 = hyp * SIN(theta)
opp2 = half_hyp * SIN(theta)

adj1 = hyp * COS(theta)
adj2 = half_hyp * COS(theta)

roty1 = av_yc - opp2
roty2 = roty1 + opp1
rotx1 = av_xc - adj2
rotx2 = rotx1 + adj1

;marsub = mardata[x1:x1+235,y1-60:y2+80,*]
im = ROT(total_binary,135,1,av_xc-300,av_yc-(y1-50),/PIVOT)
elements = WHERE(im GT 400.)       
xy_elements = xyposition(im,elements)

rot_pore_minX = MIN(xy_elements[0,*])	;For boundary locations
rot_pore_maxX = MAX(xy_elements[0,*])	;For boundary locations

rot_pore_minY = MIN(xy_elements[1,*])	;For vertical totalling 
rot_pore_maxY = MAX(xy_elements[1,*])	;For vertical totalling

;  -	 x positions: 0deg = 195,299   45deg = 203,295,   90deg = 206,275	135deg = 195,349


mult,2,1
loadct,0,/SILENT
tvim,mardata[*,*,325]>0.6<1.7,pcharsize=1.5,xtitle='X (pix)',ytitle='Y (pix)'
loadct,3,/silent
plots,[rotx1,rotx2],[roty1,roty2],thick=3,color=155
plot,SMOOTH(WAVELET_POWER2D_5MHZ[x1:x2,y1,0],2),xst=1,yst=1,thick=3,title='Filtered Power Spectrum - 5mHz',xtitle='Distance Along Slice (pix)',ytitle='Normalised Power',charsize=1.5
loadct,1,/SILENT
;verline,rot_pore_minX,line=2,color=155,thick=2
;verline,rot_pore_maxX,line=2,color=155,thick=2
verline,195,line=2,color=155,thick=2
verline,299,line=2,color=155,thick=2
loadct,0,/SILENT
mult,1,1
X2JPEG,'Filtered_power_plot_0deg_5mHz_narrow.jpg'


