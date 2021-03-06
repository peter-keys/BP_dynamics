; GROSS DIFFUSIVITY CALCULATION

FUNCTION gross_diffusivity, bin_time

restore,'/data/rosa3/oldrosa1/phk/data/10Dec2011/AR11372/AR_between_pores_diff_indices_parameters.sav',/ver
restore,'/data/rosa3/oldrosa1/phk/data/10Dec2011/AR11372/QS_diff_indices_parameters.sav',/ver
restore,'/data/rosa3/oldrosa1/phk/data/10Dec2011/AR11372/AR_mixed_polarity_indices_parameters.sav',/ver

min_existence_t = 30.					; The min existence time
;bin_time = 10.						; The time for the bins to establish the sd
;readcol,'/data/rosa3/oldrosa1/phk/data/10Dec2011/AR11372/AR_mixed_polarity_tracks.dat',f='D,D,D,D,D,x,x', btpt_num, frame, x, y, area
readcol,'/data/rosa3/oldrosa1/phk/data/10Dec2011/AR11372/QS_same_dimensions_tracking_file',f='D,D,D,D,D,x,x', btpt_num, frame, x, y, area
;readcol,'/data/rosa3/oldrosa1/phk/data/10Dec2011/AR11372/AR_between_pores.dat',f='D,D,D,D,D,x,x', btpt_num, frame, x, y, area

last_bp = MAX(btpt_num)					;The number of the last detected BP
distance = DBL(50)					;ROSA Gband Spatial Resolution (1 pixel)
time = 2.112

;no_of_bins = ROUND((MAX(indices_pore[0,*]) - MIN(indices_pore[0,*]))/bin_time)
no_of_bins = ROUND((MAX(indices_qs[0,*]) - MIN(indices_qs[0,*]))/bin_time)
;no_of_bins = ROUND((MAX(indices_mixed[0,*]) - MIN(indices_mixed[0,*]))/bin_time)

sd_array = DBLARR(no_of_bins+1)
elements_per_bin = DBLARR(no_of_bins+1)
sd_array[*] = 0.
elements_per_bin[*] = 0.

FOR i = 0, last_bp DO BEGIN &$
	BP = i  &$
	
	IF (BP EQ 2386) THEN CONTINUE &$						;10Dec2011 AR11372 QS Same Dimensions

	;IF (BP EQ 1983) THEN CONTINUE &$						;10Dec2011 AR11372 Mixed Polarity Same Dimensions

 	 ;IF (BP EQ 51) THEN CONTINUE &$			;THESE ARE BPs that have 0 gradient 
	 ;IF (BP EQ 351) THEN CONTINUE &$
	 ;IF (BP EQ 396) THEN CONTINUE &$			;THUS CRASHES THE CODE AND HAVE BEEN 
	 ;IF (BP EQ 459) THEN CONTINUE &$			;MANUALLY REMOVED
	 ;IF (BP EQ 586) THEN CONTINUE &$
	 ;IF (BP EQ 594) THEN CONTINUE &$			;NEED TO ESTABLISH A WAY OF RECTIFYING THIS
	 ;IF (BP EQ 655) THEN CONTINUE &$			;NEED MORE TIME THOUGH WHICH I DON'T HAVE
	 ;IF (BP EQ 1200) THEN CONTINUE &$			;AT THE MINUTE
	 ;IF (BP EQ 1677) THEN CONTINUE &$
	 ;IF (BP EQ 1716) THEN CONTINUE &$
	 ;IF (BP EQ 1719) THEN CONTINUE &$
	 ;IF (BP EQ 1821) THEN CONTINUE &$			;10Dec2011 AR11372 AR between Pores data
	 ;IF (BP EQ 1862) THEN CONTINUE &$
 	 ;IF (BP EQ 1978) THEN CONTINUE &$
	 ;IF (BP EQ 2076) THEN CONTINUE &$
	 ;IF (BP EQ 2079) THEN CONTINUE &$
	 ;IF (BP EQ 2120) THEN CONTINUE &$
	 ;IF (BP EQ 2132) THEN CONTINUE &$
	 ;IF (BP EQ 2487) THEN CONTINUE &$
	 ;IF (BP EQ 2678) THEN CONTINUE &$
	 ;IF (BP EQ 2678) THEN CONTINUE &$
	 ;IF (BP EQ 2750) THEN CONTINUE &$
	 ;IF (BP EQ 2938) THEN CONTINUE &$
	 ;IF (BP EQ 2959) THEN CONTINUE &$
	 ;IF (BP EQ 3160) THEN CONTINUE &$
	 ;IF (BP EQ 4653) THEN CONTINUE &$
	 ;IF (BP EQ 4761) THEN CONTINUE &$
	 ;IF (BP EQ 4764) THEN CONTINUE &$
	 ;IF (BP EQ 4838) THEN CONTINUE &$
	 ;IF (BP EQ 5094) THEN CONTINUE &$
	 ;IF (BP EQ 5370) THEN CONTINUE &$
	 ;IF (BP EQ 5430) THEN CONTINUE &$
	 ;IF (BP EQ 5745) THEN CONTINUE &$
	 ;IF (BP EQ 5970) THEN CONTINUE &$
	 ;IF (BP EQ 6000) THEN CONTINUE &$
	 ;IF (BP EQ 6030) THEN CONTINUE &$
	 ;IF (BP EQ 6112) THEN CONTINUE &$
	 ;IF (BP EQ 6313) THEN CONTINUE &$
	 ;IF (BP EQ 6742) THEN CONTINUE &$
	 ;IF (BP EQ 6781) THEN CONTINUE &$
	 ;IF (BP EQ 7080) THEN CONTINUE &$
	 ;IF (BP EQ 7131) THEN CONTINUE &$
	 ;IF (BP EQ 7285) THEN CONTINUE &$
 	 ;IF (BP EQ 7333) THEN CONTINUE &$
 	 ;IF (BP EQ 7438) THEN CONTINUE &$
 	 ;IF (BP EQ 7557) THEN CONTINUE &$
 	 ;IF (BP EQ 7615) THEN CONTINUE &$
 	 ;IF (BP EQ 7748) THEN CONTINUE &$
 	 ;IF (BP EQ 7781) THEN CONTINUE &$
 	 ;IF (BP EQ 8412) THEN CONTINUE &$
 	 ;IF (BP EQ 8418) THEN CONTINUE &$
	
	 IF (BP EQ 72) THEN CONTINUE &$				;10Dec2011 AR11372 QS Region
	 IF (BP EQ 1550) THEN CONTINUE &$
	 IF (BP EQ 2158) THEN CONTINUE &$
	 IF (BP EQ 2382) THEN CONTINUE &$
	 IF (BP EQ 2386) THEN CONTINUE &$
	 IF (BP EQ 2403) THEN CONTINUE &$
	 IF (BP EQ 2473) THEN CONTINUE &$
	 IF (BP EQ 2802) THEN CONTINUE &$

	 ;IF (BP EQ 230) THEN CONTINUE &$				;10Dec2011 AR11372 AR Mixed Region
	 ;IF (BP EQ 239) THEN CONTINUE &$
	 ;IF (BP EQ 344) THEN CONTINUE &$
	 ;IF (BP EQ 456) THEN CONTINUE &$
	 ;IF (BP EQ 601) THEN CONTINUE &$
	 ;IF (BP EQ 648) THEN CONTINUE &$
	 ;IF (BP EQ 772) THEN CONTINUE &$
	 ;IF (BP EQ 808) THEN CONTINUE &$
	 ;IF (BP EQ 1425) THEN CONTINUE &$
	 ;IF (BP EQ 1430) THEN CONTINUE &$
	 ;IF (BP EQ 1629) THEN CONTINUE &$
	 ;IF (BP EQ 1652) THEN CONTINUE &$
	 ;IF (BP EQ 1655) THEN CONTINUE &$
	 ;IF (BP EQ 1758) THEN CONTINUE &$
	 ;IF (BP EQ 2018) THEN CONTINUE &$
	 ;IF (BP EQ 2686) THEN CONTINUE &$
	 ;IF (BP EQ 2889) THEN CONTINUE &$
	 ;IF (BP EQ 2961) THEN CONTINUE &$
	 ;IF (BP EQ 3003) THEN CONTINUE &$
	 ;IF (BP EQ 3025) THEN CONTINUE &$
	 ;IF (BP EQ 3559) THEN CONTINUE &$
	 ;IF (BP EQ 3598) THEN CONTINUE &$
	 ;IF (BP EQ 3811) THEN CONTINUE &$
	 ;IF (BP EQ 3860) THEN CONTINUE &$
	 ;IF (BP EQ 3896) THEN CONTINUE &$
	 ;IF (BP EQ 3954) THEN CONTINUE &$
	 ;IF (BP EQ 3992) THEN CONTINUE &$
	 ;IF (BP EQ 4040) THEN CONTINUE &$
	 ;IF (BP EQ 4183) THEN CONTINUE &$
	 ;IF (BP EQ 4583) THEN CONTINUE &$
	 ;IF (BP EQ 5330) THEN CONTINUE &$
	 ;IF (BP EQ 5487) THEN CONTINUE &$
	 ;IF (BP EQ 5582) THEN CONTINUE &$
	 ;IF (BP EQ 5777) THEN CONTINUE &$
	 ;IF (BP EQ 5986) THEN CONTINUE &$
	 ;IF (BP EQ 6038) THEN CONTINUE &$
	 ;IF (BP EQ 6100) THEN CONTINUE &$
	 ;IF (BP EQ 6131) THEN CONTINUE &$
	 ;IF (BP EQ 6239) THEN CONTINUE &$
	 ;IF (BP EQ 6348) THEN CONTINUE &$
	 ;IF (BP EQ 6892) THEN CONTINUE &$
	 ;IF (BP EQ 7106) THEN CONTINUE &$
	 ;IF (BP EQ 7127) THEN CONTINUE &$
	 ;IF (BP EQ 7132) THEN CONTINUE &$
	 ;IF (BP EQ 7402) THEN CONTINUE &$
	 ;IF (BP EQ 7833) THEN CONTINUE &$
	 ;IF (BP EQ 7893) THEN CONTINUE &$
	 ;IF (BP EQ 7913) THEN CONTINUE &$
	 ;IF (BP EQ 8943) THEN CONTINUE &$
	 ;IF (BP EQ 9220) THEN CONTINUE &$
	 ;IF (BP EQ 9232) THEN CONTINUE &$
	 ;IF (BP EQ 9305) THEN CONTINUE &$
	 ;IF (BP EQ 9387) THEN CONTINUE &$
	 ;IF (BP EQ 9470) THEN CONTINUE &$
	 ;IF (BP EQ 9498) THEN CONTINUE &$
	 ;IF (BP EQ 9507) THEN CONTINUE &$
	 ;IF (BP EQ 9704) THEN CONTINUE &$
	 ;IF (BP EQ 9848) THEN CONTINUE &$
	
	xcoords = x[WHERE (btpt_num eq BP)]  &$				
	ycoords = y[WHERE (btpt_num eq BP)]  &$
	frame_num = frame[WHERE (btpt_num eq BP)]  &$
	final_frame = N_ELEMENTS(frame_num) - 1  &$
	
	d = (SQRT((((xcoords[final_frame]-xcoords[0])*distance)^2)+(((ycoords[final_frame]-ycoords[0])*distance)^2))) &$
	IF (d LT (distance * 2)) THEN velocity = 0. &$
  	IF (d GE (distance * 2)) THEN velocity = ((d - (distance * 2)) / (time * (frame_num[final_frame] - frame_num[0]))) &$
	frame_i = frame_num[0]  &$					
	xi = xcoords[0]  &$
	yi = ycoords[0]  &$
	xf = xcoords[final_frame]  &$
	yf = ycoords[final_frame]  &$
	f_frame_num = frame_num[0] + final_frame  &$
	existence_t = (N_ELEMENTS(frame_num)) * time  &$
	time_frame = (FINDGEN(N_ELEMENTS(FRAME_NUM-1))+1)*time &$
	
	IF ((velocity LE 7.) AND (existence_t GT min_existence_t)) THEN BEGIN &$
		no_of_bins_mbp = ROUND(((f_frame_num * time) - (frame_i * time) )/bin_time) &$		;THE NUMBER OF BINS OCCUPIED BY THE CURRENT MBP
		FOR k = 0, (no_of_bins_mbp - 1) DO BEGIN &$
			frame_of_bin = ROUND(k * (bin_time/time)) &$
			sd_of_bin = (((xcoords[frame_of_bin]-xcoords[0])*distance)^2)+(((ycoords[frame_of_bin]-ycoords[0])*distance)^2) &$
			sd_array[k] = sd_array[k] + sd_of_bin &$
			elements_per_bin[k] = elements_per_bin[k] + 1. &$
		ENDFOR &$
	ENDIF &$
ENDFOR

mean_sd_per_bin = sd_array/elements_per_bin

time_array_binned = FINDGEN(N_ELEMENTS(sd_array))*bin_time

resultant = [mean_sd_per_bin,elements_per_bin]

RETURN,mean_sd_per_bin
END
