;********************************************************************************
;			DIFFUSION INDEX CALCULATION
;
;********************************************************************************
;	+ Program for calculating diffusion index of MBPs
;	+ Take output of tracking code
;	+ Calculates the squared displacement (sd) for each MBP
; 	+ Outputs a file with values of sd at each time-step
;	+ Reads this in to establish the average over all MBPs
;	+ Plot as log-log and gamma is found by the slope of the graph
;	+ Gamma gives the diffusive properties
;	+ The values of the plot are the final outputs (with the file for each MBPs sd)
;
;	P.Keys 
;
;	USE:
;	.r /home/phk/idl/BP_dynamics/diff_index.pro
;	indices = DIFF_INDEX(file='AR_between_pores_diff_index.dat',30.)


FUNCTION DIFF_INDEX,  file=file, min_existence_t

;readcol,'/data/rosa3/oldrosa1/phk/data/10Dec2011/AR11372/AR_above_region.dat',f='D,D,D,D,D,x,x', btpt_num, frame, x, y, area
;readcol,'/data/rosa3/oldrosa1/phk/data/10Dec2011/AR11372/AR_between_pores.dat',f='D,D,D,D,D,x,x', btpt_num, frame, x, y, area
;readcol,'Segmentation_1_tracks.dat',f='D,D,D,D,D,x,x', btpt_num, frame, x, y, area
;readcol,'/data/rosa3/oldrosa1/phk/data/10Dec2011/AR11372/QS_same_dimensions_tracking_file',f='D,D,D,D,D,x,x', btpt_num, frame, x, y, area
readcol,'/data/rosa3/oldrosa1/phk/data/10Dec2011/AR11372/AR_mixed_polarity_tracks.dat',f='D,D,D,D,D,x,x', btpt_num, frame, x, y, area
;readcol,'Hinode_tracks_file_final.dat',f='D,D,D,D,D,x,x', btpt_num, frame, x, y, area
;readcol,'total_segmentation_track_file.dat',f='D,D,D,D,D,x,x', btpt_num, frame, x, y, area
;readcol,'Hinode_sims_tracking_file.dat',f='D,D,D,D,D,x,x', btpt_num, frame, x, y, area

last_bp = MAX(btpt_num)				;The number of the last detected BP
distance = DBL(50)				;ROSA Gband Spatial Resolution (1 pixel)
time = 2.112
;time = 2.1
;time = 2.56					;Cadence for this particular data set
;time = 0.528					;Cadence of Phil's reconstructed Gband data
;distance = DBL(25)				;Simulated Gband data
;time = 16.					;Simulated Gband data
;time = 8.7					;Simulated Gband data
;time = 6.7					;Simulated Gband data

;distance = DBL(60)				;IBIS Data for 28May2009
;time=39.703					;IBIS Data for 28May2009
;distance = DBL(79)				;Hinode Data
;time=31.94					;Hinode Data

;distance = DBL(50)				;test3 sim values
;time =2.1

GET_LUN, unit

OPENW,unit, ''+file+''
printf,unit,'       BP number       Initial Frame      Frame of Calc.      Squared Displacement       '
CLOSE,unit

OPENW,unit, ''+file+'_diff_index_vals'
printf,unit,'       BP number       Mean Sq. Disp   Gamma (grad)      Uncert. Gamma     Intercept    Uncert. Intercept       '
CLOSE,unit


;==========================================================================
;	CALCULATE THE VELOCITIES AND PUT EVERYTHING IN THE OUTPUT FILE
;==========================================================================

blounter = 0.
;flounter = 0
FOR i = 0, last_bp DO BEGIN &$
	BP = i  &$
	
	;IF (BP EQ 820) THEN CONTINUE &$						;10Dec2011 AR11372 above AR BPs
	
	;IF (BP EQ 529) THEN CONTINUE &$						;10Dec2011 AR11371 segmentation1
	;IF (BP EQ 1078) THEN CONTINUE &$
	;IF (BP EQ 3022) THEN CONTINUE &$

	;IF (BP EQ 2386) THEN CONTINUE &$						;10Dec2011 AR11372 QS Same Dimensions

	IF (BP EQ 1983) THEN CONTINUE &$						;10Dec2011 AR11372 Mixed Polarity Same Dimensions

	;IF (BP EQ 1175) THEN CONTINUE &$						;10Dec2011 AR11371 Total Segmentation 
	;IF (BP EQ 7100) THEN CONTINUE &$						
	;IF (BP EQ 8981) THEN CONTINUE &$
	;IF (BP EQ 8982) THEN CONTINUE &$
	
	
	 ;IF (BP EQ 11) THEN CONTINUE &$
	 ;IF (BP EQ 14) THEN CONTINUE &$
	 ;IF (BP EQ 23) THEN CONTINUE &$
	 ;IF (BP EQ 33) THEN CONTINUE &$
	 
	 ;gamma = 1. &$
	 
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
	
	 ;IF (BP EQ 72) THEN CONTINUE &$				;10Dec2011 AR11372 QS Region
	 ;IF (BP EQ 1550) THEN CONTINUE &$
	 ;IF (BP EQ 2158) THEN CONTINUE &$
	 ;IF (BP EQ 2382) THEN CONTINUE &$
	 ;IF (BP EQ 2386) THEN CONTINUE &$
	 ;IF (BP EQ 2403) THEN CONTINUE &$
	 ;IF (BP EQ 2473) THEN CONTINUE &$
	 ;IF (BP EQ 2802) THEN CONTINUE &$

	 IF (BP EQ 230) THEN CONTINUE &$				;10Dec2011 AR11372 AR Mixed Region
	 IF (BP EQ 239) THEN CONTINUE &$
	 IF (BP EQ 344) THEN CONTINUE &$
	 IF (BP EQ 456) THEN CONTINUE &$
	 IF (BP EQ 601) THEN CONTINUE &$
	 IF (BP EQ 648) THEN CONTINUE &$
	 IF (BP EQ 772) THEN CONTINUE &$
	 IF (BP EQ 808) THEN CONTINUE &$
	 IF (BP EQ 1425) THEN CONTINUE &$
	 IF (BP EQ 1430) THEN CONTINUE &$
	 IF (BP EQ 1629) THEN CONTINUE &$
	 IF (BP EQ 1652) THEN CONTINUE &$
	 IF (BP EQ 1655) THEN CONTINUE &$
	 IF (BP EQ 1758) THEN CONTINUE &$
	 IF (BP EQ 2018) THEN CONTINUE &$
	 IF (BP EQ 2686) THEN CONTINUE &$
	 IF (BP EQ 2889) THEN CONTINUE &$
	 IF (BP EQ 2961) THEN CONTINUE &$
	 IF (BP EQ 3003) THEN CONTINUE &$
	 IF (BP EQ 3025) THEN CONTINUE &$
	 IF (BP EQ 3559) THEN CONTINUE &$
	 IF (BP EQ 3598) THEN CONTINUE &$
	 IF (BP EQ 3811) THEN CONTINUE &$
	 IF (BP EQ 3860) THEN CONTINUE &$
	 IF (BP EQ 3896) THEN CONTINUE &$
	 IF (BP EQ 3954) THEN CONTINUE &$
	 IF (BP EQ 3992) THEN CONTINUE &$
	 IF (BP EQ 4040) THEN CONTINUE &$
	 IF (BP EQ 4183) THEN CONTINUE &$
	 IF (BP EQ 4583) THEN CONTINUE &$
	 IF (BP EQ 5330) THEN CONTINUE &$
	 IF (BP EQ 5487) THEN CONTINUE &$
	 IF (BP EQ 5582) THEN CONTINUE &$
	 IF (BP EQ 5777) THEN CONTINUE &$
	 IF (BP EQ 5986) THEN CONTINUE &$
	 IF (BP EQ 6038) THEN CONTINUE &$
	 IF (BP EQ 6100) THEN CONTINUE &$
	 IF (BP EQ 6131) THEN CONTINUE &$
	 IF (BP EQ 6239) THEN CONTINUE &$
	 IF (BP EQ 6348) THEN CONTINUE &$
	 IF (BP EQ 6892) THEN CONTINUE &$
	 IF (BP EQ 7106) THEN CONTINUE &$
	 IF (BP EQ 7127) THEN CONTINUE &$
	 IF (BP EQ 7132) THEN CONTINUE &$
	 IF (BP EQ 7402) THEN CONTINUE &$
	 IF (BP EQ 7833) THEN CONTINUE &$
	 IF (BP EQ 7893) THEN CONTINUE &$
	 IF (BP EQ 7913) THEN CONTINUE &$
	 IF (BP EQ 8943) THEN CONTINUE &$
	 IF (BP EQ 9220) THEN CONTINUE &$
	 IF (BP EQ 9232) THEN CONTINUE &$
	 IF (BP EQ 9305) THEN CONTINUE &$
	 IF (BP EQ 9387) THEN CONTINUE &$
	 IF (BP EQ 9470) THEN CONTINUE &$
	 IF (BP EQ 9498) THEN CONTINUE &$
	 IF (BP EQ 9507) THEN CONTINUE &$
	 IF (BP EQ 9704) THEN CONTINUE &$
	 IF (BP EQ 9848) THEN CONTINUE &$
	
	xcoords = x[WHERE (btpt_num eq BP)]  &$				
	ycoords = y[WHERE (btpt_num eq BP)]  &$
	frame_num = frame[WHERE (btpt_num eq BP)]  &$
	final_frame = N_ELEMENTS(frame_num) - 1  &$
	
	d = (SQRT((((xcoords[final_frame]-xcoords[0])*distance)^2)+(((ycoords[final_frame]-ycoords[0])*distance)^2))) &$
	IF (d LT (distance * 2)) THEN velocity = 0. &$
  	IF (d GE (distance * 2)) THEN velocity = ((d - (distance * 2)) / (time * (frame_num[final_frame] - frame_num[0]))) &$
	
	;velocity= (((SQRT((((xcoords[final_frame]-xcoords[0])*distance)^2)+(((ycoords[final_frame]-ycoords[0])*distance)^2))) - (2 * distance))/(time * N_ELEMENTS(frame_num))) &$
	;IF (velocity lt 0) THEN velocity = 0 &$
	
	
	frame_i = frame_num[0]  &$					
	xi = xcoords[0]  &$
	yi = ycoords[0]  &$
	xf = xcoords[final_frame]  &$
	yf = ycoords[final_frame]  &$
	f_frame_num = frame_num[0] + final_frame  &$
	existence_t = (N_ELEMENTS(frame_num)) * time  &$
	time_frame = (FINDGEN(N_ELEMENTS(FRAME_NUM-1))+1)*time &$
		
	IF ((velocity LE 7.) AND (existence_t GT min_existence_t)) THEN BEGIN &$
	sd = DBLARR(N_ELEMENTS(frame_num)-1) &$
	 FOR k = 1, (N_ELEMENTS(frame_num)-1) DO BEGIN &$
	  sd[k-1] = (((xcoords[k]-xcoords[0])*distance)^2)+(((ycoords[k]-ycoords[0])*distance)^2) &$
 	  sd_k = (((xcoords[k]-xcoords[0])*distance)^2)+(((ycoords[k]-ycoords[0])*distance)^2) &$
	  tab = DBLARR(4,1)  &$					
	  tab[0,*] = BP	&$	;btpt_num[i]
	  tab[1,*] = frame_num[0] &$
	  tab[2,*] = frame_num[k] &$
	  tab[3,*] = sd_k &$
	  OPENU,unit,''+file+'',/APPEND,WIDTH=250  &$	
	   printf,unit,tab  &$						;Results written into table named in command line
	  CLOSE,unit  &$
	 ENDFOR &$
	 stdev_sd = STDDEV(sd) &$
	 ;IF (stdev_sd LT 1.) THEN CONTINUE &$
	 log_tau = ALOG10(time_frame) &$
	 log_sd = ALOG10(sd) &$
 	FOR l = 0, N_ELEMENTS(log_sd)-1 DO BEGIN &$
		nan_changer = FINITE(log_sd[l]) &$
		IF (nan_changer EQ 1) THEN log_sd[l] = log_sd[l] &$
		IF (nan_changer EQ 0) THEN log_sd[l] = 1. &$
	ENDFOR &$ 
	
	 ;IF (BP EQ 351) THEN sd = sd[6:29] &$
	 
	; IF (MEAN(sd) LT 10000) THEN flounter = flounter + 1 &$
	 IF (MEAN(sd) LT 10000) THEN  CONTINUE &$
	 
	 ;coeff = ROBUST_LINEFIT(time_frame,sd,sfit,stdev_sd,dcoeff) &$
	 ;ON_IOERROR,JUMP &$
	 
	 coeff = ROBUST_LINEFIT(log_tau,log_sd,sfit,stdev_sd,dcoeff) &$
	 ;coeff[0] = intercept// coeff[1] = gradient// sfit = line-fit points // dcoeff[1] = uncertainty in grad
	 
	 ;LOADCT,0,/SILENT &$
	 ;window,0,xsize=1000,ysize=800 &$
	 ;PLOT,time_frame,sd,/XLOG,/YLOG,xtitle='Time (s)',ytitle='Sq. Displacement',yst=1,xst=1,title='BP Number = '+arr2str(BP,/trim),charsize=2.4,psym=3,charthick=1.5 &$
	 ;LOADCT,3,/SILENT &$
	 ;PLOTS,[MIN(time_frame),MAX(time_frame)],[MIN(WHERE[sfit GE 0]),MAX(WHERE[sfit GE 0])],color=155,thick=2,line=2 &$
	 ;OPLOT,time_frame,sfit,color=155 &$
	 LOADCT,0,/SILENT &$
	 
	 change_a_nan = FINITE(((coeff[1]))) &$
	 ;change_a_nan = FINITE((ALOG10(coeff[1]))) &$
	 ;IF (change_a_nan EQ 0) THEN coeff[1] = 0. &$ 
	 
	 IF (change_a_nan EQ 1) THEN (gamma = coeff[1]) ELSE (gamma = 1.) &$
	 IF (gamma LT 0) THEN gamma = gamma * (-1) &$
	 ;IF (coeff[1] LT 1) AND (coeff[1] GT 0) THEN gamma = 0.
	 
	 ;JUMP: print,BP &$
	 tabular_sd = DBLARR(7,1) &$
	 tabular_sd[0,*] = BP &$
	 tabular_sd[1,*] = MEAN(sd) &$
	 tabular_sd[2,*] = gamma &$
	 tabular_sd[3,*] = dcoeff[1] &$
	 tabular_sd[4,*] = coeff[0] &$
	 tabular_sd[5,*] = dcoeff[0] &$
	 tabular_sd[6,*] = existence_t &$
	 OPENW,unit,''+file+'_diff_index_vals',/APPEND,WIDTH=250
 		printf,unit,tabular_sd  &$						
	 CLOSE,unit  &$
	 ;av_sd_tau[0,blounter] = [existence_t, MEAN(sd),gamma] &$
	 blounter = blounter + 1. &$
	ENDIF &$
	;USE_ROBUST_LINEFIT.PRO TO WORK OUT GRADIENTS OF LOG-LOG GRAPH

ENDFOR
;print,flounter

readcol,''+file+'_diff_index_vals',f='x,D,D,D,D,x,D',mean_sq_disp,gamma_val,gamma_uncert,intercept,time_span,SKIPLINE=1

av_sd_tau = DBLARR(5,(N_ELEMENTS(mean_sq_disp)))

diffusion_coeff = DBLARR(N_ELEMENTS(intercept))

intercept_new = DBLARR(N_ELEMENTS(intercept))

FOR i = 0, N_ELEMENTS(intercept)-1 DO BEGIN &$
	
	;log_sd = ALOG10(mean_sq_disp[i]) &$
	;log_tau = ALOG10(time_span[i]) &$
	;coeff = ROBUST_LINEFIT(log_tau,log_sd,sfit,stdev_sd,dcoeff) &$
	;intercept_new[i] = coeff[0] &$
	;diffusion_coeff[i] = (((coeff[0]*coeff[1]/4.)*(time_span[i]^(coeff[1]-1)))) &$
	;print,coeff[1],diffusion_coeff[i] &$
	
	;intercept_new[i] = 10^(ABS(intercept[i])) &$
	;intercept_new[i] = mean_sq_disp[i]/((time_span[i])^(gamma_val[i])) &$
	;diffusion_coeff[i] = (((intercept_new[i])*gamma_val[i])/4.)*(time_span[i]^(gamma_val[i]-1)) &$
ENDFOR

FOR i = 0, N_ELEMENTS(mean_sq_disp)-1 DO BEGIN &$
	av_sd_tau[0,i] = time_span[i] &$
	av_sd_tau[1,i] = mean_sq_disp[i] &$
	av_sd_tau[2,i] = gamma_val[i] &$
	av_sd_tau[3,i] = gamma_uncert[i] &$
	av_sd_tau[4,i] = diffusion_coeff[i] &$
ENDFOR

mult,2,1
window,0,xsize=1400,ysize=800
PLOT,av_sd_tau[0,*],av_sd_tau[1,*],/XLOG,/YLOG,charsize=2.4,CHARTHICK=1.5,PSYM=3,THICK=3,TITLE='Overall MBP Diffusion'
print,' '
print,'Mean Gamma (diffusion index) = ' +arr2str(MEAN(av_sd_tau[2,*]),/trim)
print,' '
print,'Error in Gamma = ' +arr2str(STDDEV(av_sd_tau[2,*])/2.,/trim)	
print,' '
print,'Diffusion coefficient (D) = ' +arr2str(MEAN(av_sd_tau[4,*]),/trim)
print,' ' 
print,'Mean Intercept = '+arr2str(MEAN(intercept),/trim)
print,' '
print,'Error in D = ' +arr2str(STDDEV(av_sd_tau[4,*])/2.,/trim)	
	

index_hist = HISTOGRAM(av_sd_tau[2,*],binsize=0.1)
index_hist = index_hist/(MAX(index_hist))

PLOT,index_hist,psym=10,XTITLE='Diffusion Index (gamma)',YTITLE='Normalised Freq.',thick=2.5,charsize=2.4,charthick=1.5

mult,1,1
FREE_LUN,unit
RETURN, av_sd_tau

END
