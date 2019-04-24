;+
;ROUTINE:	Velocity_table
;
;PURPOSE:	Establishes the average velocity of each BP in the tracking file and exports it to a table file
;
;USEAGE:	velocity = Velocity_table(file='output_file.dat')
;
;INPUT:		output_file.dat   =>  the desired name of the file that the calculated velocities will be exported to
;		min_existence_t  => the length of time for it to be considered as an MBP (set by user).
;
;OUTPUT:	The absolute velocity for each detected BP throughout their entire observed existence
;
;AUTHOR:	Peter H. Keys
;
;-



FUNCTION Velocity_table, file=file, min_existence_t

;readcol,'sim_track.dat', f='D,D,D,D,D,x,x', btpt_num, frame, x, y, area
;readcol,'obs_track_new.dat', f='D,D,D,D,D,x,x', btpt_num, frame, x, y, area
;readcol,'subfield_tracking', f='D,D,D,D,D,x,x', btpt_num, frame, x, y, area
;readcol,'tracking_subfield_b', f='D,D,D,D,D,x,x', btpt_num, frame, x, y, area
;readcol,'gband_4170_tracking_file', f='D,D,D,D,D,x,x', btpt_num, frame, x, y, area
;readcol,'AR_region2_tracking_table', f='D,D,D,D,D,x,x', btpt_num, frame, x, y, area
;readcol,'gb_test1_tracking_file', f='D,D,D,D,D,x,x', btpt_num, frame, x, y, area
;readcol,'ig_test2_tracks',f='D,D,D,D,D,x,x', btpt_num, frame, x, y, area
;readcol,'400G_tracks_x3',f='D,D,D,D,D,x,x', btpt_num, frame, x, y, area
;readcol,'reverse_tracking_file',f='D,D,D,D,D,x,x', btpt_num, frame, x, y, area
;readcol,'Sergiy_Gband_tracking_file',f='D,D,D,D,D,x,x', btpt_num, frame, x, y, area
;readcol,'4170_tracking_file',f='D,D,D,D,D,x,x', btpt_num, frame, x, y, area
;readcol,'sodium_tracking_file',f='D,D,D,D,D,x,x', btpt_num, frame, x, y, area
;readcol,'test3_tracking_file_sigma2',f='D,D,D,D,D,x,x', btpt_num, frame, x, y, area
;readcol,'test3_tracking_file_sigma1',f='D,D,D,D,D,x,x', btpt_num, frame, x, y, area
;readcol,'binned_tio_band_tracking_file',f='D,D,D,D,D,x,x', btpt_num, frame, x, y, area
;readcol,'IBIS_tracking_new.dat',f='D,D,D,D,D,x,x', btpt_num, frame, x, y, area
;readcol,'hinode_mbp_tracks.dat',f='D,D,D,D,D,x,x', btpt_num, frame, x, y, area
;readcol,'hinode_mbp_tracks_new.dat',f='D,D,D,D,D,x,x', btpt_num, frame, x, y, area
;readcol,'AR_above_region.dat',f='D,D,D,D,D,x,x', btpt_num, frame, x, y, area
;readcol,'AR_between_pores.dat',f='D,D,D,D,D,x,x', btpt_num, frame, x, y, area
;readcol,'Segmentation_1_tracks.dat',f='D,D,D,D,D,x,x', btpt_num, frame, x, y, area
;readcol,'QS_same_dimensions_tracking_file',f='D,D,D,D,D,x,x', btpt_num, frame, x, y, area
;readcol,'AR_mixed_polarity_tracks.dat',f='D,D,D,D,D,x,x', btpt_num, frame, x, y, area
;readcol,'Hinode_tracks_file_final.dat',f='D,D,D,D,D,x,x', btpt_num, frame, x, y, area
;readcol,'total_segmentation_track_file.dat',f='D,D,D,D,D,x,x', btpt_num, frame, x, y, area
;readcol,'Hinode_sims_tracking_file.dat',f='D,D,D,D,D,x,x', btpt_num, frame, x, y, area
;readcol,'NE_sub_tracking_file.dat',f='D,D,D,D,D,x,x', btpt_num, frame, x, y, area
readcol,'/data/rosa3/oldrosa1/phk/data/27Jul2014/tracking_file2.dat',f='I,I,I,I,I,x', btpt_num, frame, x, y, area		;NO MERGERS IN THIS DATA SET SO ONE LESS X


last_bp = MAX(btpt_num)				;The number of the last detected BP
;distance = DBL(50)				;ROSA Gband Spatial Resolution (1 pixel)
;time = 2.112
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

distance = DBL(43)				;SST WB DATA
time = 303.					;SST WB DATA


velocity_array = DBLARR((last_bp + 1))

GET_LUN, unit

OPENW,unit, ''+file+''
printf,unit,'       BP number       Initial Frame   Initial x       Initial y       Final Frame       Final x        Final y      Lifetime(s)   Velocity(km/s)'
CLOSE,unit

;==========================================================================
;	CALCULATE THE VELOCITIES AND PUT EVERYTHING IN THE OUTPUT FILE
;==========================================================================
counter = 0.

FOR i = 0, last_bp DO BEGIN &$
	BP = i  &$
	;IF (BP EQ 6042) THEN CONTINUE				;Missing BP Gband-4170 data
	;IF (BP EQ 6367) THEN CONTINUE
	;IF (BP EQ 7495) THEN CONTINUE
	;IF (BP EQ 8782) THEN CONTINUE
	;IF (BP EQ 15242) THEN CONTINUE
	;IF (BP EQ 16027) THEN CONTINUE
	;IF (BP EQ 17674) THEN CONTINUE

	;IF (BP EQ 3228) THEN CONTINUE &$			;Missing BP num in subfield search a
	
	;IF (BP EQ 9384) THEN CONTINUE &$			;Missing BP nums in 2sec data search
	;IF (BP EQ 10078) THEN CONTINUE &$
	;IF (BP EQ 16359) THEN CONTINUE &$
	;IF (BP EQ 22839) THEN CONTINUE &$

	;IF (BP EQ 7753) THEN CONTINUE &$			;Missing in reverse data set
	;IF (BP EQ 14978) THEN CONTINUE &$				
	;IF (BP EQ 20215) THEN CONTINUE &$				
	;IF (BP EQ 21469) THEN CONTINUE &$
	;IF (BP EQ 21779) THEN CONTINUE &$
	;IF (BP EQ 21820) THEN CONTINUE &$
	;IF (BP EQ 23227) THEN CONTINUE &$
	;IF (BP EQ 23882) THEN CONTINUE &$
	
	;IF (BP EQ 1602) THEN CONTINUE &$			;4170 data set for 28May2009
	;IF (BP EQ 1834) THEN CONTINUE &$
	;IF (BP EQ 3351) THEN CONTINUE &$
	;IF (BP EQ 3986) THEN CONTINUE &$
	;IF (BP EQ 4291) THEN CONTINUE &$
	;IF (BP EQ 4322) THEN CONTINUE &$
	;IF (BP EQ 6979) THEN CONTINUE &$
	;IF (BP EQ 10882) THEN CONTINUE &$
	;IF (BP EQ 11625) THEN CONTINUE &$
	;IF (BP EQ 12788) THEN CONTINUE &$
	;IF (BP EQ 18478) THEN CONTINUE &$
	;IF (BP EQ 18528) THEN CONTINUE &$
	;IF (BP EQ 20648) THEN CONTINUE &$

	;IF (BP EQ 167) THEN CONTINUE &$					;test3 sims missing BPs sigma 2
	;IF (BP EQ 206) THEN CONTINUE &$
	;IF (BP EQ 1283) THEN CONTINUE &$
	;IF (BP EQ 2207) THEN CONTINUE &$
	;IF (BP EQ 3107) THEN CONTINUE &$
	;IF (BP EQ 3218) THEN CONTINUE &$
	
	;IF (BP EQ 1213) THEN CONTINUE &$					;test3 sims missing BPs sigma 1
	
	;IF (BP EQ 216) THEN CONTINUE &$					;tio sims missing BPs sigma 2
	;IF (BP EQ 231) THEN CONTINUE &$
	
	;IF (BP EQ 820) THEN CONTINUE &$						;10Dec2011 AR11372 above AR BPs
	
	;IF (BP EQ 529) THEN CONTINUE &$						;10Dec2011 AR11371 segmentation1
	;IF (BP EQ 1078) THEN CONTINUE &$
	;IF (BP EQ 3022) THEN CONTINUE &$

	;IF (BP EQ 2386) THEN CONTINUE &$						;10Dec2011 AR11372 QS Same Dimensions

	;IF (BP EQ 1983) THEN CONTINUE &$						;10Dec2011 AR11372 Mixed Polarity Same Dimensions

	;IF (BP EQ 1175) THEN CONTINUE &$						;10Dec2011 AR11371 Total Segmentation 
	;IF (BP EQ 7100) THEN CONTINUE &$						
	;IF (BP EQ 8981) THEN CONTINUE &$
	;IF (BP EQ 8982) THEN CONTINUE &$
	
	IF (BP EQ 4) THEN CONTINUE &$							;10Dec2011 AR11371 NE SUB 
	IF (BP EQ 49) THEN CONTINUE &$
	IF (BP EQ 131) THEN CONTINUE &$
	IF (BP EQ 202) THEN CONTINUE &$
	IF (BP EQ 918) THEN CONTINUE &$
	IF (BP EQ 3076) THEN CONTINUE &$
	IF (BP EQ 3796) THEN CONTINUE &$
	IF (BP EQ 4087) THEN CONTINUE &$
	IF (BP EQ 6324) THEN CONTINUE &$
	IF (BP EQ 6521) THEN CONTINUE &$
	IF (BP EQ 12030) THEN CONTINUE &$
	IF (BP EQ 12043) THEN CONTINUE &$
	IF (BP EQ 12961) THEN CONTINUE &$
	IF (BP EQ 14401) THEN CONTINUE &$
	
	xcoords = x[WHERE (btpt_num eq BP)]  &$				;Key parameters for each BP found using MBP_tracking.pro
	ycoords = y[WHERE (btpt_num eq BP)]  &$
	frame_num = frame[WHERE (btpt_num eq BP)]  &$
	final_frame = N_ELEMENTS(frame_num) - 1  &$
	
	d = (SQRT((((xcoords[final_frame]-xcoords[0])*distance)^2)+(((ycoords[final_frame]-ycoords[0])*distance)^2))) &$
	IF (d LT (distance * 2)) THEN velocity = 0. &$
  	IF (d GE (distance * 2)) THEN velocity = ((d - (distance * 2)) / (time * (frame_num[final_frame] - frame_num[0]))) &$
	
	;velocity= (((SQRT((((xcoords[final_frame]-xcoords[0])*distance)^2)+(((ycoords[final_frame]-ycoords[0])*distance)^2))) - (2 * distance))/(time * N_ELEMENTS(frame_num))) &$
	;IF (velocity lt 0) THEN velocity = 0 &$
	
	
	frame_i = frame_num[0]  &$					;Useful numbers which will be written to a table
	xi = xcoords[0]  &$
	yi = ycoords[0]  &$
	xf = xcoords[final_frame]  &$
	yf = ycoords[final_frame]  &$
	f_frame_num = frame_num[0] + final_frame  &$
	existence_t = (N_ELEMENTS(frame_num)) * time  &$
	
	;print,btpt_num[i]  &$
	;print,frame_i  &$
	;print,xi  &$
	;print,yi  &$
	;print,f_frame_num  &$
	;print,xf  &$
	;print,yf  &$
	;print,existence_t  &$
	;print,velocity  &$
	
	IF ((velocity LE 7.) AND (existence_t GT min_existence_t)) THEN BEGIN &$
	tab = DBLARR(9,1)  &$					;Creates dimensions of table that the  results are wrote to
	tab[0,*] = BP	&$	;btpt_num[i]
	tab[1,*] = frame_i  &$
	tab[2,*] = xi  &$
	tab[3,*] = yi  &$
	tab[4,*] = f_frame_num  &$
	tab[5,*] = xf  &$
	tab[6,*] = yf  &$
	tab[7,*] = existence_t  &$
	tab[8,*] = velocity  &$
	;IF ((velocity GT 5.) AND (velocity LT 7.) AND (existence_t GT min_existence_t)) THEN print,BP,velocity &$ ;TEST LINE

		OPENU,unit,''+file+'',/APPEND,WIDTH=250  &$	

		IF (existence_t LE min_existence_t) THEN BEGIN  &$
			 CLOSE,unit &$
		ENDIF  &$

		IF (existence_t GT min_existence_t) THEN BEGIN &$
			printf,unit,tab  &$						;Results written into table named in command line
			CLOSE,unit  &$
		ENDIF  &$
	ENDIF &$

	velocity_array[i] = velocity  &$				;The array of the velocity values which is used for the histogram

	IF ((velocity GE 2) AND (velocity LE 7)) THEN counter = counter + 1 &$
	average_velocity = TOTAL(velocity_array)/N_ELEMENTS(velocity_array)  &$
ENDFOR

;==========================================================================
;	CREATE A HISTOGRAM AND CALCULATE VITAL STATISTICS
;==========================================================================

readcol,''+file+'',f='x,x,x,x,x,x,x,D,D',Lifetime,BP_velocity,SKIPLINE=1

!p.background=255.
!p.color=0.

vel_histogram = HISTOGRAM(BP_velocity,binsize=5,max=20)		;Makes array of the density function of a set of values

column_names=STRARR(5)						;Names for the Histogram's columns
column_names[0]='0'
column_names[1]='1-5'
column_names[2]='5-10'
column_names[3]='10-15'
column_names[4]='15-20'

colour_index = INTARR(5)					;Colour index for the bars of the histogram
colour_index[0]=35   
colour_index[1]=120
colour_index[2]=170
colour_index[3]=220
colour_index[4]=240

loadct,3,/silent

bar_plot,vel_histogram,background=255,barnames=column_names,barspace=0.1,color=colour_index,title='Distribution of BP Velocities',xtitle='Velocity (km/s)',ytitle='Frequency'			;Plots the histogram using the density function of vel_histogram

av = MEAN(BP_velocity)
mean_deviation = MEANABSDEV(BP_velocity)
median_value = MEDIAN(BP_velocity)
median_deviation = MEANABSDEV(BP_velocity,/MEDIAN)
max_vel = MAX(BP_velocity)
min_vel = MIN(BP_velocity)
range_vel = max_vel - min_vel
av_lifetime = MEAN(Lifetime)
max_lifetime = MAX(Lifetime)
min_lifetime = MIN(Lifetime)
num_of_bps = N_ELEMENTS(BP_velocity)
percentage_bps_studied = (num_of_bps/(last_bp+1))*100.

;==========================================================================
;	PRINT THE STATISTICS TO THE SCREEN OF THE IDL SESSION
;==========================================================================

print,'  '
print,'  '
print,'Number of MBPs Found = ' +arr2str(num_of_bps,/trim)
print,'  '
print,'  '

print,'Average Velocity = ' + arr2str(av,/trim)
print,'Mean Deviation = ' + arr2str(mean_deviation,/trim)
print,'Median = ' + arr2str(median_value,/trim)
print,'Deviation in the Median = ' + arr2str(median_deviation,/trim)

print,'  '

print,'Max Velocity = ' + arr2str(max_vel,/trim)
print,'Min Velocity = ' + arr2str(min_vel,/trim)
print,'Velocity Range = ' + arr2str(range_vel,/trim)

print,'  '
print,'  '
print,'Mean MBP Lifetime = ' + arr2str(av_lifetime,/trim)
print,'Max MBP Lifetime = ' + arr2str(max_lifetime,/trim)
print,'Min MBP Lifetime = ' + arr2str(min_lifetime,/trim)
print,'  '
print,'  '
print,'Number of BPs with velocity > 2km/s = ' + arr2str(counter,/trim)
print,'% of BPs with velocity > 2km/s = ' + arr2str((counter/num_of_bps)*100,/trim)
print,'  '
print,'  '
loadct,0,/silent

;==========================================================================
;	PRINT THE STATISTICS TO A FILE WITH A SIMILAR FILE PATH
;==========================================================================

OPENW,unit,''+file+'_statistics'
printf,unit,'Total Number of MBPs tracked with MBP_tracking.pro = ' +arr2str(last_bp+1,/trim)
printf,unit,'  '
printf,unit,'Number of MBPs Found which exist longer than '+arr2str(min_existence_t,/trim) +'s = ' +arr2str(num_of_bps,/trim)
printf,unit,'  '
printf,unit,'% of MBPs studied out of Total MBPs found = ' +arr2str(percentage_bps_studied,/trim) +'%'
printf,unit,'  '
printf,unit,'------------------------------------------------------------------------------------------------------'
printf,unit,'STATISTICS based on the ' +arr2str(num_of_bps,/trim) +' MBPs analysed:'
printf,unit,'  '

printf,unit,'Average Velocity = ' + arr2str(av,/trim)
printf,unit,'Mean Deviation = ' + arr2str(mean_deviation,/trim)
printf,unit,'Median = ' + arr2str(median_value,/trim)
printf,unit,'Deviation in the Median = ' + arr2str(median_deviation,/trim)

printf,unit,'  '

printf,unit,'Max Velocity = ' + arr2str(max_vel,/trim)
printf,unit,'Min Velocity = ' + arr2str(min_vel,/trim)
printf,unit,'Velocity Range = ' + arr2str(range_vel,/trim)

printf,unit,'  '
printf,unit,'  '
printf,unit,'Mean MBP Lifetime = ' + arr2str(av_lifetime,/trim)
printf,unit,'Max MBP Lifetime = ' + arr2str(max_lifetime,/trim)
printf,unit,'Min MBP Lifetime = ' + arr2str(min_lifetime,/trim)
printf,unit,'  '
printf,unit,'  '
printf,unit,'Number of BPs with velocity > 2km/s = ' + arr2str(counter,/trim)
printf,unit,'% of BPs with velocity > 2km/s = ' + arr2str((counter/num_of_bps)*100,/trim)
printf,unit,'  '
printf,unit,'  '
CLOSE,unit

FREE_LUN,unit
RETURN, average_velocity

END
