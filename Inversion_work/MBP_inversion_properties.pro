;+
;ROUTINE:	MBP_inversion_properties.pro
;
;PURPOSE:	Reads in the outputs of NICOLE inversions and 
;		the tracking algoritm results to make a file
; 		with all the info in the one place (i.e., file)		
;
;USEAGE:	MBP_inversion_properties,file=file
;
;INPUT:		Code reads in the output of the Inversions 
;		and the tracking code. You need to augment 
;		the actual code to get it so that you have 
; 		the right paths to the relevant data sets 
;		etc... 
;		
;		You need to input a file name to store 
;		the info out to.
;
;OUTPUT:	file containg all the relevant info from 
; 		the Inversions and tracking arrays. File 
;		will contain the following variables:
;			- MBP Num	: Number 				(Tracking)
;			- First Frame	: Of MBP 				(Tracking)	
;			- Last Frame	: Of MBP 				(Tracking)	
;			- Lifetime 	: Of MBP 				(Tracking)
;			- First X/Y 	: Initial Centre point 			(Tracking)
;			- Final X/Y 	: Final centre point 			(Tracking)
;			- Velocity 	: Av Transverse Velocity 		(Tracking)
;			- MAX B 	: Max LOS B field in lifetime 		(Inversions)
;			- MIN B 	: Min LOS B field in lifetime 		(Inversions)
;			- AV B 		: Average B field over lifetime 	(Inversions)
;			- Av Area 	: Average MBP area 			(Tracking)
;			- LOS VEL 	: Average LOS velocity 			(Inversions)
;			- TEMP 		: Average Temperature 			(Inversions)
;
;NOTES:			
;
;AUTHOR:	Peter H. Keys
;
;-

; :: MBP Num :: First Frame :: Last Frame :: Lifetime :: First X/Y :: Final X/Y :: Velocity :: MAX B :: MIN B :: AV B :: Av Area :: LOS VEL :: TEMP :: 

PRO MBP_inversion_properties,file=file
GET_LUN, unit

OPENW,unit, ''+file+''
printf,unit,'       BP number     Initial Frame     Initial x       Initial y      Final Frame       Final x         Final y       Lifetime(s)   Velocity(km/s)       Max B          Min B            Av. B       Max LOS Vel      Min LOS Vel     Av. LOS Vel      Av. Temp       Av. Area (pix)'
CLOSE,unit


; > Inversions files directory.... (these are stored as directories by the MBP number)

inverted_dir = '/data/solarstore3/rlh/27Jul2014SST/final_MBP/evolimages/'

; > Location of the tracking and velocity files... (velocity file most important) 

;tracking_file = FILE_SEARCH('/data/solarstore3/phk/data/27Jul2014/tracking_large_21sig') 		; File I sent Becky (wrong file in the end...)
alt_tracking_file  = FILE_SEARCH('/data/solarstore3/phk/data/27Jul2014/final_tracking_file') 		; More recent file 

; > Load in the Names of the 300 MBPs
;NAMES           DOUBLE    = Array[300]

restore,'/data/solarstore3/rlh/27Jul2014SST/final_MBP/namesof300.sav',/ver
total_MBP_num = N_ELEMENTS(names)

; > Resolutions of the data sets
distance = 50.					
time = 33.

; > Read in the tracking file...
readcol,alt_tracking_file[0],f='D,D,D,D,D,x', btpt_num, frame, x, y, area

; > Open a loop to read in the info for each MBP
FOR n = 0, total_MBP_num -1 DO BEGIN &$
	
	; >> Tracking file properties to be used....
	MBP = names[n] &$
	xcoords = x[WHERE (btpt_num eq MBP)]  &$
	ycoords = y[WHERE (btpt_num eq MBP)]  &$				
	initialx = Xcoords[0] &$
	initialy = Ycoords[0] &$	
	Finalx = Xcoords[N_ELEMENTS(Xcoords)-1] &$
	Finaly = Ycoords[N_ELEMENTS(Ycoords)-1] &$
	frame_num = frame[WHERE (btpt_num eq MBP)] &$
	initial_frame = frame[0] &$
	final_frame = frame_num[N_ELEMENTS(frame_num) - 1]  &$
	lifetime = N_ELEMENTS(frame_num)*time &$
	disp = SQRT((((Finalx - initialx)*distance)^2.) + (((Finaly - initialy)*distance)^2.)) &$
	vel = disp/lifetime &$
	;vel_arr = N_ELEMENTS(xcoords)-1 &$
	;FOR v = 0, N_ELEMENTS(vel_arr)-1 DO BEGIN &$
	;	disp = SQRT((((xcoords[v+1] - xcoords[v])*distance)^2.) + (((ycoords[v+1] - ycoords[v])*distance)^2.)) &$
	;	vel_arr[v] = disp / time &$
	;ENDFOR &$
	;vel = MEAN(vel_arr) &$
	Av_area = MEAN(area[WHERE (btpt_num eq MBP)])  &$
	
	; >> Inversion files properties
	; Need to get the MBP num as a string for directories...
	bp_name = STRING(names[n],FORMAT='(D7.0)') &$
	bp_name = bp_name.compress() &$
	bp_name = bp_name.remove(-1) &$
	restore,inverted_dir+'MBP'+bp_name+'/evolution.sav' &$	; Remove the ,/ver command so it doesn't over populate the screen
	MAX_B = MAX(LOSB,/ABSOLUTE) &$
	MIN_B = MIN(LOSB,/ABSOLUTE) &$
	AV_B = MEAN(LOSB) &$
	MAX_LOSV = MAX(LOSV,/ABSOLUTE) &$
	MIN_LOSV = MIN(LOSV,/ABSOLUTE) &$
	Av_LOSV = MEAN(LOSV) &$
	Av_temp = MEAN(TEMP) &$

; COLUMN NAMES:	
;BP number       
;Initial Frame   
;Initial x       
;Initial y       
;Final Frame       
;Final x        
;Final y      
;Lifetime(s)   
;Velocity(km/s)   
;Max B   
;Min B   
;Av. B   
;Max LOS Vel   
;Min LOS Vel  
;Av. LOS Vel  
;Av. Temp   
;Av. Area (pix)
	
	tab = DBLARR(17,1)  &$					;Creates dimensions of table that the  results are wrote to
	tab[0,*] = MBP	&$	;btpt_num[i]
	tab[1,*] = initial_frame  &$
	tab[2,*] = initialx  &$
	tab[3,*] =  initialy &$
	tab[4,*] = final_frame  &$
	tab[5,*] =  Finalx &$
	tab[6,*] =  Finaly &$
	tab[7,*] =  lifetime &$
	tab[8,*] =  vel &$
	tab[9,*] = MAX_B  &$
	tab[10,*] = MIN_B  &$
	tab[11,*] = AV_B  &$
	tab[12,*] = MAX_LOSV  &$
	tab[13,*] = MIN_LOSV  &$
	tab[14,*] = Av_LOSV  &$
	tab[15,*] = Av_temp  &$
	tab[16,*] = Av_area  &$
	OPENU,unit,''+file+'',/APPEND,WIDTH=650  &$	
		printf,unit,tab  &$						
	CLOSE,unit  &$
ENDFOR

FREE_LUN,unit

END
