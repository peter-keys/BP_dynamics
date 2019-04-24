;+
;
; ROUTINE: 	MBP_detect
;
; PURPOSE: 	TO DETECT MBPS IN TWO DIMENSIONAL IMAGES
; USEAGE:  	DATA, START_IM, END_IM
; INPUT:  	DATA	 - STORAGE DIRECTORY OF THE DATA	
;		START_IM - NUMBER OF FIRST IMAGE IN SEQUENCE
;		END_IM	 - NUMBER OF LAST IMAGE IN SEQUENCE
; OUTPUT:  	DETECT -  BINARY IMAGE OF DETECTED MBPS
; CALLING:	(FITS FORMAT) result = MBP_detect(data='/home/xx/data/',file='/home/xx/results/',start_im=0,end_im=1000)
;		(.SAV FORMAT) result = MBP_detect(data='/home/xx/data/cube.sav',start_im=0,end_im=1000)
;		NB: If keyword FILE is set output shall be in .FITS format to designated directory
;		NB: In .SAV format code shall automatically operate across entire datacube.

FUNCTION MBP_detect, data=data, start_im=start_im, end_im=end_im,imdim=imdim
;----------------------------------------------------------------------------------------------------------------------------
;SECTION ONE: 	LINES 38 TO 57
;PURPOSE: 	To setup some intial parameters, read in data and prepare required output
;		files or array
;PARAMETERS: 	FILES - the original data images in FITS [line **] or .sav format [line **]
;		IM_SIZE - stores the image [x,y] dimensions and number of images in .sav format.
;		EDGE - removes the boundary pixels from the image to permit object growing.
;		EDG - highlights the edge of the image to facilate determination of turning points
;		      in close proximity to image boundary. 				
;		DETECT - binary map of MBPs that have been detected.
;DESCRIPTION:	
;The section checks the data and autoamtically determines its, (FITS or .SAV).  The data is 
;restores to parameter FILES.  If using .SAV format, the user shall be prompted to input the  
;number of the required variable from the list that has been restored. In .SAV format, the    
;code shall execute across the entire datacube automatically. To facilitate the growning of 
;objects close to the image boundary by the REGION_GROW function, parameter EDGE shall remove 
;all pixels that are adjacent to the image boundary, whilst EDG is employed to signal a 
;failure in determining turning points close to the image boundary.
;----------------------------------------------------------------------------------------------------------------------------
IF (N_ELEMENTS(STRSPLIT(''+data+'','.',/EXTRACT))) eq (1) THEN (gofits) = 1. ELSE (gofits)=0.

sl = string(byte([27, 91, 68, 27, 91, 68, 27, 91, 68]))
allofthem = N_ELEMENTS(images)

IF (gofits eq (1.)) THEN BEGIN
	files = FILE_SEARCH(''+data+'*.fits')
	IF (KEYWORD_SET(imdim)) THEN BEGIN
		image = READFITS(files[0],/sil) 
		image = image[imdim[0]:imdim[1],imdim[2]:imdim[3]]
		im_size = SIZE(image)
		;read, paus
	ENDIF ELSE BEGIN
		im_size =SIZE(READFITS(files[0],/sil))
	ENDELSE
ENDIF ELSE BEGIN
	restore,''+data+'',/ver
	sobj= OBJ_NEW('IDL_Savefile',''+data+'')
	READ,inp,PROMPT='PLEASE ENTER NUMBER OF DESIRED INPUT:'
	files = TEMPORARY(SCOPE_VARFETCH(((REVERSE(sobj->Names()))(inp)),level=2)) 
	im_size = SIZE(files)
	imst = 0.
	imend = im_size[3]
ENDELSE	

edge =INTARR(im_size[1],im_size[2])
edg  = MAKE_ARRAY(im_size[1],im_size[2],VALUE=1)  
edge[1:(im_size[1]-2),1:(im_size[2]-2)] =1.
edg[1:(im_size[1]-2),1:(im_size[2]-2)] =0.
detect = INTARR(im_size[1], im_size[2],(end_im-start_im+1))
;----------------------------------------------------------------------------------------------------------------------------
;SECTION TWO: 	LINES 83 TO 94
;PURPOSE: 	To setup some intial parameters, read in data and prepare required output
;		files or array
;PARAMETERS: 	OBJS - binary array that stores all the objects in the image as a value of one
;		SPLT - array that stores the objects that have been separated by the algorithm
;		       so that all objects are investigated
;		TRAN - Transfer array to move the detected objects in a single image to the
;		       output DETECT array [line 350]
;		LIM -  the statistically detemined gradient threshold objects must attain to
;		       to be identified as an MBP. 				
;		IMGST & OBJST - used to store both image and objects in correct orientation.
;DESCRIPTION:	
;This section initiates a for loop to process through each individual image.  If using FITS
;format the images are read in individually or in .SAV format the image is selected from
;the data cube [line 87]  The data image is normalised to the mean [line 88] prior the
;calculation of a sigma value.  This sigma value is employed to set an intensity threshold,
;below pixels shall be considered as inter-granular lanes.  Pixels above this threshold are
;proscibed a value of one in the OBJS array, providing all the objects that shall be
;considered by the code.  As the NORM_IM and OBJS arrays are rotated later in the code,
;parameters STORE and IMGST are set to the correct orientiation of these arrays to facilitate
;the growning and placement of detected MBPs in [SECTION 4] of the code.  Finally the
;procedure MBP_threhold is called to calculate an image specific gradient threshold that
;objects must meet to be considered an MBP
;----------------------------------------------------------------------------------------------------------------------------
for l = (start_im), (end_im) DO BEGIN
	pixel=l
	one_quater = 0.
	half = 0.
	three_quater = 0.
	
	objs = INTARR(im_size[1],im_size[2])
	splt = INTARR(im_size[1],im_size[2])
	tran = INTARR(im_size[1],im_size[2])
	IF (gofits) eq (1.) THEN (im)= READFITS(files[l],/SIL) ELSE (im) = files[*,*,(l)]
	IF (KEYWORD_SET(imdim)) THEN (norm_im) = (im[imdim[0]:imdim[1],imdim[2]:imdim[3]])/(MEAN(im[imdim[0]:imdim[1],imdim[2]:imdim[3]])) ELSE (norm_im) = (im)/(MEAN(im))
	;read, paus
	sigma = SQRT((MOMENT(norm_im[WHERE(norm_im gt (0.) AND (norm_im) lt (10.))]))(1))
	objs[WHERE(norm_im ge ((MEAN(norm_im))+(0.5*sigma)))] = 1.
	objs = (objs)*(edge)
	imgst = norm_im
	store = objs
	lim = MBP_THRESHOLD(norm_im,objs,edg,splt,im_size)
;----------------------------------------------------------------------------------------------------------------------------
;SECTION THREE: LINES 137 TO 224
;PURPOSE: 	To examine each object and determine the intensity gradient across them in 4
;		separate directions [i.e. 0,45,90,135 degrees].  If objects attain the
;		requried threshold in all four directions they pass into [SECTION FOUR] 
;PARAMETERS: 	GRAD - an interger array used to signal if an object has been successfully in
;		attaining the required gradient threshold [SUCCESS->VALUE=1 FAIL ->VALUE=0]
;		OBJ - array that stores the individual object underconsideration
;		INTIM - array holding the intensity values of the considered object.
;		SPLGRP - Binary array that stores an object after it has been separated by the
;			 function MBP_SPLIT.
;		DIST -  four element array holding the distance between the turning points of
;			the intensity profiles				
;		CONST - Array holding the extension of the slice [1st Element], rotation angle
;			[2nd ELement] and image rotation alert [3rd Element].
;		TPT - Two element array storing the values of the detemined turning points
;		INTER - an array employed to highlight internal intensity structure of the
;			considered object, which is ignored when searching for turning points.
;		LINE - the slice array
;DESCRIPTION:	
;Section three enters a WHILE loop to examine objects individually until every object has been
;investigated.  The algorithm grows each indiviudal obeject [lines 140:141] and removes it from
;further consideration [line 142]  Prior to intensity profiling, if the object has a large 
;internal intensity range it shall be split apart [lines 144 to 154]  This is to allow all
;objects to be examined individually. The function MBP_split performs this separation [line
;148].  Once separated the object, if it conforms to certain dimensions [line 156], shall be 
;profiled in four directions, (0,45,90,135 degrees).  The objects Center of Gravity [CofG] is 
;determined [lines 161:164].  A slice is placed and rotated through the objects CofG 
;[lines 171:178] and any internal intensity structure is identified in the slice [lines 180:
;184]  If the slice reaches the edge of the image without two turning points being
;detected, then the process is quit and the object is ignored.  Major tunring points are 
;indentified [lines 186:191] and stored in array TPT.  If two turing points are not found the 
;slice/line is automatically extended and the profile retaken [line 193]. Once two turning
;points are found the maximum gradient, on either side of the maximum is taken.  If either
;fails to reach the required threshold the process is stopped [line 201] and the object is not
;considered an MBP.  If the gradients are above the threshold then the GRAD array is given a
;value 1 for each of the four directions and the distance between the turning points is also
;record[lines 203:206].  Once past 90 degrees, IDL reads in the slice in the wrong
; thus, to correct this, the image and object are rotated by 90 degrees and the slice rotated a
;further 45 degrees to emulate 135 degree rotation.  Finally, if 2 of the recorded distances
;are large, then the object is rejected as an MBP. Lines 244:248 need to be augmented for the 
;spatial resolution of the data being analysed.
;----------------------------------------------------------------------------------------------------------------------------
	WHILE (MAX(objs)) gt (0.) DO BEGIN
		grad = INTARR(4)
		obj = INTARR(im_size[1],im_size[2])
		ROI = REGION_GROW(objs,((WHERE(objs eq (1)))(0)))
		obj[ROI]=objs[ROI]
		objs = (objs)-(obj)
			
		IF (N_ELEMENTS(WHERE(obj gt (0.)))) ge (4.) AND (MAX(splt+obj)) eq (1.) THEN BEGIN &$
			splt=(splt)+(obj) &$
			intim = (obj)*(norm_im) &$
			siga = SQRT((MOMENT(intim[WHERE((intim) gt (0.) AND (intim) lt (10.))]))(1))
			IF (((MAX(intim))-(MIN(intim[WHERE(intim gt 0.)])))/siga) gt (3.) THEN BEGIN
				splgrp = MBP_SPLIT(intim,im_size,siga) &$
				store  = ((store-obj)+splgrp) &$
				objs   = (objs)+(splgrp) &$
				obj    = (obj)-(obj) &$
				reg    = REGION_GROW(splgrp,((WHERE(splgrp eq (1)))(0))) 
				obj[reg]=splgrp[reg] &$
				objs = (objs)-(obj)
			ENDIF &$
		ENDIF

		IF (N_ELEMENTS(WHERE(obj eq (1.)))) lt (300) AND (N_ELEMENTS(WHERE(obj eq (1.)))) ge (4.)THEN BEGIN &$
			dist = INTARR(4)
			const = INTARR(3)
			const[0] = 20
			
			s = SIZE(obj, /Dimensions)       
 	 		totalMass = TOTAL(obj) 
   			x = ROUND(TOTAL(TOTAL(obj,2)*INDGEN(s[0]))/totalMass)    
   			y = ROUND(TOTAL(TOTAL(obj,1)*INDGEN(s[1]))/totalMass)
			
			WHILE (const[1]) lt (4) DO BEGIN &$
				tpt = MAKE_ARRAY(2,VALUE=-1) &$
				inter = INTARR(im_size[1], im_size[2]) &$
				line = INTARR(im_size[1],im_size[2]) &$ 
				IF (y-const[0]) lt (0.) OR (y+const[0]) gt (im_size[2]-1) THEN BREAK  &$
				line[x,(y-const[0]):(y+const[0])] = 1. &$
				
				IF (const[1]) lt (3) THEN (rot) = ROT(line,(const[1]*45),1,x,y,MISSING = 0.,/PIVOT) ELSE (rot) = ROT(line,(45),1.0,x,y,MISSING=0.,/PIVOT) &$
				IF (MAX(rot+edg)) eq (2) THEN BREAK &$
				slice = (norm_im)*(rot) &$
				arr = slice[WHERE(slice gt 0)] &$
				dx = MBP_GRADIENT(SMOOTH(arr,2)) &$
				
				IF (MAX(obj+rot)) gt (1.) THEN (inter[WHERE((obj+rot) eq (2))]) = (1.) &$
				inter = (inter)+(rot) &$
				ins = inter[WHERE(inter gt (0.))] &$
				IF (MAX(obj+rot)) gt (1.) THEN (maxpt) = WHERE(arr eq (MAX(arr[WHERE(ins eq 2)]))) ELSE (maxpt) = WHERE(arr eq (MAX(arr[((N_ELEMENTS(arr)/2)-5):((N_ELEMENTS(arr)/2)+5)])))
				IF ((WHERE(ins gt (1.)))(0)) gt (-1) THEN (ins[(MIN(WHERE(ins eq (2)))):(MAX(WHERE(ins eq (2))))]) = 2. &$
				
				for j = (3.), (N_ELEMENTS(dx)-4) DO BEGIN &$
					IF (ins[j]) lt (2.) AND (dx[j-1]) lt (0.) AND (dx[j]) gt (0.) AND  (TOTAL(dx[(j-3):(j-1)])) le (0.) AND (TOTAL(dx[(j+1):(j+3)])) ge (0.) THEN BEGIN &$
						IF (j) lt (maxpt[0]) AND (j) gt (tpt[0]) THEN (tpt[0]) = (j) &$
						IF (j) gt (maxpt[0]) AND (tpt[1]) eq (-1) THEN (tpt[1]) = (j) &$
						IF (tpt[0]) gt (-1) AND (tpt[1]) gt (-1) THEN BREAK
					ENDIF &$
				ENDFOR &$
				
				IF (tpt[0]) eq (-1) OR (tpt[1]) eq (-1) THEN (const[0]) = (const[0]) + 20. &$
				
				IF (tpt[0]) gt (-1) AND (tpt[1]) gt (-1) THEN BEGIN &$
					leftdx  = MBP_GRADIENT(arr[tpt[0]:maxpt[0]]) &$
					rightdx = MBP_GRADIENT(arr[maxpt[0]:tpt[1]]) &$
					IF (MAX(leftdx, /ABSOLUTE)) lt (0.) THEN (maxleft) =(MAX(leftdx, /ABSOLUTE)) *(-1) ELSE (maxleft) =MAX(leftdx,/ABSOLUTE)  &$
					IF (MAX(rightdx,/ABSOLUTE)) lt (0.) THEN (maxright)=(MAX(rightdx,/ABSOLUTE))*(-1) ELSE (maxright)=MAX(rightdx,/ABSOLUTE)  &$

					IF (maxleft) lt (lim) OR (maxright) lt (lim) THEN BREAK &$
					
					grad[const[1]] = 1.  &$
					dist[const[1]] = (tpt[1]) - (tpt[0]) &$
					diffl = (maxleft) -(lim) &$
					diffr = (maxright)-(lim) &$
					
					const[0] = 20. &$
					const[1] = (const[1])+(1) &$

					IF (const[1]) eq (3.) THEN BEGIN &$
						const[2]=1. &$
						obj = ROTATE(obj,1) &$
						edg = ROTATE(edg,1) &$
						norm_im = ROTATE(norm_im,1) &$
						im_size = SIZE(norm_im) &$
						s = SIZE(obj, /Dimensions) &$       
 	 					totalMass = TOTAL(obj)  &$
   						x = ROUND(TOTAL(TOTAL(obj,2)*INDGEN(s[0]))/totalMass)  &$   
   						y = ROUND(TOTAL(TOTAL(obj,1)*INDGEN(s[1]))/totalMass) &$
					ENDIF &$
				ENDIF &$
			ENDWHILE &$
			
			IF ((WHERE(dist ge 164))(0)) gt (-1.) THEN (grad[*]) = 0.    ;HiFi GREGOR approx 0.0378"/pix
			;IF ((WHERE(dist ge 180))(0)) gt (-1.) THEN (grad[*]) = 0.    ;Sim data 0.0345"/pix
			;IF ((WHERE(dist ge 145))(0)) gt (-1.) THEN (grad[*]) = 0.	;BIFROST simulations (31km/pix)
			;IF ((WHERE(dist ge 105))(0)) gt (-1.) THEN (grad[*]) = 0.	;TEST OF THE SST PARAMS
			;IF ((WHERE(dist ge 104.9))(0)) gt (-1.) THEN (grad[*]) = 0.    ;AR SST data 0.0592"/pix			
			;IF ((WHERE(dist ge 57))(0)) gt (-1.) THEN (grad[*]) = 0.    ;Hinode data 57 0.108"/pix		
			;IF ((WHERE(dist ge 90))(0)) gt (-1.) THEN (grad[*]) = 0.    ;Observed ROSA data 0.069"/pix
			;IF ((WHERE(dist ge 75))(0)) gt (-1.) THEN (grad[*]) = 0.    ;IBIS ROSA data 0.0828"/pix
;----------------------------------------------------------------------------------------------------------------------------
;SECTION FOUR: 	LINES 244 TO 342
;PURPOSE: 	To grow the objects that have achieved the required gradient threshold in all
;		directions from [SECTION 3]
;PARAMETERS: 	LIMIT_CACHE - stores the 36 individual turning point intensity levels at each
;			      of the 36 considered angles 
;		GROWN - binary array to store the grown MBP
;DESCRIPTION:	
;[SECTION 4] begins by correcting for image and object rotation [lines 249:252] Section four
;employs the same basic method as [SECTION 3], to search for turning points in the intensity 
;profiles of the MBPs [lines 262:310].  However, [SECTION 4] rotates the line through 360 
;degrees in 5 degree steps and looks to store the highest turning point intensity at each angle
;[line 294]  Finally the highest turning point intensity is employed as a limit, above which
; conjoining pixels are considered part of the MBP [line 312].  The grown MBP is placed into
;the binary array GROWN [line 315:317] and the code checks that no other objects has been grown
;[lines 319:326].  If it is found that another object has been grown,then the threshold is
;varied to restrict the growth to a single object[line 320]  Finally, if the final gradient
;failed in [SECTION 3] the code corrects the image rotation [lines 333 to 338]. Lines 361:365 
;need to be altered depending on the data. This is effectively the cut-off to stop granules being detected.
;----------------------------------------------------------------------------------------------------------------------------
			IF (N_ELEMENTS(WHERE(grad eq 1))) eq (4.) THEN BEGIN
				const[0] = 20
				const[1:2]=0
				limit_cache = FLTARR(36)
								
				obj = ROTATE(obj,3)
				edg = ROTATE(edg,3)
				norm_im = ROTATE(norm_im,3)
				im_size = SIZE(norm_im)
				
				objst = obj
				store =(store)-(objst)
				
				s = SIZE(obj, /Dimensions)       
 	 			totalMass = TOTAL(obj) 
   				x = ROUND(TOTAL(TOTAL(obj,2)*INDGEN(s[0]))/totalMass)    
   				y = ROUND(TOTAL(TOTAL(obj,1)*INDGEN(s[1]))/totalMass)
				xst = x
				yst = y
				
				WHILE (const[1]) lt (36) DO BEGIN &$
					tpt  = MAKE_ARRAY(2,VALUE=-1) &$
					inter = INTARR(im_size[1], im_size[2]) &$
					line = INTARR(im_size[1],im_size[2]) &$ 
					IF (y-const[0]) lt (0.) THEN (ylow) = 0. ELSE (ylow) = (y-const[0]) &$ 
					IF (y+const[0]) gt (im_size[2]-1) THEN (yhig) = (im_size[2]-1) ELSE (yhig) = (y+const[0]) &$
					line[x,(ylow):(yhig)] = 1. &$
					
					IF (const[1]) lt (19) THEN (rot) = ROT(line,(const[1]*5),1,x,y,missing = 0.,/PIVOT) ELSE (rot) = ROT(line,((const[1]*5)-90),1.0,x,y,MISSING=0.,/PIVOT) &$
					slice = (norm_im)*(rot) &$
					arr = slice[WHERE(slice gt 0)] &$
					dx = MBP_GRADIENT(SMOOTH(arr,2)) &$
					
					IF (MAX(obj+rot)) gt (1.) THEN (inter[WHERE((obj+rot) eq (2))]) = (1.) &$
					inter = (inter)+(rot) &$
					ins = inter[WHERE(inter gt (0.))] &$
					IF (MAX(obj+rot)) gt (1.) THEN (maxpt) = WHERE(arr eq (MAX(arr[WHERE(ins eq 2)]))) ELSE (maxpt) = WHERE(arr eq (MAX(arr[((N_ELEMENTS(arr)/2)-5):((N_ELEMENTS(arr)/2)+5)])))
					IF ((WHERE(ins gt (1.)))(0)) gt (-1) THEN (ins[(MIN(WHERE(ins eq (2)))):(MAX(WHERE(ins eq (2))))]) = 2. &$

					for j = (3.), (N_ELEMENTS(dx)-4) DO BEGIN &$
						IF (ins[j]) lt (2.) AND (dx[j-1]) lt (0.) AND (dx[j]) gt (0.) AND  (TOTAL(dx[(j-3):(j-1)])) le (0.) AND (TOTAL(dx[(j+1):(j+3)])) ge (0.) THEN BEGIN &$
							IF (j) lt (maxpt[0]) AND (j) gt (tpt[0]) THEN (tpt[0]) = (j) &$
							IF (j) gt (maxpt[0]) AND (tpt[1]) eq (-1) THEN (tpt[1]) = (j) &$
							IF (tpt[0]) gt (-1) AND (tpt[1]) gt (-1) THEN BREAK
						ENDIF &$
					ENDFOR &$
					
					IF ((WHERE((rot+edg) eq (2)))(0)) eq ((WHERE(rot eq 1))(0)) THEN (tpt[0]) = 0.  &$
					IF ((WHERE((rot+edg) eq (2)))((N_ELEMENTS(WHERE((rot+edg) eq (2))))-1)) eq ((WHERE(rot eq 1))((N_ELEMENTS(WHERE(rot eq 1)))-1)) THEN (tpt[1]) = 0.  &$
					IF (WHERE((line+edg) eq (2))) eq ((WHERE(line eq 1))(0)) AND (tpt[0]) eq (-1.) THEN (tpt[0]) = 0.  &$
					IF (WHERE((line+edg) eq (2))) eq ((WHERE(line eq 1))((N_ELEMENTS(WHERE(line eq 1)))-1)) AND (tpt[1]) eq (-1.) THEN (tpt[1]) = 0.  &$

				 	IF (tpt[0]) eq (-1) OR (tpt[1]) eq (-1) THEN (const[0]) = (const[0]) + 20 &$

					IF (tpt[0]) gt (-1) AND (tpt[1]) gt (-1) THEN BEGIN &$
						IF (arr[tpt[0]]) gt (arr[tpt[1]]) AND (tpt[0]) gt (0.) THEN (limit_cache[const[1]]) = (arr[tpt[0]]) ELSE (limit_cache[const[1]])=(arr[tpt[1]]) &$
						const[0] = 20 &$
						const[1] = (const[1])+(1) &$

						IF (const[1]) eq (19.) THEN BEGIN &$
							const[2]=1. &$
							obj = ROTATE(obj,1) &$
							edg = ROTATE(edg,1) &$
							norm_im = ROTATE(norm_im,1) &$
							im_size = SIZE(norm_im) &$
							s = SIZE(obj, /Dimensions) &$       
 	 						totalMass = TOTAL(obj)  &$
   							x = ROUND(TOTAL(TOTAL(obj,2)*INDGEN(s[0]))/totalMass)  &$   
   							y = ROUND(TOTAL(TOTAL(obj,1)*INDGEN(s[1]))/totalMass) &$
						ENDIF &$
					ENDIF &$
				ENDWHILE &$
			
				RofI = REGION_GROW(imgst,(WHERE(objst eq (1))),THRESHOLD=[MAX(limit_cache),MAX(imgst)])
				
				IF (RofI[0]) gt (-1.) THEN BEGIN
					;IF (im_size[2]) eq (704) THEN read, paus
					grown = INTARR(im_size[2],im_size[1])
					grown[RofI]=imgst[RofI]
					grown[WHERE(grown gt (0.))]=1.
					vary=0.001
					
					WHILE (MAX(grown+store)) gt (1.) DO BEGIN &$ 
						RofI = REGION_GROW(imgst, RofI, THRESHOLD=[(MAX(limit_cache)+vary),(MAX(imgst*objst))]) &$
						IF (RofI[0]) eq (-1.) THEN BREAK &$
						grown = INTARR(im_size[2],im_size[1]) &$
						grown[RofI]=imgst[RofI] &$
						grown[WHERE(grown gt (0.))]=1. &$
						grown = MBP_close(grown,xst,yst) &$ 
						vary = vary+0.001 &$
					ENDWHILE
					
					IF (AVG(diffl)) lt (0.1) OR (AVG(diffr)) lt (0.1) AND (N_ELEMENTS(WHERE(grown eq 1.))) gt (50) THEN (grown[*,*]) = 0.
					;CHANGE THE VALUE BELOW BASED ON THE RES: (eg. ROSA & SIM): (0.069"/0.0345")^2 * 200.
					IF (N_ELEMENTS(WHERE(grown eq 1.))) le (666.) AND (N_ELEMENTS(WHERE(grown eq 1.))) ge (4)  THEN (tran) = (tran)+(grown)  ;HiFI GREGOR					
					;IF (N_ELEMENTS(WHERE(grown eq 1.))) le (90.) AND (N_ELEMENTS(WHERE(grown eq 1.))) ge (4)  THEN (tran) = (tran)+(grown)  ;TESTING FOR SST WB					
					;IF (N_ELEMENTS(WHERE(grown eq 1.))) le (272.) AND (N_ELEMENTS(WHERE(grown eq 1.))) ge (4)  THEN (tran) = (tran)+(grown)  ;Aaron SST data					
					;IF (N_ELEMENTS(WHERE(grown eq 1.))) le (80.) AND (N_ELEMENTS(WHERE(grown eq 1.))) ge (4)  THEN (tran) = (tran)+(grown)  ;Hinode data					
					;IF (N_ELEMENTS(WHERE(grown eq 1.))) le (800.) AND (N_ELEMENTS(WHERE(grown eq 1.))) ge (4)  THEN (tran) = (tran)+(grown)  ;MuRAM Sim data
					;IF (N_ELEMENTS(WHERE(grown eq 1.))) le (521.) AND (N_ELEMENTS(WHERE(grown eq 1.))) ge (4)  THEN (tran) = (tran)+(grown)  ;BIFROST Simulations					
					;IF (N_ELEMENTS(WHERE(grown eq 1.))) le (200.) AND (N_ELEMENTS(WHERE(grown eq 1.))) ge (4)  THEN (tran) = (tran)+(grown)  ;Observed ROSA data
					;IF (N_ELEMENTS(WHERE(grown eq 1.))) le (139.) AND (N_ELEMENTS(WHERE(grown eq 1.))) ge (4)  THEN (tran) = (tran)+(grown)  ;IBIS data

					IF (MAX(tran)) gt (1) THEN (tran[WHERE(tran gt 1)]) = 1.
				ENDIF
				store = (store)+(objst)
			ENDIF

			IF (const[2]) eq (1.) THEN BEGIN &$
				obj = ROTATE(obj,3) &$
				edg = ROTATE(edg,3) &$
				norm_im = ROTATE(norm_im,3) &$
				im_size = SIZE(norm_im) &$
				const[2] = 0.
			ENDIF
		
		IF (N_ELEMENTS(WHERE(objs[*,0:(im_size[2]*0.25)] eq 0.))) eq (N_ELEMENTS(objs[*,0:(im_size[2]*0.25)])) AND (one_quater) eq (0.) THEN print, '1/4 of  frame    ', sl+string(pixel, format="(i5)"), ' processed'  
		IF (N_ELEMENTS(WHERE(objs[*,0:(im_size[2]*0.5)] eq 0.)))  eq (N_ELEMENTS(objs[*,0:(im_size[2]*0.5)])) AND (half) eq (0.) THEN print, '1/2 of  frame    ', sl+string(pixel, format="(i5)"), ' processed'  
		IF (N_ELEMENTS(WHERE(objs[*,0:(im_size[2]*0.75)] eq 0.))) eq (N_ELEMENTS(objs[*,0:(im_size[2]*0.75)])) AND (three_quater) eq (0.) THEN print, '3/4 of  frame    ', sl+string(pixel, format="(i5)"), ' processed'  
		IF (N_ELEMENTS(WHERE(objs[*,0:(im_size[2]*0.25)] eq 0.))) eq (N_ELEMENTS(objs[*,0:(im_size[2]*0.25)])) THEN (one_quater) = 1.  
		IF (N_ELEMENTS(WHERE(objs[*,0:(im_size[2]*0.5)] eq 0.)))  eq (N_ELEMENTS(objs[*,0:(im_size[2]*0.5)])) THEN (half) = 1.  
		IF (N_ELEMENTS(WHERE(objs[*,0:(im_size[2]*0.75)] eq 0.))) eq (N_ELEMENTS(objs[*,0:(im_size[2]*0.75)])) THEN (three_quater) = 1.  

		ENDIF
		
	ENDWHILE
;----------------------------------------------------------------------------------------------------------------------------
;SECTION FIVE: 	LINES  TO 
;PURPOSE: 	Stores the binary detection image as either a .FITS file or as part of an array
;PARAMETERS: 	NONE
;DESCRIPTION:	
;Saves the detection image, a binary map, in either .FITS format or as part of an output array.
;----------------------------------------------------------------------------------------------------------------------------
	detect[*,*,(l-start_im)] = tran
ENDFOR

RETURN, detect

end
;----------------------------------------------------------------------------------------------------------------------------
;----------------------------------------------------------------------------------------------------------------------------
