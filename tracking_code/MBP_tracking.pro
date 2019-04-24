;@/home/pjc/idl/revised_vers/version_two/two_five/dotr.bat 
;tracking=MBP_tracking(data='/data/rosa3/oldrosa1/phk/data/28May2009/gband/gband_quiet_2sec.sav',array=detect,file='test.dat',imst=0,imend=20,cad=2.)
;+
;
; ROUTINE:  MBP_tracking
;
; PURPOSE:  TO TRACK MBPS THROUGHOUT THEIR LIFETIME, STABILISING ANY WHICH ARE MISSING
;	    FROM THE DETECTION FRAMES
; USEAGE:   DATA, ARR, TABLE, IMST, IMEND, CAD
; INPUT:    DATA - DIRECTORY OF ORIGINAL DATA.  IF DATA IN .SAV FORMAT INCLUDE THE NAME
;	    OF THE DATA CUBE IN THE CALLING PROCESS
;	    ARR - THE DETECTION ARRAY FROM THE SEARCH ALGORITHM
;	    FILE - DIRECTORY AND NAME OF FINAL OUTPUT TABLE FILE
;	    IMST - THE FIRST IMAGE FROM WHICH TO TRACK
;	    IMEND - FINAL IMAGE OF THE TIME SERIES
;	    CAD - CADENCE OF THE DATA
; OUTPUT:   STABLE - ARRAY WHERE THE LONG LIVED FEATURES ARE PRESENT THROUGHOUT THEIR
;	    LIVETIME
; CALLING:  (FITS FORMAT) RESULT = MBP_tracking(data='/home/xx/data/',array=array,file='/home/xx/result.dat',imst=0,imend=1000,cad=2)
;	    (.SAV FORMAT) RESULT = MBP_tracking(data='/home/xx/data/file.sav',array=array,file='/home/xx/result.dat',imst=0,imend=1000,cad=2)			   	   
;	    NB: In .SAV format the code shall operate across the entire datacube. 
;	    NB: If keyword OUTPUT set, the resulting images shall be written to the designated
;		file in .FIT
;AUTHOR:    Version 2.3 - Philip. J. Crockett, QUB, 30 NOV 2010
;   	    	    	      (Email: pcrockett02@qub.ac.uk)
;		Adapted by P. H. Keys, QUB, 2011 - present...

;AMENDMENT HISTORY
;21/01/2011 - (v2.1)-Higher initial intensity threshold employed to set a seed region.  Now potential seed regions must
;		     overlap original object by 1/2 the number of pixels prior to being separated.  This stops the inclusion
;		     of single pixel seeds.
;21/01/2011 - (v2.1)-Adjustable growing threshold to ensure every chance is given to grow a missed MBP.
;21/01/2011 - (v2.1)-Grown array is now set to zero during varying stage to prevent erronous inclusion of small MBPs
;24/01/2011 - (v2.2)-Addition of MBP_seed program to store multiple seed regions if object is observed to separate
;17/02/2011 - (v2.3)-The mergin of two objects found with an associated object is now recognised
;17/02/2011 - (v2.3)-THe merging of objects is now noted with the addition of a 7th column to table using MBP_merger v2.1


FUNCTION MBP_tracking, data=data, array=array, file=file,imst=imst,imend=imend,cad=cad
;----------------------------------------------------------------------------------------------------------------------------
;SECTION ONE: 	LINES 46 TO 72
;PURPOSE: 	To setup some intial parameters, read in data and prepare required output files.
;PARAMETERS: 	FILES - the original data images in FITS [line **] or .sav format [line **]
;		IM_SIZE - stores the image [x,y] dimensions and number of images in .sav format.
;		IMGS - number of images that constitute 30secs of data depending on cadence. Employed
;		       to search successive frames if MBP is found to be missing
;		EDG - highlights the edge of the image to facilate determination of turning points
;		      in close proximity to image boundary. 				
;		STABLE - binary map of MBPs that have been detected or stabilised.
;		LABEL - number used to identify individual MBPs
;DESCRIPTION:	
;The section checks if output files exist and deletes them, prior to opening new files.
;(IMPORTANT NOTE:CHANGE FILENAME '/home/pjc/table.dat' [line 51&56] TO SUIT INDIVIDUAL SYSTEM).  The
;code autoamtically determines the format of the data,(FITS or .SAV), and restores the required data
;to parameter FILES.  If using .SAV format, the user shall be prompted to input the number of the
;required variable from the list that has been restored. In .SAV format, the code shall execute
;across the entire datacube automatically. The defined IMGS parameter shall be employed to hunt
;missing MBPs across 30secs of data [section 5].  Finally, to facilitate the growning of objects
;close to the image boundary by the REGION_GROW function, parameter EDGE shall remove all pixels that are
;adjacent to the image boundary, whilst EDG is employed to signal a failure in determining
;turning points close to the image boundary.  LABEL is emplyed to identify each detected MBP		
;----------------------------------------------------------------------------------------------------------------------------
sl = string(byte([27, 91, 68, 27, 91, 68, 27, 91, 68]))
IF (FILE_TEST(''+file+'')) THEN FILE_DELETE,''+file+'' 
IF (FILE_TEST('table.dat')) THEN FILE_DELETE,'table.dat'
IF (N_ELEMENTS(STRSPLIT(''+data+'','.',/EXTRACT))) eq (1) THEN (gofits) = 1. ELSE (gofits)=0.

OPENW,1,''+file+''
CLOSE,1
OPENW,1,'table.dat'
CLOSE,1

IF (gofits eq (1.)) THEN BEGIN
	files = FILE_SEARCH(''+data+'*.fits')
	im_size = SIZE(READFITS(files[0],/SIL))
	allofthem=N_ELEMENTS(files)-1
ENDIF ELSE BEGIN
	restore,''+data+'',/ver
	sobj= OBJ_NEW('IDL_Savefile',''+data+'')
	READ,inp,PROMPT='PLEASE ENTER NUMBER OF DESIRED INPUT:'
	files = TEMPORARY(SCOPE_VARFETCH(((REVERSE(sobj->Names()))(inp)),level=2)) 
	im_size = SIZE(files)
	allofthem=im_size[3]
ENDELSE

;imgs = ((ROUND(192./cad))+(1))				;Augmented for the reduced cadence of HINODE
;imgs = ((ROUND(117./cad))+(1))				;Augmented for the reduced cadence of IBIS
;imgs = ((ROUND(30./cad))+(1))				;Looks 30s to see if it is still there
imgs = ((ROUND(606./cad))+(1))				;Augmented for SST WB data (303s cadence)
edg  = MAKE_ARRAY(im_size[1],im_size[2],VALUE=1)  
edg[1:(im_size[1]-2),1:(im_size[2]-2)] =0.
stable = INTARR(im_size[1],im_size[2],((imend-imst)+1))
label = -1
;----------------------------------------------------------------------------------------------------------------------------
;SECTION TWO: 	LINES 103 TO 157
;PURPOSE: 	Investigates each individual MBP that has been detected, records their
;		parameters and prints them to file. 
;PARAMETERS: 	NEXT - binary array set up to store the associated MBP in the next frame.
;		FRAME - binary array set up to store all MBPs that have an associated MBP in the
;			following detection frame.
;		OBJS - the detected objects that are being investigated.
;		OBJ - the current MBP under investigate.
;		FOLLOW - the successive/following detection frame.				
;		TABLE - storage of information from the previous two frames.
;		INFO -  an array storing a unique BTPT number, frame number, [x,y] position and
;			area of the objects.  The final column indicates what frame a missing
;			MBP is found in.
;		ASSOC - array storing the associated MBP		
;DESCRIPTION:	
;This section reads in the succeding data image, either from .FITS format or .SAV format.  The
;image is normalised to the mean and the sigma value determined.  Initally this section records 
;all of the information concerning the MBPs detected in the first image, OBJS(line[113]).  
;The following detection image is arrange in parameter FOLLOW (line[***]) to allow object association. 
;After this the code enters a WHILE loop until all objects have been investigated.  Each object is 
;grown (lines[137:138]) prior to the calculation of the center of gravity and the placement of the 
;information into INFO (lines[149:154]). After the intial image, this section writes all the information 
;from the previous two frames to file, throught the parameter''+FILE+''(lines[117:13Û]). The function 
;MBP_tabulating.pro, takes the information from ''+FILE+'' and searches for associated MBPs in the 
;following image.  Any new MBPs that appear in the image shall be inputted to OBJS (line[132]) and 
;searched for in the following sections.
;----------------------------------------------------------------------------------------------------------------------------
for l = (imst), (imend-1) DO BEGIN
	pixel=l
	next = INTARR(im_size[1],im_size[2])
	
	IF (gofits) eq (1.) THEN (imnx)= READFITS(files[l+1],/SIL) ELSE (imnx) = files[*,*,(l+1)]
	norm = (imnx)/(MEAN(imnx))
	sigi = SQRT((MOMENT(norm[WHERE((norm gt (0.)) AND (norm lt (10.)))]))(1))

	IF (l) eq (imst) THEN BEGIN &$
		frame =INTARR(im_size[1],im_size[2]) &$
		objs  =array[*,*,(l)] &$
	ENDIF ELSE BEGIN &$
		;readcol,'table.dat',f='D,D,D,D,D',bptno,frnum,x,y,area,exist,merged,/sil &$	        ;MERGING
		readcol,'table.dat',f='D,D,D,D,D',bptno,frnum,x,y,area,exist,/sil &$	       ;NON MERGING
		;tab = INTARR(7,N_ELEMENTS(bptno)) &$							;MERGING
		tab = INTARR(6,N_ELEMENTS(bptno)) &$							;NON MERGING
		tab[0,*] = bptno &$
		tab[1,*] = frnum&$
		tab[2,*] = x &$
		tab[3,*] = y &$
		tab[4,*] = area &$
		tab[5,*] = exist &$
		;tab[6,*] = merged &$
		index = SORT(tab[1,WHERE(tab[1,*] le (l))]) &$ 
		;FOR j=0,6 DO (tab[j,*]) = tab[j,(index)] &$					;MERGING
		FOR j=0,5 DO (tab[j,*]) = tab[j,(index)] &$					;NON MERGING
		OPENU,1,''+file+'',/APPEND &$
		printf,1, tab &$
		CLOSE,1 &$
		OPENW,1,'table.dat'
		CLOSE,1
		btpts =(array[*,*,(l)])+(follow) &$ 
		IF (MAX(btpts gt (1.))) THEN (btpts[WHERE(btpts gt (1.))]) = 1. &$
		strut = MBP_tabulating(''+file+'',norm,btpts,array,next,edg,imgs,l,imend,sigi,im_size) &$
		objs = strut.objs
		next = strut.next
		frame = stable[*,*,(l-imst)] &$
	ENDELSE

	IF (l) eq (imst) THEN (follow) = array[*,*,(l+1)] ELSE (follow) = (array[*,*,(l+1)])+(strut.follow)
	IF (MAX(follow)) gt (1.) THEN (follow[WHERE(follow gt 0)]) = 1.	
	
	IF (FILE_SIZE(''+file+'')) eq (0.) THEN BEGIN
		prev = -1.
		xprev = -1.
		yprev = -1.
	ENDIF
	
	IF (FILE_SIZE(''+file+'')) gt (0.) THEN readcol,''+file+'',f='I,X,I,I,X,X,X',prev,xprev,yprev,/SIL

	WHILE (MAX(objs)) gt (0.) DO BEGIN &$
		
		obj     = INTARR(im_size[1],im_size[2]) 
		info    = MAKE_ARRAY(7,2,/INT,VALUE=-1)
		assoc   = INTARR(im_size[1],im_size[2])
		grown   = FLTARR(im_size[1],im_size[2])
		growing = 0.
		ROI = REGION_GROW(objs, ((WHERE(objs eq (1.)))(0))) 
		obj[ROI] = objs[ROI] 
		objs = (objs)-(obj)
				
		IF (N_ELEMENTS(WHERE(obj gt (0.)))) ge (4.) THEN BEGIN &$
			label = (label)+(1) &$		
			frame = (frame)+(obj) &$
			IF (MAX(frame)) gt (1.) THEN (frame[WHERE(frame gt 1)])=1 &$
			sobj = SIZE(obj, /Dimensions) &$	
    			mobj = TOTAL(obj) &$	
    			xobj = ROUND(TOTAL( TOTAL (obj,2) * Indgen(sobj[0]) ) / mobj) &$	
    			yobj = ROUND(TOTAL( TOTAL (obj,1) * Indgen(sobj[1]) ) / mobj) &$	
			info[0,0] = label &$	
			info[1,0] = (l) &$	
			info[2,0] = (xobj) &$	
			info[3,0] = (yobj) &$	
			info[4,0] = N_ELEMENTS(WHERE(obj gt (0.))) &$	
			info[5,0] = (l) &$
;----------------------------------------------------------------------------------------------------------------------------
;SECTION THREE: LINES 167 TO 170
;PURPOSE: 	Search for associated MBP in the following image
;PARAMETERS: 	NONE
;DESCRIPTION:
;This section takes each object from the frame under investigation (OBJ) and searchs for any objects
;in the following frame (FOLLOW) that overlap the object
;----------------------------------------------------------------------------------------------------------------------------
			IF (N_ELEMENTS(WHERE((obj+follow) gt (1.)))) ge (2.) THEN BEGIN &$
				RofI  = REGION_GROW(follow, WHERE((obj+follow) gt (1.))) &$
				assoc[RofI] = follow[RofI] &$
				assoc = MBP_close(assoc,xobj,yobj)&$
				IF (N_ELEMENTS(WHERE(assoc gt 0))) lt (4) THEN (assoc)=INTARR(im_size[1],im_size[2])
			ENDIF 
;----------------------------------------------------------------------------------------------------------------------------
;SECTION FOUR: 	LINES 179 TO 191
;PURPOSE: 	If an associated MBP is found this section records the details about it
;PARAMETERS: 	NONE
;DESCRIPTION:
;If the ASSOC MBP is found the section records all the information concerning it and place it
;into the info array
;----------------------------------------------------------------------------------------------------------------------------
			IF (MAX(assoc)) gt (0.) THEN BEGIN   &$
      				next = (next)+(assoc)&$
				sass = SIZE(assoc, /Dimensions) &$
    				mass = TOTAL(assoc) &$
    				xass = ROUND(TOTAL(TOTAL(assoc,2)*Indgen(sass[0]))/mass) &$
    				yass = ROUND(TOTAL(TOTAL(assoc,1)*Indgen(sass[1]))/mass) &$
				info[0,1] = label &$
				info[1,1] = (l+1) &$
				info[2,1] = xass &$
				info[3,1] = yass &$
				info[4,1] = N_ELEMENTS(WHERE(assoc gt (0.))) &$
				info[5,1] = (l+1) &$
    			ENDIF
;----------------------------------------------------------------------------------------------------------------------------
;SECTION FIVE: 	LINES 204 TO 222
;PURPOSE:   	If an associated MBP is not found, this section searchs 30seconds of data for
;		an associated MBP.
;PARAMETERS: 	NONE
;
;DESCRIPTION:
;This section sets up a search area across which the code shall attempt to find an associated
;MBP in the following 30 seconds of data.  The funciton MBP_search.pro performs this search and
;shall return a positive binary image if an MBP is found, whilst a null array indicates that the
;obj no longer exists.
;----------------------------------------------------------------------------------------------------------------------------
			IF (MAX(assoc)) eq (0.) THEN BEGIN &$
				growing = 1. 
				source = INTARR(im_size[1],im_size[2])
				xy = MBP_xyposition(obj, WHERE(obj gt (0.)))
   				IF ((MIN(xy[0,*]))-5) lt (0.) THEN (xlo) = 0. ELSE (xlo) = ((MIN(xy[0,*]))-5)
    				IF ((MIN(xy[1,*]))-5) lt (0.) THEN (ylo) = 0. ELSE (ylo) = ((MIN(xy[1,*]))-5)
				IF ((MAX(xy[0,*]))+5) gt (im_size[1]-1.) THEN (xhi) = (im_size[1])-(1.) ELSE (xhi) = ((MAX(xy[0,*]))+5)
    				IF ((MAX(xy[1,*]))+5) gt (im_size[2]-1.) THEN (yhi) = (im_size[2])-(1.) ELSE (yhi) = ((MAX(xy[1,*]))+5)
			
				IF (l+(imgs)) ge (imend) THEN (up) = ((imend)-2.) ELSE (up) = (l+(imgs))  &$
				
				for i = (l+2),(up) DO BEGIN &$
					huntim = array[*,*,(i)] &$
					IF (N_ELEMENTS(WHERE((obj+huntim) gt (1.)))) gt (2.) THEN BEGIN &$
						RofI  = REGION_GROW(huntim, WHERE((obj+huntim) gt (1.))) &$
						assoc[RofI] = huntim[RofI] &$
					ENDIF ELSE BEGIN &$
						assoc = MBP_search(huntim,xlo,xhi,ylo,yhi,xobj,yobj) &$
					ENDELSE &$
    					IF (MAX(assoc)) eq (1.) AND (N_ELEMENTS(WHERE(assoc gt 0))) ge (4) THEN BREAK &$
				ENDFOR
;----------------------------------------------------------------------------------------------------------------------------
;SECTION SIX: 	LINES 247 TO 362
;PURPOSE:   	If an associated MBP is found in the next 30 seconds then the code shall grow
;		the MBP in the frame it was missed
;PARAMETERS: 	BIN - binary image to store objects from the frame where the MBP was missed.
;		SEED - provides a seed region from which to grow the MBP
;		CACHE - stores the intensity of the turning points.
;		CONST - a 3 element array. 1st ELEMENT - Extension of the slice
;					   2nd ELEMENT - Rotation Angle Control
;					   3rd ELEMENT - Rotated image control.
;		OBJST & IMGST - Stores the objects and image in correct orienation to allow for
;				quick and realistic growing later
;		GROWN - stores the grown MBP.
;
;DESCRIPTION:
;If an associated MBP is found in the 30 seconds of data then an MBP is grown in the frame it 
;was missed. A seed region is set up from the image where the MBP was missed (lines[249:263]). 
;This region is then employed find turning points in intensity profiles from the image. Profiles 
;in 5 degree steps around the seed region are taken, providing 36 intensities.  The highest 
;turning point intensity is used as a threshold above which any conjoining pixels are consider 
;to be part of the MBP.  To restrict the growth of MBPs, so as not to include separate objects,
;the code shall increase the threshold to a level where the grown MBP does not overlap other 
;objects (lines[337:344]).
;----------------------------------------------------------------------------------------------------------------------------			
				IF (MAX(assoc)) eq (1.) AND (N_ELEMENTS(WHERE(assoc gt 0))) ge (4)  THEN BEGIN
					for q = 0., 20. DO BEGIN &$
						bin  = INTARR(im_size[1], im_size[2]) &$
						bin[WHERE(norm gt ((MEAN(norm))+((2.-(q/10.))*sigi)))] = 1. &$
						IF (N_ELEMENTS(WHERE((bin+obj) gt (1)))) ge ((N_ELEMENTS(WHERE(obj gt 0)))/2) AND (N_ELEMENTS(WHERE((bin+obj) gt (1)))) ge (4) THEN BEGIN &$
							reg  = REGION_GROW(bin,WHERE(bin+obj gt (1.))) &$
							source[reg]=bin[reg] &$
							intim = (source)*(norm) &$
							sig = SQRT((MOMENT(intim[WHERE((intim) gt (0.) AND (intim) lt (10.))]))(1)) &$
							IF (((MAX(intim))-(MIN(intim[WHERE(intim gt (0))])))/sig) gt 3 THEN BEGIN &$
								splgrp = MBP_split(intim,im_size,sig) &$
								bin = ((bin-source)+splgrp) &$
								source = MBP_seed(obj,splgrp,xobj,yobj,im_size) &$
							ENDIF &$
						ENDIF &$
						IF (MAX(source)) gt (0) THEN BREAK &$
					ENDFOR
										
					IF (N_ELEMENTS(WHERE((source[*,*,0]+follow) gt (1.)))) gt (2.) THEN BEGIN &$
						ovlp = INTARR(im_size[1], im_size[2]) &$
						RofI  = REGION_GROW(follow, WHERE((source[*,*,0]+follow) gt (1.))) &$
						ovlp[RofI] = follow[RofI] &$
						ovlp = MBP_close(ovlp,xobj,yobj)&$
						follow = (follow)-(ovlp) &$
 					ENDIF

					IF (MAX(source)) gt (0) THEN BEGIN
						IF (N_ELEMENTS(SIZE(source))) eq (6) THEN (sources)=((SIZE(source))(3)) ELSE (sources) = 1.

						for g = 0, (sources-1) DO BEGIN
							seed = source[*,*,(g)]
							ss = SIZE(seed, /Dimensions) &$	
    							ms = TOTAL(seed) &$	
    							xs = ROUND(TOTAL( TOTAL (seed,2) * Indgen(ss[0]) ) / ms) &$	
    							ys = ROUND(TOTAL( TOTAL (seed,1) * Indgen(ss[1]) ) / ms) &$	

							cache    = FLTARR(36)
							const    = INTARR(3)
							const[0] = 20
							objst    = seed
							imgst	 = norm
							bin	 = (bin)-(seed)

							WHILE (const[1]) lt (36) DO BEGIN &$
								tpt  = MAKE_ARRAY(2,VALUE=-1) &$
								inter = INTARR(im_size[1], im_size[2]) &$
								line = INTARR(im_size[1],im_size[2]) &$ 
								IF (ys-const[0]) lt (0.) THEN (ylow) = 0. ELSE (ylow) = (ys-const[0]) &$ 
								IF (ys+const[0]) gt (im_size[2]-1) THEN (yhig) = (im_size[2]-1) ELSE (yhig) = (ys+const[0]) &$
								line[xs,(ylow):(yhig)] = 1. &$

								IF (const[1]) lt (19) THEN (rot) = ROT(line,(const[1]*5),1,xs,ys,MISSING = 0.,/PIVOT) ELSE (rot) = ROT(line,((const[1]*5)-90),1.0,xs,ys,MISSING=0.,/PIVOT) &$
								slice = (norm)*(rot) &$
								arr = slice[WHERE(slice gt 0)] &$
								dx = MBP_GRADIENT(SMOOTH(arr,2)) &$

								IF (MAX(seed+rot)) gt (1.) THEN (inter[WHERE((seed+rot) eq (2))]) = (1.) &$
								inter = (inter)+(rot) &$
								ins = inter[WHERE(inter gt (0.))] &$
								IF (MAX(seed+rot)) gt (1.) THEN (maxpt) = WHERE(arr eq (MAX(arr[WHERE(ins eq 2)]))) ELSE (maxpt) = WHERE(arr eq (MAX(arr[((N_ELEMENTS(arr)/2)-5):((N_ELEMENTS(arr)/2)+5)]))) &$
								IF ((WHERE(ins gt (1.)))(0)) gt (-1) THEN (ins[(MIN(WHERE(ins eq (2)))):(MAX(WHERE(ins eq (2))))]) = 2. &$

								for j = (3.), (N_ELEMENTS(dx)-4) DO BEGIN &$
									IF (ins[j]) lt (2.) AND (dx[j-1]) lt (0.) AND (dx[j]) gt (0.) AND  (TOTAL(dx[(j-3):(j-1)])) le (0.) AND (TOTAL(dx[(j+1):(j+3)])) ge (0.) THEN BEGIN &$
										IF (j) lt (maxpt[0]) AND (j) gt (tpt[0]) THEN (tpt[0]) = (j) &$
										IF (j) gt (maxpt[0]) AND (tpt[1]) eq (-1) THEN (tpt[1]) = (j) &$
										IF (tpt[0]) gt (-1) AND (tpt[1]) gt (-1) THEN BREAK &$
									ENDIF &$
								ENDFOR &$

								IF ((WHERE((rot+edg) eq (2)))(0)) eq ((WHERE(rot eq 1))(0)) AND (tpt[0]) eq (-1.) THEN (tpt[0]) = 0.  &$
								IF ((WHERE((rot+edg) eq (2)))((N_ELEMENTS(WHERE((rot+edg) eq (2))))-1)) eq ((WHERE(rot eq 1))((N_ELEMENTS(WHERE(rot eq 1)))-1)) AND (tpt[1]) eq (-1.) THEN (tpt[1]) = 0.  &$
								IF (WHERE((line+edg) eq (2))) eq ((WHERE(line eq 1))(0)) AND (tpt[0]) eq (-1.) THEN (tpt[0]) = 0.  &$
								IF (WHERE((line+edg) eq (2))) eq ((WHERE(line eq 1))((N_ELEMENTS(WHERE(line eq 1)))-1)) AND (tpt[1]) eq (-1.) THEN (tpt[1]) = 0.  &$

				 				IF (tpt[0]) eq (-1) OR (tpt[1]) eq (-1) THEN (const[0]) = (const[0]) + 20. &$

								IF (tpt[0]) gt (-1) AND (tpt[1]) gt (-1) THEN BEGIN &$
									IF (arr[tpt[0]]) gt (arr[tpt[1]]) AND (tpt[0]) gt (0.) THEN (cache[const[1]]) = (arr[tpt[0]]) ELSE (cache[const[1]])=(arr[tpt[1]]) &$
									const[0] = 20. &$
									const[1] = (const[1])+(1) &$

									IF (const[1]) eq (19.) THEN BEGIN &$
										const[2]=1. &$
										seed = ROTATE(seed,1) &$
										edg  = ROTATE(edg,1) &$
										norm = ROTATE(norm,1) &$
										im_size = SIZE(norm) &$
										ss = SIZE(seed, /Dimensions) &$       
 	 									ms = TOTAL(seed)  &$
   										xs = ROUND(TOTAL(TOTAL(seed,2)*INDGEN(ss[0]))/ms)  &$   
   										ys = ROUND(TOTAL(TOTAL(seed,1)*INDGEN(ss[1]))/ms) &$
									ENDIF &$
								ENDIF &$
							ENDWHILE &$

							RofI = REGION_GROW(imgst,(WHERE(objst eq (1))),THRESHOLD=[MAX(cache),MAX(imgst)])
							WHILE (N_ELEMENTS(RofI)) lt (4) DO BEGIN &$
								cache[WHERE(cache eq MAX(cache))] = 0. &$
								RofI = REGION_GROW(imgst,(WHERE(objst eq (1))),THRESHOLD=[MAX(cache),MAX(imgst)]) &$
							ENDWHILE 

							IF (RofI[0]) gt (-1.) THEN BEGIN
								grown[RofI]=imgst[RofI]
								grown[WHERE(grown gt (0.))]=1.
								vary=0.001
								
								WHILE (MAX(grown+bin)) gt (1.) OR (MAX(grown+next)) gt (1) OR (MAX(grown+follow)) gt (1) DO BEGIN &$
									grown = FLTARR(im_size[2],im_size[1]) &$
									RofI = REGION_GROW(imgst,(WHERE(objst eq (1))), THRESHOLD=[(MAX(cache)+vary),(MAX(imgst))]) &$
									IF (RofI[0]) eq (-1.) THEN BREAK &$
									grown[RofI]=imgst[RofI] &$
									grown[WHERE(grown gt (0.))]=1. &$
									grown = MBP_close(grown,xobj,yobj) &$
									vary = vary+0.001 &$
								ENDWHILE
								
								;NEED TO CHANGE FOR OTHER DATA SETS (200 for ROSA) (272 for SST) (400 for Sims)
								IF (N_ELEMENTS(WHERE(grown gt 0))) ge (4.) AND (N_ELEMENTS(WHERE(grown gt 0))) lt (272)  THEN BEGIN &$
									follow = (follow)+(grown) &$
									IF (g) eq (0) THEN BEGIN &$
										next = (next)+(grown) &$
										sgr = SIZE(grown,/DIMENSIONS) &$
										mgr = TOTAL(grown) &$
										xgr = ROUND(TOTAL( TOTAL(grown, 2) * Indgen(sgr[0]) )/mgr) &$	
										ygr = ROUND(TOTAL( TOTAL(grown, 1) * Indgen(sgr[1]) )/mgr) &$
										info[0,1] = (label) &$
										info[1,1] = (l+1) &$
										info[2,1] = xgr &$
										info[3,1] = ygr &$
										info[4,1] = N_ELEMENTS(WHERE(grown gt (0.))) &$
										info[5,1] = (i) &$
									ENDIF &$
								ENDIF
							ENDIF
							
							IF (const[2]) eq (1.) THEN BEGIN &$
								seed = ROTATE(seed,3) &$
								edg = ROTATE(edg,3) &$
								norm = ROTATE(norm,3) &$
								im_size = SIZE(norm) &$
								const[2] = 0. &$
							ENDIF
							
							bin = (bin)+(seed)
						ENDFOR
					ENDIF
				ENDIF
;----------------------------------------------------------------------------------------------------------------------------
;SECTION FIVE: 	LINES 374 TO 393
;PURPOSE:   	If an associated MBP is not found in the next 30 seconds then the code shall
;		consider the MBP to be dead
;PARAMETERS: 	NONE
;
;DESCRIPTION:
;If an associated MBP is not found in the follwoing 30secs of data then the MBP is considered to
;be dead and is removed from the information store.  This section also rotates the images back
;to the correct orientaiton if they had been rotated (line[380:385]).
;----------------------------------------------------------------------------------------------------------------------------
				;NEED TO CHANGE VALUE BASED ON RES (e.g. 200 for ROSA 272 for SST 400 for sims)
				IF (MAX(assoc)) eq (0.) OR (MAX(source)) eq (0.) OR (N_ELEMENTS(WHERE(assoc gt (0.)))) lt (4.) OR (N_ELEMENTS(WHERE(grown gt (0.)))) lt (4.) OR (N_ELEMENTS(WHERE(grown gt (0.)))) gt (272.)  THEN BEGIN &$
				;IF (MAX(assoc)) eq (0.) OR (MAX(source)) eq (0.) OR (N_ELEMENTS(WHERE(assoc gt (0.)))) lt (4.) OR (N_ELEMENTS(WHERE(grown gt (0.)))) lt (4.) OR (N_ELEMENTS(WHERE(grown gt (0.)))) gt (400.)  THEN BEGIN &$
					info[*,*] = -1. &$
					label     = (label)-(1) &$
					frame 	  = (frame)-(obj) &$
					IF (growing) eq (0) AND (MAX(next+assoc)) gt (1) THEN (next) = (next)-(assoc)
					IF (growing) eq (1) AND (MAX(next+grown)) gt (1) THEN (next) = (next)-(grown)
				ENDIF
			ENDIF
		ENDIF
		
		IF (MIN(info[0:5,*])) gt (-1.) THEN BEGIN &$
			OPENU, 1, 'table.dat',/APPEND &$
			printf,1,info &$
			CLOSE,1 &$
		ENDIF
;----------------------------------------------------------------------------------------------------------------------------
;SECTION SIX: 	LINES 405 TO 418
;PURPOSE:   	To identify any merging objects and remove track those that are closest to the
;		final outcome
;PARAMETERS: 	TEMP - stores the NEXT parameter, which stores the associated MBPs.
;
;DESCRIPTION:
;The procedure merger_check, shall  use the information from the ''+FILE+'' to search for objects
;that are tracked to the same location and determine which of these are closest to the final
;mergered object and consider it as the continuation of the obect 
;----------------------------------------------------------------------------------------------------------------------------
		;IF (MAX(next)) gt (1.) THEN BEGIN
		;	temp = next
		;	IF (MAX(grown)) eq (0.) THEN MBP_merger,''+file+'',next,assoc,prev,xprev,yprev,(xobj),(yobj),(l)
		;	IF (MAX(grown)) eq (1.) THEN MBP_merger,''+file+'',next,grown,prev,xprev,yprev,(xobj),(yobj),(l)
		;	next = temp
		;	next[WHERE(next gt (1.))] = 1.
		;ENDIF
	ENDWHILE
	stable[*,*,(l-imst)] = frame
	stable[*,*,((l-imst)+1)] = next
	
	cls
	writeu, -1, 'Processed frame    ', sl+string(pixel, format="(i4)"), ' of    ', sl+string(allofthem, format="(i4)"),' '
ENDFOR

IF (FILE_TEST('table.dat')) THEN FILE_DELETE,'table.dat'

RETURN,stable

end
