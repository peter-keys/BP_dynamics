;+
;
; ROUTINE:	MBP_tabulating
;
; PURPOSE:    	TO TRACK MBPS FROM AN IMAGE TO THE SUCCEEDING IMAGE AND STORE THE ASSOCIATED
;		MBPS PROPERTIES 
;		
; USEAGE:   	FILE, NORM, OBJS, ARR, IMGS, IMNUM, IMEND, SIGI	  	
; INPUT:	FILE - FILENAME OF DATA
;		NORM - NORMALISED FRAME OF THE FOLLOWING IMAGE
;		OBJS - BINARY ARRAY OF THE OBJS TO BE INVESTIGATED
;		ARRAY - DETECTION ARRAY TO SEARCH ACROSS FOR MISSING OBJECTS
;		IMGS - THE NUMBER OF IMGS TO MAKE UP 30 SECONDS OF DATA
;		IMNUM - THE NUMBER OF THE IMAGE CURRENTLY UNDER INVESTIGATION
;		IMEND - THE FINAL IMAGE TO BE CONSIDERED
;		SIGI - THE SIGMA VALUE OF IMAGE NORM
; OUTPUT:	OBJS- RETURNS THE NEW OBJECTS THAT APPEAR IN THE IMAGE TO BE TRACKED
; CALLING	CALLED FROM MBP_tracking.pro [line 132]
; AUTHOR:     	Version 2.3 
;
;AMENDMENT HISTORY
; (v2.1)-	Higher initial intensity threshold employed to set a seed region.  Now potential seed regions must
;		     overlap original object by 1/2 the number of pixels prior to being separated.  This stops the inclusion
;		     of single pixel seeds.
; (v2.1)-	Adjustable growing threshold to ensure every chance is given to grow a missed MBP.
; (v2.1)-	Grown array is now set to zero during varying stage to prevent erronous inclusion of small MBPs
; (v2.2)-	Addition of MBP_seed program which handles multiple seed regions if object is observed to separate.
; (v2.3)-	The mergin of two objects found with an associated object is now recognised
; (v2.3)-	THe merging of objects is now noted with the addition of a 7th column to table.


FUNCTION MBP_tabulating, file, norm, objs, array, next,edg, imgs, imnum, imend,sigi,im_size
;----------------------------------------------------------------------------------------------------------------------------
;SECTION ONE: 	LINES 30 TO 42
;PURPOSE: 	To read in the previously tabulated data and track the detected MBPs to
;		the next image.
;PARAMETERS: 	OLD_INFO - array storing the tabulated data from the previous frame
;		NEXT - binary array to store the objects tracked to the next frame.
;		FOLLOW - the following detection frame
;DESCRIPTION:	
;Section one reads in the information about the previously detected objects from file
;[line 39] and records the data from the last frame in parameter OLD_INFO[lines 40:45]
;The [x,y] size of the image is calculated and employed to set up array NEXT [line 49],
;which shall store all associated objects from the following frame.  FOLLOW stores the
;next binary detection image and is used to search for associated MBPs.
;----------------------------------------------------------------------------------------------------------------------------
readcol,''+file+'',f='I,I,I,I,x,I',num,frmnum,xco,yco,fondin,/sil
old_info = INTARR(5,N_ELEMENTS(WHERE(frmnum eq (imnum))))
old_info[0,*] = num[WHERE(frmnum eq (imnum))]
old_info[1,*] = frmnum[WHERE(frmnum eq (imnum))]
old_info[2,*] = xco[WHERE(frmnum eq (imnum))]
old_info[3,*] = yco[WHERE(frmnum eq (imnum))]
old_info[4,*] = fondin[WHERE(frmnum eq (imnum))]

next = INTARR(im_size[1],im_size[2])
objsto = objs
follow = array[*,*,(imnum+1)]
;----------------------------------------------------------------------------------------------------------------------------
;SECTION TWO: 	LINES 56 TO 69
;PURPOSE: 	To identify and grow the MBP that have already been detected and label
;PARAMETERS: 	OBJ - array storing the individual object that is being considered
;		ASSOC - array to store the associated object if found
;		INFO - array to hold the information about the associated
;		       object.[ID_num,frame_num, xcoord,ycoord,area,frame found in]
;DESCRIPTION:	
;Section two goes through each object from the previous image and reads in the old [x,y] 
;co-ordinate from the previous frame and grows the MBP that was located there.  If the 
;[x,y] position is outside the object it shall grow the closest object.  Function 
;MBP_close [line 75] facilitates this.
;----------------------------------------------------------------------------------------------------------------------------
for z = 0, ((N_ELEMENTS(WHERE(frmnum eq (imnum))))-1) DO BEGIN
	obj=INTARR(im_size[1],im_size[2]) 
	assoc = INTARR(im_size[1],im_size[2])
	grown = FLTARR(im_size[1],im_size[2])
 	growing = 0
	info = MAKE_ARRAY(7,1,/INT,VALUE=-1)
	pix = MBP_xyconverter(objs, old_info[2:3,(z)])
	
	IF (objs[pix]) eq (1.) THEN BEGIN &$
		ROFI = REGION_GROW(objs,pix)  &$
		obj[ROFI] = objs[ROFI] &$
	ENDIF ELSE BEGIN   &$
		obj = MBP_close(objs, old_info[2,(z)], old_info[3,(z)]) &$
	ENDELSE
		
	objs = (objs) - (obj)	
	
	IF (N_ELEMENTS(WHERE(obj gt (0)))) ge (4) THEN BEGIN
;----------------------------------------------------------------------------------------------------------------------------
;SECTION THREE:	LINES 38 TO 57
;PURPOSE: 	To search for associated (overlapping) objects in the following image.
;PARAMETERS: 	NONE
;DESCRIPTION:	
;Section three searches the following detection image, stored in array FOLLOW, for any
;overlapping objects.  These are then considered the associated object and grown and
;stored in array ASSOC
;----------------------------------------------------------------------------------------------------------------------------
		IF (N_ELEMENTS(WHERE((obj+follow) gt (1.)))) ge (2.) THEN BEGIN &$
			RofI  = REGION_GROW(follow, WHERE((obj+follow) gt (1.))) &$
			assoc[RofI] = follow[RofI] &$
			assoc = MBP_close(assoc,old_info[2,(z)],old_info[3,(z)])&$
			IF (N_ELEMENTS(WHERE(assoc gt 0))) lt (4) THEN (assoc)=INTARR(im_size[1],im_size[2])
		ENDIF
;----------------------------------------------------------------------------------------------------------------------------
;SECTION FOUR:	LINES 97 TO 109
;PURPOSE: 	If an associated object is found in the succeeding image, its information
;		is stored.
;PARAMETERS: 	NONE
;DESCRIPTION:	
;If an assocaited MBP is found in the immediately succedding detect image, section four
;shall calculated the associated [x,y] [i.e. xass, yass paratmeters[lines 108:111]], and
;store them in array INFO.  The elements of INFO are [BTPT ID_NUM,
;FRAME_NUM,XCORD,YCORD,AREA,FRAME FOUND IN].  Thus, the BTPT ID_NUM shall be the same as
;the previous ID_NUM [line 112] the FRAME_NUM shall increment by one [line 112] AREA is
;the number of elements/pixels that make up the object [line 116] and the frame it was
;found in is the same as FRAME_NUM [line 117]
;----------------------------------------------------------------------------------------------------------------------------
		IF (MAX(assoc)) gt (0.) THEN BEGIN   &$
      			next = (next)+(assoc)&$
			sass = SIZE(assoc, /Dimensions) &$
    			mass = TOTAL(assoc) &$
    			xass = ROUND(TOTAL(TOTAL(assoc,2)*Indgen(sass[0]))/mass) &$
    			yass = ROUND(TOTAL(TOTAL(assoc,1)*Indgen(sass[1]))/mass) &$
			info[0,0] = (old_info[0,(z)]) &$
			info[1,0] = ((old_info[1,(z)])+1) &$
			info[2,0] = xass &$
			info[3,0] = yass &$
			info[4,0] = N_ELEMENTS(WHERE(assoc gt (0.))) &$
			info[5,0] = (imnum+1)  &$		
		ENDIF
;----------------------------------------------------------------------------------------------------------------------------
;SECTION FIVE:	LINES 124 TO 147
;PURPOSE: 	If an associated object is not found in the succeeding image, the
;		algorithm shall search for it across the following 30 seconds of data
;PARAMETERS: 	NONE
;DESCRIPTION:	
;If an assocaited MBP is not found in the immediately succedding detect image, section
;five shall hun across the following 30seconds of data for an associated object.  It
;performs this task in two ways.  Firstly it looks for any overlapping objects in the
;images [lines 149:152].  Secondly, if none are found, it scans an 10*10 pixel area
;located around the MBP.  This area is setup in [lines 134:138] and the function
;MBP_search, performs the scan.  Note the parameter IMGS [line 145] is employed here to 
;stiplulate the number of images across which to scan.
;----------------------------------------------------------------------------------------------------------------------------
		IF (MAX(assoc)) eq (0.) THEN BEGIN &$ 
			source  = INTARR(im_size[1], im_size[2])
			growing = 1.
			IF (old_info[4,(z)]) gt (imnum) THEN (assoc[*,*]) = 1. ELSE BEGIN
				xy = MBP_xyposition(obj, WHERE(obj gt (0.)))
   				IF ((MIN(xy[0,*]))-5) lt (0.) THEN (xlo) = 0. ELSE (xlo) = ((MIN(xy[0,*]))-5)
    				IF ((MIN(xy[1,*]))-5) lt (0.) THEN (ylo) = 0. ELSE (ylo) = ((MIN(xy[1,*]))-5)
				IF ((MAX(xy[0,*]))+5) gt (im_size[1]-1.) THEN (xhi) = (im_size[1])-(1.) ELSE (xhi) = ((MAX(xy[0,*]))+5)
    				IF ((MAX(xy[1,*]))+5) gt (im_size[2]-1.) THEN (yhi) = (im_size[2])-(1.) ELSE (yhi) = ((MAX(xy[1,*]))+5)

				IF (imnum+(imgs)) ge (imend) THEN (up) = ((imend)-2.) ELSE (up) = (imnum+(imgs))  &$

				for i = (imnum+2),(up) DO BEGIN &$
					huntim = array[*,*,(i)] &$
					IF (N_ELEMENTS(WHERE((obj+huntim) gt (1.)))) gt (2.) THEN BEGIN &$
						RofI  = REGION_GROW(huntim, WHERE((obj+huntim) gt (1.))) &$
						assoc[RofI] = huntim[RofI] &$
					ENDIF ELSE BEGIN &$
						assoc = MBP_search(huntim,xlo,xhi,ylo,yhi,old_info[2,(z)],old_info[3,(z)]) &$
					ENDELSE &$
    					IF (MAX(assoc)) eq (1.) AND (N_ELEMENTS(WHERE(assoc gt 0))) ge (4) THEN BREAK &$
				ENDFOR
			ENDELSE
;----------------------------------------------------------------------------------------------------------------------------
;SECTION SIX:	LINES 176 TO 291
;PURPOSE: 	If an associated object is found in the following 30 seconds of data then
;		the algorithm locates a seed region in the image it was missed and grows
;		the missing MBP
;PARAMETERS: 	BIN - 	stores the objects in the frame where the MBP was missed. This
;		      	array holds the candidates to be seeds
;		SEED -	a seed object from which the missng MBP is grown
;		CACHE - stores the 36 individual thresholds to be used as a limit in
;			growing the object
;		CONST - a three element array storing the slice extenstion [Element 1],
;			the rotation angle [Element 2] and the image rotation alert
;			[Element 3]
;		OBJST & IMGST - store the object and the image in the original
;				orientation
;		TPT -	two element array to store the turning points location
;		INTER - array used to store the internal intensity structure across the
;			slice
;		GROWN - array to store the grown object
;DESCRIPTION:
;If an associated MBP is found in the 30 seconds of data then an MBP is grown in the frame it 
;was missed. A seed region is set up from the image where the MBP was missed (lines[188:202]). 
;This region is then employed find turning points in intensity profiles from the image. Profiles 
;in 5 degree steps around the seed region are taken, providing 36 intensities.  The highest 
;turning point intensity is used as a threshold above which any conjoining pixels are consider 
;to be part of the MBP.  To restrict the growth of MBPs, so as not to include separate objects,
;the code shall increase the threshold to a level where the grown MBP does not overlap other 
;objects (lines[276:283]).  The parameters of the grown MBP are stored in INFO [lines 286:299]
;----------------------------------------------------------------------------------------------------------------------------
			IF (MAX(assoc)) eq (1.) AND (N_ELEMENTS(WHERE(assoc gt 0))) ge (4) THEN BEGIN
				for q = 0., 20. DO BEGIN &$
					bin  = INTARR(im_size[1], im_size[2]) &$
					bin[WHERE(norm gt ((MEAN(norm))+((2.-(q/10.))*sigi)))] = 1. &$
					IF (N_ELEMENTS(WHERE((bin+obj) gt (1)))) ge ((N_ELEMENTS(WHERE(obj gt 0)))/2) THEN BEGIN &$
						reg  = REGION_GROW(bin,WHERE(bin+obj gt (1.))) &$
						source[reg]=bin[reg] &$
						intim = (source)*(norm) &$
						sig = SQRT((MOMENT(intim[WHERE((intim) gt (0.) AND (intim) lt (10.))]))(1)) &$
						IF (((MAX(intim))-(MIN(intim[WHERE(intim gt (0))])))/sig) gt 3 THEN BEGIN &$
							splgrp = MBP_split(intim,im_size,sig) &$
							bin = ((bin-source)+splgrp) &$
							source = MBP_seed(obj,splgrp,old_info[2,(z)],old_info[3,(z)],im_size) &$
						ENDIF &$
					ENDIF &$
					IF (MAX(source)) gt (0) THEN BREAK &$
				ENDFOR

				IF (N_ELEMENTS(WHERE((source[*,*,0]+follow) gt (1.)))) gt (2.) THEN BEGIN &$
					ovlp  = INTARR(im_size[1], im_size[2])
					RofI  = REGION_GROW(follow, WHERE((source[*,*,0]+follow) gt (1.))) &$
					ovlp[RofI] = follow[RofI] &$
					ovlp = MBP_close(ovlp,old_info[2,(z)],old_info[3,(z)])&$
					follow = (follow)-(ovlp)
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

							IF (const[1]) lt (19) THEN (rot) = ROT(line,(const[1]*5),1,xs,ys,missing = 0.,/PIVOT) ELSE (rot) = ROT(line,((const[1]*5)-90),1.0,xs,ys,MISSING=0.,/PIVOT) &$
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
								RofI = REGION_GROW(imgst,(WHERE(objst eq (1))),THRESHOLD=[(MAX(cache)+vary),(MAX(imgst))]) &$
								IF (RofI[0]) eq (-1.) THEN BREAK &$
								grown[RofI]=imgst[RofI] &$
								grown[WHERE(grown gt (0.))]=1. &$
								grown = MBP_close(grown,old_info[2,(z)],old_info[3,(z)]) &$ 
								vary = vary+0.001 &$
							ENDWHILE

							IF (N_ELEMENTS(WHERE(grown gt 0))) ge (4.) AND (N_ELEMENTS(WHERE(grown gt 0))) lt (200) THEN BEGIN &$
								follow = (follow)+(grown) &$
								IF (g) eq (0) THEN BEGIN &$
									next = (next)+(grown) &$
									sgr = SIZE(grown,/DIMENSIONS) &$
									mgr = TOTAL(grown) &$
									xgr = ROUND(TOTAL( TOTAL(grown, 2) * Indgen(sgr[0]) )/mgr) &$	
									ygr = ROUND(TOTAL( TOTAL(grown, 1) * Indgen(sgr[1]) )/mgr) &$
									info[0,0] = (old_info[0,(z)]) &$
									info[1,0] = (imnum+1) &$
									info[2,0] = xgr &$
									info[3,0] = ygr &$
									info[4,0] = N_ELEMENTS(WHERE(grown gt (0.))) &$
									IF (old_info[4,(z)]) gt (imnum) THEN (info[5,0]) = (old_info[4,(z)]) ELSE (info[5,0]) = (i) &$
								ENDIF
							ENDIF &$

							IF (const[2]) eq (1.) THEN BEGIN &$
								seed = ROTATE(seed,3) &$
								edg = ROTATE(edg,3) &$
								norm = ROTATE(norm,3) &$
								im_size = SIZE(norm) &$
								const[2]=0.
							ENDIF
						
							bin = (bin)+(seed)
						ENDIF
					ENDFOR
				ENDIF
			ENDIF
;----------------------------------------------------------------------------------------------------------------------------
;SECTION SEVEN:	LINES 305 TO 306
;PURPOSE: 	If an associated object is not found in the following 30 seconds of data then
;		the algorithm assumes that the object no longer exists.
;PARAMETERS: 	NONE
;DESCRIPTION:
;If an associated MBP is not found in the 30 seconds of data, or the grown MBP is too small, 
;then the object is considere to no longer exist and the information is not recorded [lines
;301:304] and [lines 316:320].  If the image was rotated, this section also realigns the images
;to the original orientation.  Finally, any objects that were not tracked from the previous
;frame, (i.e. any objects that are new objects), are returned to the tracking process [line
;328].
;----------------------------------------------------------------------------------------------------------------------------	
			IF (MAX(assoc)) eq (0.) OR (MAX(source)) eq (0.) OR (N_ELEMENTS(WHERE(assoc gt (0.)))) lt (4.) OR (N_ELEMENTS(WHERE(grown gt (0.)))) lt (4.) OR (N_ELEMENTS(WHERE(grown gt 0))) gt (200) THEN BEGIN
				info[*,*] = -1. &$
				IF (growing) eq (0) AND (MAX(next+assoc)) gt (1) THEN (next) = (next)-(assoc)
				IF (growing) eq (1) AND (MAX(next+grown)) gt (1) THEN (next) = (next)-(grown)
			ENDIF
		ENDIF
	ENDIF
	
	IF (MIN(info[0:5])) gt (-1.) THEN BEGIN &$
		OPENU, 1, 'table.dat',/APPEND &$
		printf,1,info &$
		CLOSE,1 &$
	ENDIF
	
	IF (MAX(next)) gt (1.) THEN BEGIN
		temp = next
		IF (MAX(grown)) eq (0.) THEN MBP_merger,''+file+'',next,assoc,(old_info[0,*]),(old_info[2,*]),(old_info[3,*]),(old_info[2,(z)]), (old_info[3,(z)]), (imnum)
		IF (MAX(grown)) eq (1.) THEN MBP_merger,''+file+'',next,grown,(old_info[0,*]),(old_info[2,*]),(old_info[3,*]),(old_info[2,(z)]), (old_info[3,(z)]), (imnum)
		next = temp
		next[WHERE(next gt (1.))] = 1.
	ENDIF
ENDFOR
str = {objs:objs,follow:follow,next:next}
RETURN, str

end
;----------------------------------------------------------------------------------------------------------------------------
;----------------------------------------------------------------------------------------------------------------------------
