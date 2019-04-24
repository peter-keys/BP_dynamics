;+
;
; ROUTINE: 	MBP_threshold
;
; PURPOSE: 	FINDS A LIMIT ON THE INTENSITY GRADIENT REQUIRED FOR MBP
;		
; USEAGE:  	NORM_IM,OBJS, EDG, SPLT, IM_SIZE
; INPUT:   	NORM_IM - IMAGE BEING INVESTIGATED NORMLISED TO THE MEAN 
;		OBJS - THE OBJECTS ACROSS WHICH TO TAKE INTENSITY PROFILES
;		EDG - TO REMOVE THOSE OBJECTS CLOSE TO THE EDGE WHICH FAIL TO 
;		      PROVIDE TWO TURING POINTS IN THEIR PROFILES
;		SPLT - ARRAY TO STORE THE SEPARATED OBJECTS SO THEY SHALL NOT BE
;		       DIVIDED AGAIN	  	
;		IM_SIZE - PROVIDES THE REQUIRED [X,Y] DIMENSIONS TO SET UP 
;			  VARIOUS PARAMETERS
; OUTPUT:  	THRES - A THRESHOLD INTENSITY GRADIENT EMPLOYED TO IDENTIFY
;			MBPS IN MBP_detect
; CALLING:	CALLED IN MBP_detect.pro
; AUTHOR:  	Version 2.0 - Philip. J. Crockett, QUB, 25 Nov 2010
;   	    	         (Email: pcrockett02@qub.ac.uk)

FUNCTION MBP_threshold,norm_im,objs,edg,splt,im_size

first = -1
count = -1
splst = INTARR(im_size[1], im_size[2])
;====================================================================================================================================
;LOOK FOR GRADIENTS
;====================================================================================================================================			
WHILE (count) lt (300) AND (MAX(objs)) gt (0.) DO BEGIN
	grad = FLTARR(4)
	obj = INTARR(im_size[1],im_size[2])
	ROI = REGION_GROW(objs,(WHERE(objs*norm_im eq (MAX(objs*norm_im)))))
	obj[ROI]=objs[ROI]
	splst = (splst)+(obj)
	objs = (objs)-(obj)

	IF (N_ELEMENTS(WHERE(objs gt (0.)))) ge (4.) AND (MAX(splt+obj)) eq (1.) THEN BEGIN
		splt=(splt)+(obj)
		intim = (obj)*(norm_im)
		siga = SQRT((MOMENT(intim[WHERE((intim) gt (0.) AND (intim) lt (10.))]))(1))
		IF (((MAX(intim))-(MIN(intim[WHERE(intim gt 0.)])))/siga) gt (3.) THEN BEGIN
			splgrp = MBP_SPLIT(intim,im_size,siga)
			splst  = ((splst-obj)+splgrp)
			objs   = (objs)+(splgrp)
			obj    = (obj)-(obj)
			reg    = REGION_GROW(splgrp,((WHERE((splgrp*norm_im) eq (MAX(splgrp*norm_im))))(0)))	
			obj[reg]=splgrp[reg]
			objs = (objs) - (obj)
		ENDIF
	ENDIF

	IF (N_ELEMENTS(WHERE(obj eq (1.)))) lt (300) AND (N_ELEMENTS(WHERE(obj eq (1.)))) ge (4.)THEN BEGIN &$
		dist = INTARR(4)
		const = INTARR(3)
		const[0] = 20.

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

			IF (const[1]) lt (3) THEN (rot) = ROT(line,(const[1]*45),1,x,y,missing = 0,/PIVOT) ELSE (rot) = ROT(line,(45),1.0,x,y,MISSING=0.,/PIVOT) &$
			IF (MAX(rot+edg)) eq (2) THEN BREAK &$
			slice = (norm_im)*(rot) &$
			arr = slice[WHERE(slice gt 0)] &$
			maxpt = WHERE(arr[(const[0]-5):(const[0]+5)] eq (MAX(arr[(const[0]-5):(const[0]+5)])))+(const[0]-5) &$
			dx = MBP_GRADIENT(SMOOTH(arr,2)) &$
			
			IF (MAX(obj+rot)) gt (1.) THEN inter[WHERE((obj+rot) eq (2))] = (1.) &$
			inter = (inter)+(rot) &$
			ins = inter[WHERE(inter gt (0.))] &$
					
			IF ((WHERE(ins gt (1.)))(0)) gt (-1) THEN (ins[(MIN(WHERE(ins eq (2)))):(MAX(WHERE(ins eq (2))))]) = 2. &$

			for j = (3.), (N_ELEMENTS(dx)-4) DO BEGIN &$
				IF (ins[j]) lt (2.) AND (dx[j-1]) lt (0.) AND (dx[j]) gt (0.) AND (TOTAL(dx[(j-3):(j-1)])) le (0.) AND (TOTAL(dx[(j+1):(j+3)])) ge (0.) THEN BEGIN &$
					IF (j) lt (maxpt[0]) AND (j) gt (tpt[0]) THEN (tpt[0]) = (j) &$
					IF (j) gt (maxpt[0]) AND (tpt[1]) eq (-1) THEN (tpt[1]) = (j) &$
					IF (tpt[0]) gt (-1) AND (tpt[1]) gt (-1) THEN BREAK &$
				ENDIF &$
			ENDFOR &$

			IF (tpt[0]) eq (-1) OR (tpt[1]) eq (-1) THEN (const[0]) = (const[0]) + 20. &$
		
			IF (tpt[0]) gt (-1) AND (tpt[1]) gt (-1) THEN BEGIN &$
				maxdx = MAX(MBP_GRADIENT(arr[tpt[0]:tpt[1]]),/ABSOLUTE)
				IF (maxdx) lt (0.) THEN (maxdx)=(maxdx)*(-1.)
				grad[const[1]] = maxdx  &$
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
		IF (N_ELEMENTS(WHERE(grad gt (0.)))) eq (4) THEN (count) = (count)+(1)
		IF (N_ELEMENTS(WHERE(grad gt (0.)))) eq (4) THEN (first) = (first)+(1) 
		IF (N_ELEMENTS(WHERE(grad gt (0.)))) eq (4) AND (first) eq (0.) THEN (grads) = grad &$
		
		IF (N_ELEMENTS(WHERE(grad gt (0.)))) eq (4) AND (first) gt (0.) THEN BEGIN &$
			mem = grads &$
			grads = FLTARR(N_ELEMENTS(grad)+N_ELEMENTS(mem)) 
			grads[(0):((N_ELEMENTS(mem))-1)] = mem 
			grads[(N_ELEMENTS(mem)):*] = grad 
		ENDIF &$
	
		IF (const[2]) eq (1.) THEN BEGIN &$
			obj = ROTATE(obj,3) &$
			edg = ROTATE(edg,3) &$
			norm_im = ROTATE(norm_im,3) &$
			im_size = SIZE(norm_im) &$
			const[2] = 0.
		ENDIF
	ENDIF
;====================================================================================================================================
;IF BTPT PASSES THE GRADIENT THEN GROW THE BTPT
;-You have to test the value 'sigi' for different data sets to get best results.
;-Different values work best for different data sets. You have to visullay inspect 
; different sample values to make sure that most MBPs are observed and that spurious 
; detections in the form of granules are not detected.s
;====================================================================================================================================  	
ENDWHILE	

objs = (objs)+(splst)
objs[WHERE(objs gt (0.))] =1.
;sigi=SQRT((MOMENT(grads))(1))
;sigi=SQRT((MOMENT(grads))(1)) *2   ;SIMS corrections 3
;sigi=SQRT((MOMENT(grads))(1)) *1.3   ;ROSA/Hinode
;sigi=SQRT((MOMENT(grads))(1)) *0.8    ;ROSA DATA
;sigi=SQRT((MOMENT(grads))(1)) *1.5   ;SST WB DATA TEST 1 
sigi=SQRT((MOMENT(grads))(1)) *2.0   ;GREGOR HiFI TESTs (continuum)

;thres = (sigfig(MEAN(grads)+(0.5*sigi),2)) ; SST H_ALPHA
;thres = (sigfig(MEDIAN(grads)-(0.5*sigi),2)) ; ROSA
;thres = (sigfig(MEAN(grads)-(sigi),2)) + ((sigfig(MEAN(grads)-(sigi),2))*0.04); ROSA ERRORS TEST
;thres = (sigfig(MEAN(grads)-(sigi),2)) ; G_BAND - simulations
;thres = (sigfig(MEAN(grads)-(0.75*sigi),2)) ; New SST WB Data
thres = (sigfig(MEDIAN(grads)-(0.75*sigi),2)) ; HiFi Gregor Test


;0.8 ROSA 1.5 Sims ;Tried 1.1 for rebinned sims, gives 30,000km average its the best
RETURN, thres

end
