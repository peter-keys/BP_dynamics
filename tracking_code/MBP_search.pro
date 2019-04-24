;+
;
; ROUTINE:    	MBP_search
;
; PURPOSE:    	SETS-UP AND EXAMINES A SCAN AREA TO SEARCH FOR AN ASSOCIATED MBP 
; USEAGE:     	NEXTIM, LOX, HIX, LOY, HIY,X,Y
; INPUT:	NEXTIM - THIS IS THE ARRAY OF DETECTED MBPS IN THE FRAME UNDER
;			 INVESTIGATION
;		LOX - THE LOWER X BOUNDARY OF THE SEARCH AREA
;		HIX - THE UPPPER X BOUNDARY OF THE SEARCH AREA
;		LOY - THE LOWER Y BOUNDARY OF THE SEARCH AREA
;		HIY - THE UPPER Y BOUNDARY OF THE SEARCH AREA
;		X - THIS IS THE CENTER OF GRAVITY X-COORD OF THE PREVIOUS MBP
;		Y - THIS IS THE CENTER OF GRAVITY Y_COORD OF THE PREVIOUS MBP
; OUTPUT:	CONFIRM - ARRAY TO STORE THE CONFIRMED ASSOCATIED MBP
; CALLING:	CALLED IN MBP_tracking & MBP_tabulating
; AUTHOR:     	Version 3.0 
;
; AMENDMENT HISTORY
; (v2) - Annotated detail was added
; (v3) - Full rewrite of the code.  Thought of a new method to detect the
;	       closest associted object and wrote version 3

FUNCTION MBP_search, nextim, lox, hix, loy, hiy, x, y
;----------------------------------------------------------------------------------------------------------------------------
;SECTION ONE: 	LINE 35 TO 38
;PURPOSE: 	To setup some intial parameters
;PARAMETERS: 	SNX - the [x,y] dimensitons of the array under consideration
;		AREA - the area under investigation
;		BTPT_CACHE - array to store up to 10 objects that exist within the
;			     search area
;		BTPT_STORE - stores the objects so that they are only recorded once
;			     during the scanning process		
;DESCRIPTION:	
;Set up the scan area across which to search for associated objects and initilises a
;distance parameter so that the closest object is chosen as the associated object.
;----------------------------------------------------------------------------------------------------------------------------
snx = SIZE(nextim)
area = INTARR(snx[1],snx[2])
confirm =INTARR(snx[1], snx[2])	
area[lox:hix,loy:hiy] = 1.
dist = 1000
;----------------------------------------------------------------------------------------------------------------------------
;SECTION TWO: 	LINES 55 TO 81
;PURPOSE: 	To scan the area and investigate each MBP that exists within it
;		and selecting an associate MBP
;PARAMETERS: 	CANDS - stores all possible candidates that exist in the scan area 
;		CAND - stores each individual candidate for investigation
;		COFG - stores the centre of gravity for each candidate
;DESCRIPTION:	
;Secion 2 scans the area for all possible candidates and investigates them
;individually.  Each candidate's centre of gravity is calculated.  The object must pass
;three stipulations to be considered an assocatied object.  Firstly, an objects center
;of gravity must exist within the search area.  Secondly, the center of gravity of the associated
;object must not have moved more than 5 pixels.  Thirdly, the object must be the closet
;object to the original MBP.  If all of these stipulations are met, then the object is
;deemed the associated MBP.
;----------------------------------------------------------------------------------------------------------------------------
IF (MAX(area+nextim)) gt (1.) THEN BEGIN	
	cands = INTARR(snx[1],snx[2])
	RofBP = REGION_GROW(nextim,WHERE((area+nextim) gt (1.)))
	cands[RofBP] = nextim[RofBP]
	
	WHILE (MAX(cands)) gt (0.) DO BEGIN &$
		cand = INTARR(snx[1], snx[2]) &$
		region = REGION_GROW(cands,((WHERE(cands gt 0))(0))) &$
		cand[region] = cands[region] &$
		
		scand = SIZE(cand, /Dimensions) &$        
		mcand = TOTAL(cand) &$         			
		xcand = (TOTAL(TOTAL(cand,2)*Indgen(scand[0]))/mcand) &$    
		ycand = (TOTAL(TOTAL(cand,1)*Indgen(scand[1]))/mcand) &$
		
		CofG = INTARR(snx[1],snx[2]) &$
		CofG[xcand,ycand] = 1. &$
		
		cands = (cands)-(cand) &$
		
		IF (x) gt (xcand) THEN (xdist)= (x)-(xcand) ELSE (xdist)=(xcand)-(x) &$
		IF (y) gt (ycand) THEN (ydist)= (y)-(ycand) ELSE (ydist)=(ycand)-(y) &$
		dist_2old = SQRT(((xdist)^2)+((ydist)^2)) &$

		IF (MAX(CofG+area)) eq (2.) AND (dist_2old) le (5.) AND (dist_2old) lt (dist) THEN (confirm) = cand  &$ 
		IF (dist_2old) lt (dist) THEN (dist) = (dist_2old) &$
	ENDWHILE
;----------------------------------------------------------------------------------------------------------------------------
;SECTION THREE:	LINES ** TO **
;PURPOSE: 	If no associated MBP is confirmed, to return a null array
;PARAMETERS: 	CONFIRM - a null array indicating a failure to detect an associated
;DESCRIPTION:	
;Returns a null area if no associated object is found.
;----------------------------------------------------------------------------------------------------------------------------
ENDIF ELSE (confirm) = INTARR(snx[1], snx[2])

RETURN, confirm

end
