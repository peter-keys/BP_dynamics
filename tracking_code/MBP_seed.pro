;+
;
; ROUTINE:	MBP_seed
;
; PURPOSE:    	TO DETERMINE A SEED AREA FOR MBPS THAT HAVE BEEN MISSED
;		
; USEAGE:   	OB,GRP,X,Y	  	
; INPUT:	OBJ - ORIGINAL OBJECT TO BE TRACKED
;		SPLGRP - THE GROUP OF POSSIBLE SEEDS
;		X - THE X-COORDINATE OF THE ORIGINAL OBJECT
;		Y - THE Y-COORDINATE OF THE ORIGINAL OBJECT
;		si - SIZE OF IMAGE
; OUTPUT:	POSS - RETURNS THE POSSIBLE SEED REGION
; CALLING	CALLED FROM MBP_tracking.pro & MBP_tabulating [line 132]
; AUTHOR:     	Version 1.0 
;
;AMENDMENT HISTORY
; (v2.1)-	Higher initial intensity threshold employed to set a seed region.  Now potential seed regions must
;		     overlap original object by 1/2 the number of pixels prior to being separated.  This stops the inclusion
;		     of single pixel seeds.
; (v2.1)-	Adjustable growing threshold to ensure every chance is given to grow a missed MBP.
; (v2.1)-	Grown array is now set to zero during varying stage to prevent erronous inclusion of small MBPs

FUNCTION MBP_seed, ob, grp, x, y, si

poss = INTARR(si[1], si[2])
count = 0.

IF (MAX(ob+grp)) eq (1.) THEN BEGIN
	IF (x-49.) lt (0.) THEN (xdn) = 0. ELSE (xdn) = (x-49.) 
	IF (x+49.) gt (si[1]-1.) THEN (xup) = (si[1]-1.) ELSE (xup) = (x+49.) 
	IF (y-49.) lt (0.) THEN (ydn) = 0. ELSE (ydn) = (y-49.) 
	IF (y+49.) gt (si[2]-1.) THEN (yup) = (si[2]-1.) ELSE (yup) = (y+49.) 

	temp = grp[xdn:xup,ydn:yup]
	temps = SIZE(temp)
	reobj = INTARR(temps[1],temps[2])
	tmped = INTARR(temps[1],temps[2])
	tmped[(1):(temps[1]-2),(1):(temps[2]-2)] = 1.

	temp = (temp)*(tmped)
	
	dist = 5.
	xnw = (x)-(xdn)
	ynw = (y)-(ydn)

	WHILE (MAX(temp)) gt (0.) DO BEGIN &$
		newbp = INTARR(temps[1],temps[2]) &$
		regcls =REGION_GROW(temp,((WHERE(temp eq (1)))(0))) &$
		newbp[regcls] = temp[regcls]  &$    

		temp = (temp)-(newbp) &$

		bps = SIZE(newbp,/DIMENSIONS) &$ 
		bpm = TOTAL(newbp) &$ 
		xbp = ROUND(TOTAL(TOTAL(newbp,2)*Indgen(bps[0]))/bpm) &$ 
		ybp = ROUND(TOTAL(TOTAL(newbp,1)*Indgen(bps[1]))/bpm) &$ 

		IF (xbp) ge (xnw) THEN (x_dist) = (xbp) - (xnw) ELSE (x_dist) = (xnw) - (xbp) &$
		IF (ybp) ge (ynw) THEN (y_dist) = (ybp) - (ynw) ELSE (y_dist) = (ynw) - (ybp) &$
		IF (SQRT((x_dist)^2+(y_dist)^2)) lt (dist) THEN (reobj) = newbp &$
		IF (SQRT((x_dist)^2+(y_dist)^2)) lt (dist) THEN (dist) = (SQRT((x_dist)^2+(y_dist)^2)) &$
	ENDWHILE
	IF (MAX(reobj)) gt (0) AND (N_ELEMENTS(WHERE(reobj gt 0))) ge (4) THEN (poss[xdn:xup,ydn:yup]) = (poss[xdn:xup,ydn:yup])+(reobj)
ENDIF

IF (MAX(ob+grp)) gt (1.) THEN BEGIN
	WHILE (MAX(ob+grp)) gt (1.) DO BEGIN &$
		bud = INTARR(si[1], si[2]) &$
		RegS = REGION_GROW(grp, ((WHERE((ob+grp) gt (1)))(0))) &$ 
		bud[RegS] = grp[RegS] &$
		
		IF (N_ELEMENTS(WHERE(bud gt 0))) ge (4) THEN BEGIN &$
			buds = SIZE(bud,/DIMENSIONS) &$ 
			budm = TOTAL(bud) &$ 
			xbud = ROUND(TOTAL(TOTAL(bud,2)*Indgen(buds[0]))/budm) &$ 
			ybud = ROUND(TOTAL(TOTAL(bud,1)*Indgen(buds[1]))/budm) &$ 

			IF (xbud) ge (x) THEN (x_dist) = (xbud) - (x) ELSE (x_dist) = (x) - (xbud) &$
			IF (ybud) ge (y) THEN (y_dist) = (ybud) - (y) ELSE (y_dist) = (y) - (ybud) &$

			IF (count) eq (0.) THEN BEGIN &$
				(poss) = (bud) &$ 
				dist = (SQRT((x_dist)^2+(y_dist)^2)) &$
			ENDIF ELSE BEGIN &$
				temp = poss &$
				poss = INTARR(si[1], si[2], count+1) &$
				IF (SQRT((x_dist)^2+(y_dist)^2)) lt (dist) THEN (inlo) = (1.) ELSE (inlo) = (0.) &$ 
				IF (SQRT((x_dist)^2+(y_dist)^2)) lt (dist) THEN (inhi) = (count) ELSE (inhi) = (count-1.) &$
				poss[*,*,inlo:inhi] = temp &$
				IF (SQRT((x_dist)^2+(y_dist)^2)) lt (dist) THEN (inset) = (0.) ELSE (inset) = (count) &$ 
				poss[*,*,inset] = bud &$
			ENDELSE &$
			(count)=(count)+(1) &$
		ENDIF &$
		grp = (grp) - (bud) &$
	ENDWHILE
ENDIF

RETURN, poss

end


	
