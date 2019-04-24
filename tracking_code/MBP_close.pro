;+
;
; ROUTINE:      MBP_close
;
; PURPOSE:      DETERMINES WHICH OBJECT IS CLOSEST TO COFG
; USEAGE:       IMG, X,Y
; INPUT:        IMG - THE BINARY IMAGE TO FIND CLOSEST OBJECT FROM
;               X - THE X-COORDINAT TO BE CLOSEST TO
;               Y - THE Y-COORDINAT TO BE CLOSEST TO.  
; OUTPUT:      	RES - THE CLOSEST OBJECT
; CALLING:	CALLED WITHIN MBP_tracking.pro & MBP_tabulation.pro
; AUTHOR:       Version 2.0 
;

FUNCTION MBP_close, img, x, y

arrs = SIZE(img)
res = INTARR(arrs[1],arrs[2])

IF (x-49.) lt (0.) THEN (xdn) = 0. ELSE (xdn) = (x-49.) 
IF (x+49.) gt (arrs[1]-1.) THEN (xup) = (arrs[1]-1.) ELSE (xup) = (x+49.) 
IF (y-49.) lt (0.) THEN (ydn) = 0. ELSE (ydn) = (y-49.) 
IF (y+49.) gt (arrs[2]-1.) THEN (yup) = (arrs[2]-1.) ELSE (yup) = (y+49.) 

temp = img[xdn:xup,ydn:yup]
temps = SIZE(temp)
reobj = INTARR(temps[1],temps[2])
tmped = INTARR(temps[1],temps[2])
tmped[(1):(temps[1]-2),(1):(temps[2]-2)] = 1.

temp = (temp)*(tmped)

dist = 1000.
xnw = (x)-(xdn)
ynw = (y)-(ydn)

WHILE (MAX(temp)) gt (0.) DO BEGIN &$
	newbp = INTARR(temps[1],temps[2])
	regcls= REGION_GROW(temp,((WHERE(temp eq (1)))(0))) &$
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

res[xdn:xup,ydn:yup] = (res[xdn:xup,ydn:yup])+(reobj)

RETURN, res

end 
