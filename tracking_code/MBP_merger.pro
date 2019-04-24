;+
; ROUTINE:    	MBP_merger
;
; PURPOSE:    	CHECKS WHICH MBP WAS CLOSER TO A MERGED MBP TO TRACK IT
;		
; USEAGE:     	NEXT_FRAME
; INPUT:	FOLLWOING FRAME WHERE TWO MBPS HAVE MERGED
; OUTPUT:	CLOSEST MERGER TO FOLLOW
; AUTHOR:     	Version 2.1 
;

;AMENDMENT HISTORY
;(v2.1) - included merging numbers in the output table
;(v2.1) - if only one merged number then this means that MBP
;	  has been grown to touch a previously found MBP, thus
;	  code shall overwrite previous MBP with new information 
;	  on size and (x, y) position	
;+

PRO MBP_merger, file, mergd,separ,prev, xprev,yprev, x, y, imno

imsize= SIZE(mergd)
contin=-1

readcol,'table.dat',f='I,I,I,I,I,I,I',id,frm,xCoG,yCoG,ar,frmin,merg,/SIL
inst = INTARR(7, N_ELEMENTS(id))
inst[0,*] = id
inst[1,*] = frm
inst[2,*] = xCoG
inst[3,*] = yCoG
inst[4,*] = ar
inst[5,*] = frmin
inst[6,*] = merg

temp = (mergd)-(separ)
same = INTARR(imsize[1],imsize[2])
RofM = REGION_GROW(temp,WHERE(mergd eq (2.)))
same[RofM] = temp[RofM]

sm= SIZE(same, /Dimensions) &$      
mm = TOTAL(same) &$
xm = ROUND(TOTAL( TOTAL(same, 2) * Indgen(sm[0]) ) / mm) &$   
ym = ROUND(TOTAL( TOTAL(same, 1) * Indgen(sm[1]) ) / mm) &$

IF (xm) ge (x) THEN (xdis)=(xm)-(x) ELSE (xdis)=(x)-(xm) &$
IF (xm) ge (y) THEN (ydis)=(ym)-(y) ELSE (ydis)=(y)-(ym) &$
mgdgap =  (SQRT((xdis)^2+(ydis)^2))

mergd_nums = inst[0,(WHERE((inst[2,*]) eq (xm) AND (inst[3,*]) eq (ym) AND (inst[1,*]) eq (imno+1)))] 

IF (N_ELEMENTS(mergd_nums)) eq (1) THEN BEGIN &$
	mergd_nums = INTARR(2) &$
	ssep= SIZE(separ, /Dimensions) &$      
	msep = TOTAL(separ) &$
	xsep = ROUND(TOTAL( TOTAL(separ, 2) * Indgen(ssep[0]) ) / msep) &$   
	ysep = ROUND(TOTAL( TOTAL(separ, 1) * Indgen(ssep[1]) ) / msep) &$
	asep = N_ELEMENTS(WHERE(separ gt (0))) &$
	mergd_nums[0] = inst[0,(WHERE((inst[2,*]) eq (xm) AND (inst[3,*]) eq (ym)))]  &$
	mergd_nums[1] = inst[0, (N_ELEMENTS(inst[0,*])-1)] &$
	inst[2,WHERE(inst[0,*] eq (mergd_nums[0]))]=xsep &$
	inst[3,WHERE(inst[0,*] eq (mergd_nums[0]))]=ysep &$
	inst[4,WHERE(inst[0,*] eq (mergd_nums[0]))]=asep &$
ENDIF

for i = 0, (N_ELEMENTS(mergd_nums)-1) DO BEGIN &$
	IF ((WHERE(prev eq (mergd_nums[i])))((N_ELEMENTS(WHERE(prev eq mergd_nums[i])))-1)) gt (-1) THEN BEGIN &$
		x_new = xprev[((WHERE(prev eq mergd_nums[i]))((N_ELEMENTS(WHERE(prev eq mergd_nums[i])))-1))]  &$
		y_new = yprev[((WHERE(prev eq mergd_nums[i]))((N_ELEMENTS(WHERE(prev eq mergd_nums[i])))-1))]  &$ 
	ENDIF ELSE BEGIN &$
		x_new = inst[2,WHERE(inst[0,*] eq mergd_nums[i] AND (inst[1,*]) eq (imno))]  &$
		y_new = inst[3,WHERE(inst[0,*] eq mergd_nums[i] AND (inst[1,*]) eq (imno))]  &$ 
	ENDELSE &$

	IF (xm) ge (x_new) THEN (xdis) = (xm) - (x_new) ELSE (xdis) = (x_new) - (xm) &$
	IF (ym) ge (y_new) THEN (ydis) = (ym) - (y_new) ELSE (ydis) = (y_new) - (ym) &$
	IF (SQRT((xdis)^2+(ydis)^2)) gt (mgdgap) OR (contin gt (-1.) AND (SQRT((xdis)^2+(ydis)^2)) eq (mgdgap)) THEN (inst[2, WHERE(inst[0,*] eq (mergd_nums[i]) AND (inst[1,*]) eq (imno+1))])=(-1) &$
	IF (SQRT((xdis)^2+(ydis)^2)) lt (mgdgap) OR (contin eq (-1.) AND (SQRT((xdis)^2+(ydis)^2)) eq (mgdgap)) THEN (contin) = (mergd_nums[i]) &$
	IF (SQRT((xdis)^2+(ydis)^2)) lt (mgdgap) THEN (mgdgap) = (SQRT((xdis)^2+(ydis)^2)) &$
ENDFOR

inst[6,(WHERE(inst[0,*] eq (contin) AND (inst[1,*]) eq (imno+1)))]=inst[0,(WHERE(inst[2,*] eq (-1.)))]
inst = inst[*,WHERE(inst[2,*] gt (-1.))]

OPENW,1,'table.dat'
printf,1, inst
CLOSE,1 

end
