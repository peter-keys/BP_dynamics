;+
;
; ROUTINE:	MBP_count
; PURPOSE:	TO DETEMINE HOW MANY OBJECTS EXIST IN AN IMAGE
;		
; USAGE: 	ARRAY
; INPUT: 	ARRAY - THE IMAGE UNDER INVESTIGATION
; OUTPUT: 	CNT - THE NUMBER OF OBJECTS COUNTED
; CALLING:	CALLED IN MBP_tracking.pro & MBP_tabulation.pro
; AUTHOR: 	Version 2..0 
;

FUNCTION MBP_count, array

cnt = 0.
inbp = array
siz = SIZE(array)

WHILE (MAX(inbp)) gt (0) DO BEGIN
	multi = INTARR(siz[1],siz[2])
	ROI_multi = REGION_GROW(inbp,((WHERE(inbp gt 0.))(0)))
	multi[ROI_multi] = inbp[ROI_multi]
	inbp = (inbp) - (multi)
	cnt = (cnt)+(1.)	
ENDWHILE

RETURN, cnt

end

