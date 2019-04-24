;+
;
; ROUTINE:      MBP_gradient
;
; PURPOSE:      FIND THE FIRST DERIVATIVE OF A ONE DIMENSIONAL ARRAY
; USEAGE:       ARR
; INPUT:	ARR - ARRAY TO CALCULATE FIRST DERIVATIVE FROM                  
; OUTPUT:       FIRDER - THE FIRST DERIVATIVE OF THE ARRAY
; CALLING:	CALLED IN MBP_detect, MBP_tracking & MBP_tabulation
; AUTHOR:       Version 2.0 
;

FUNCTION MBP_gradient,arr
sz = size(arr)
if sz(0) ne 1 then goto, finishup
firder = arr(1:sz(1)-1)-arr(0:sz(1)-2)
return,firder
finishup:
end
