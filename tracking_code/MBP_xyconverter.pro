;+
; ROUTINE:    xyposition
;
; PURPOSE:    Prints the x & y coords, given an index
;
; USEAGE:     flags = xypositionconverter(values, array)
; INPUT:      values -> array of x & y which need converted to 1D image values
;   	      array -> 2D image from which the values are taken
;
; OUTPUT:    
; 
; Example:

; AUTHOR:     D. B. Jess
;
;-

FUNCTION MBP_xyconverter, array, values

   S = SIZE(array)
   flag = SIZE(values)
   IF N_ELEMENTS(flag) lt 5. THEN flag = 1. ; STOPS SINGLE ELEMENTS ARRAYS FROM BEING USED WRONGLY
   IF N_ELEMENTS(flag) ge 5. THEN flag = flag[2] ; ESTABLISHES HOW MANY PAIRS OF XY COORDINATES THEIR ARE
   flag = DBLARR(flag)
   FOR i = 0.,(N_ELEMENTS(flag)-1) DO flag[i] = (((values[1,i] * S[1]) + values[0,i]))
   
RETURN,flag

END
