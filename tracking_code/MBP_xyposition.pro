;+
; ROUTINE:    xyposition
;
; PURPOSE:    Prints the x & y coords, given an index
;
; USEAGE:     printxy,Array, Index, x,y
; INPUT:
;
; OUTPUT:    
; 
; Example:

; AUTHOR:     
;
;-

FUNCTION MBP_xyposition, array, flag

   xy = INTARR(2,N_ELEMENTS(flag));sets up an integer array with 2 x (number of elements that are in the range set)
   S = SIZE(array);gives the size information of the time series array
   
   FOR i = 0.,(N_ELEMENTS(flag)-1) DO BEGIN ;loops around all the elements in the flag []
       y = FIX(flag[i] / S(1))
       x = flag[i] - (y * S(1))
       
       ;y = FLOAT(FIX(FLOAT(flag[i])/S(1)));determines the y-coord
       ;x = FLOAT(flag[i] - (y * S(1))); determines the x-coord
       xy[0,i] = x ;sets column 0 of xy[] to value of x coord
       xy[1,i] = y ;sets column 1 of xy[] to value of y-coord
   ENDFOR

RETURN,xy ; returns the xy[] to the function

END
