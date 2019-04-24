;+
;ROUTINE:	Zeeman_split_estimate
;
;PURPOSE:	A way of calculating the Zeeman Spitting B-field for Stokes V plots as a check  
;
;USEAGE:	result =ZEEMAN_SPLIT_ESTIMATE(profiles, grids, nicole_b)
;
;INPUT:		profiles = The Stokes V profiles you want to make an estimate on
;		grid = The wavelength grid of your interpolated Stokes V profiles
;		nicole_b = The LOSB evolution estimates (from Becky) found with NICOLE for comparison later
;
;OUTPUT:	Array of values for estimates for the MBP (for 6301, 6302, and the actual Nicole value
;		for comparison...)
;		
;AUTHOR:	Peter H. Keys, August 2011
;
;-

FUNCTION zeeman_split_estimate,profiles,grid,nicole_b



sizing = SIZE(profiles,/dimensions)
IF N_ELEMENTS(sizing) eq 2. THEN frames = sizing[1]

!p.background=255.
!p.color=0.

window,0,title='Zeeman Split estimations',xsize=1000,ysize=800
mult,1,1

print,' '
print,'------------------------------------'
Print,' You will make estimates of '
print,' Zeeman splitting based on '
print,' Stokes V profiles. This requires'
print,' you to click the peaks in each to  '
print,' work out the split.' 
print,' '
print,' Start with 6301. Then cycle through'
print,' again for 6302 profiles.....'
print,'------------------------------------'
print,' '

B_fields = DBLARR(3, frames)

FOR i = 0, frames - 1 DO BEGIN &$
	loadct,0,/SILENT &$
	plot,grid,profiles[*,i], xtitle='Wavelength (A)',ytitle='Stokes V',Title='Click Peaks, Left to Right, Both 6301 and 6302',subtitle='Frame '+arr2str(i,/trim)+' of '+arr2str(frames-1,/trim),thick=2,charsize=1.4
	print,'Pick 6301 peak 1 & 2 '&$
	peak_1_6301 = pick_elements() &$
	tek_color
	verline,peak_1_6301[0],line=2,color=2 &$
	peak_2_6301 = pick_elements() &$
	verline,peak_2_6301[0],line=2,color=2 &$
	diff_6301 = (ABS(peak_1_6301[0]-peak_2_6301[0]))/2. &$ 		;Sigma value (i.e., average distance from non-split) 
	B_6301 = diff_6301 / ((4.67E-13) * 1.67 * ((6301.5)^2.)) &$ 	;B field from Zeeman (Spruit 1981?79?) for 6301 (g = 1.67)
	
	print,' ' &$
	print,'Pick 6302 peak 1 & 2 '&$
	peak_1_6302 = pick_elements() &$
	verline,peak_1_6302[0],line=2,color=2 &
	peak_2_6302 = pick_elements() &$
	verline,peak_2_6302[0],line=2,color=2 &$	
	print,' ' &$
	diff_6302 = (ABS(peak_1_6302[0]-peak_2_6302[0]))/2. &$ 		;Sigma value (i.e., average distance from non-split) 
	B_6302 = diff_6301 / ((4.67E-13) * 2.5 * ((6302.5)^2.)) &$ 	;B field from Zeeman (Spruit 1981?79?) for 6302 (g = 1.67)	
	
	;NB B field is in Gauss for this calculation....
	
	B_fields[*,i] = [B_6301,B_6302,ABS(nicole_b[i])] &$	
ENDFOR
 

RETURN, B_fields

END

