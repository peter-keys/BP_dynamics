;+
;
; ROUTINE: 	SPLIT_PLAYABOUT
;
; PURPOSE: 	PLAYABOUT SPLITTING THE MBPS
; USEAGE:  	BTPT_GROUP
; INPUT:   	BTPT_GROUP - GROUP TO BE SPLIT
;		  	
; OUTPUT:  	SPLIT_BTPTS
; AUTHOR:  	Version 1.0 - Philip. J. Crockett, QUB, 12 May 2008
;   	    	         (Email: pcrockett02@qub.ac.uk)

FUNCTION MBP_split,group,im_size,btpt_sigma
sto = group
objspl = INTARR(im_size[1], im_size[2])
group[WHERE(group lt ((MIN(group[WHERE(group gt 0.)]))+(0.5*btpt_sigma)))] = 0.

WHILE (MAX(group)) gt (0.) DO BEGIN &$
	binary = INTARR(im_size[1], im_size[2]) &$
	binary[WHERE(group gt (0.))] = 1. &$
	reg = REGION_GROW(binary,((WHERE(binary gt 0))(0)))  &$
	
	splobj = FLTARR(im_size[1], im_size[2]) &$ 
	splobj[reg]=group[reg] &$
	group = (group) - (splobj) &$
		
	IF (((MAX(splobj))-(MIN(splobj[WHERE(splobj gt 0.)])))/btpt_sigma) le (3.) THEN BEGIN &$
		IF (N_ELEMENTS(WHERE(splobj gt (0.)))) ge (4.) THEN (objspl[WHERE(splobj gt (0.))]) = 1. &$
	ENDIF ELSE BEGIN &$ 
		splobj[WHERE(splobj lt ((MIN(splobj[WHERE(splobj gt 0.)]))+(0.5*btpt_sigma)))] = 0. &$
		group = (group)+(splobj) &$
	ENDELSE	 &$ 
ENDWHILE

IF (MAX(objspl)) eq (0.) THEN (objspl[WHERE(sto gt ((MAX(sto))-(3*btpt_sigma)))]) = 1.

RETURN, objspl

end
