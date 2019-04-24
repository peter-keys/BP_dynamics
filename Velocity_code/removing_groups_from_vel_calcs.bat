readcol,'/data/rosa3/oldrosa1/phk/data/28May2009/gband/tracking_files/obs_track_velocity_table',f='x,x,D,D,x,D,D,x,D',ix,iy,xf,yf,BP_velocity,SKIPLINE=1
readcol,'/data/rosa3/oldrosa1/phk/data/Hinode/Hinode_tracks_file_final_velocities',f='x,x,D,D,x,D,D,D,D',ix,iy,xf,yf,lifetime,BP_velocity,SKIPLINE=1

print,N_ELEMENTS(ix)

print,N_ELEMENTS(WHERE(ix GT 375 AND ix LT 420))
print,N_ELEMENTS(WHERE(xf GT 375 AND xf LT 420))

print,N_ELEMENTS(WHERE(iy GT 40 AND iy LT 85))
print,N_ELEMENTS(WHERE(yf GT 40 AND yf LT 85))

;Remove number of elements that pass into the box from final vel array

vel_array = DBLARR(3056-466-526-457-367)
vel_array = DBLARR(N_ELEMENTS(ix))
lifetime_array = DBLARR(N_ELEMENTS(ix))


FOR i = 0, (N_ELEMENTS(ix) - 1) DO BEGIN &$
	IF ix[i] GT 530 AND ix[i] LT 615 THEN CONTINUE &$
	IF xf[i] GT 530 AND xf[i] LT 615 THEN CONTINUE &$
	IF iy[i] GT 165 AND iy[i] LT 205 THEN CONTINUE &$
	IF yf[i] GT 165 AND yf[i] LT 205 THEN CONTINUE &$
	
	IF ix[i] GT 755 AND ix[i] LT 815 THEN CONTINUE &$
	IF xf[i] GT 755 AND xf[i] LT 815 THEN CONTINUE &$
	IF iy[i] GT 200 AND iy[i] LT 250 THEN CONTINUE &$
	IF yf[i] GT 200 AND yf[i] LT 250 THEN CONTINUE &$
	
	IF ix[i] GT 350 AND ix[i] LT 415 THEN CONTINUE &$
	IF xf[i] GT 350 AND xf[i] LT 415 THEN CONTINUE &$
	IF iy[i] GT 255 AND iy[i] LT 315 THEN CONTINUE &$
	IF yf[i] GT 255 AND yf[i] LT 315 THEN CONTINUE &$
	
	IF ix[i] GT 375 AND ix[i] LT 420 THEN CONTINUE &$
	IF xf[i] GT 375 AND xf[i] LT 420 THEN CONTINUE &$
	IF iy[i] GT 40 AND iy[i] LT 85 THEN CONTINUE &$
	IF yf[i] GT 40 AND yf[i] LT 85 THEN CONTINUE &$
	
	lifetime_array[i] = lifetime[i] &$
	vel_array[i] = BP_velocity[i] &$
ENDFOR

print,MEAN(lifetime_array[WHERE(lifetime_array GT 0)])
print,MEAN(vel_array[WHERE(vel_array GT 0)])

print,(MEANABSDEV(lifetime_array[WHERE(lifetime_array GT 0)],/MEDIAN)/2)/60
print,MEANABSDEV(vel_array[WHERE(vel_array GT 0)])/2
