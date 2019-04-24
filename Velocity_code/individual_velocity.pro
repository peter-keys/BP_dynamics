PRO individual_velocity, BP

;READCOL,'test3_tracking_file_sigma1',f='D,D,D,D,D,x,x',BP_num,frame,x,y,area
READCOL,'/data/rosa3/oldrosa1/phk/data/28May2009/gband/tracking_files/obs_track_new.dat',f='D,D,D,D,D,x,x',BP_num,frame,x,y,area

frame_num = frame[WHERE (BP_num eq BP)]
first_frame = MIN(frame_num)
last_frame = MAX(frame_num)
xcoords = x[WHERE (BP_num eq BP)]  
ycoords = y[WHERE (BP_num eq BP)]  
distance = 25
time = 2.1

print,'Frame		x		y		Velocity'
vel = DBLARR(N_ELEMENTS(frame_num))
FOR i = first_frame, last_frame DO BEGIN &$
	FOR j = 1, (N_ELEMENTS(frame_num)-1) DO BEGIN &$
	  vel[0] = 0 &$
	  d = (SQRT((((xcoords[j]-xcoords[j-1])*distance)^2)+(((ycoords[j]-ycoords[j-1])*distance)^2))) &$
	  IF (d LT (distance * 2)) THEN velocity = 0. &$
  	  IF (d GE (distance * 2)) THEN velocity = ((d - (distance * 2)) / (time * (frame_num[j] - frame_num[j-1]))) &$
	  vel[j] = velocity &$

	ENDFOR &$
	print,arr2str(i,/trim)+'	'+arr2str(xcoords[i-first_frame],/trim)+'	'+arr2str(ycoords[i-first_frame],/trim)+'	'+arr2str(vel[i-first_frame],/trim) &$

ENDFOR
print, '  '
print,'Average velocity = '+arr2str(MEAN(vel),/trim)
print,'Number of frames = '+arr2str(N_ELEMENTS(frame_num),/trim)

END
