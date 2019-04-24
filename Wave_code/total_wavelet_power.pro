FUNCTION total_wavelet_power, datacube, cadence, first_period_in_sec, last_period_in_sec, sub_oct

lightcurve = REFORM(datacube[0,0,*])
wave = WAVELET(lightcurve,cadence,PERIOD=period,SCALE=scale,COI=coi,/PAD,$
       SIGNIF=signif,siglvl=0.99,FFT_THEOR=fft_theor,mother='morlet',DJ=1./sub_oct)
elements = WHERE(period ge first_period_in_sec AND period le last_period_in_sec)
start_period_element = MIN(elements)
end_period_element = MAX(elements)
num_periods = end_period_element - start_period_element + 1

a=0
b=0
dt = cadence/60.        ; time in minutes, of cadence
pad = 1
s0 = 2*dt   
dj = 1/sub_oct ; this will do 16 sub-octaves per octave
j1 = 1./dj  ; for the period range
mother = 'Morlet'
sl = string(byte([27, 91, 68, 27, 91, 68, 27, 91, 68]))
countup = 0.
t0 = SYSTIME(1) 
t1 = t0

xsize = N_ELEMENTS(datacube[*,0,0])
ysize = N_ELEMENTS(datacube[0,*,0])
zsize = N_ELEMENTS(datacube[0,0,*])

timelength = zsize * dt
print,''
print,'Time series lasts for :',timelength,'min'
print,''
wait,5.

lightcurve = REFORM(datacube[0,0,*])
;lightcurve = (lightcurve - AVERAGE(lightcurve)) 
time = FINDGEN(zsize) * dt
xrange = [0,timelength]
recon_lightcurve = wave_recon(lightcurve,cadence)
wave = WAVELET(recon_lightcurve,cadence,PERIOD=period,SCALE=scale,COI=coi,/PAD,$
       SIGNIF=signif,siglvl=0.99,FFT_THEOR=fft_theor,mother=mother,DJ=dj)
all_periods = period
save_period = period[start_period_element:end_period_element]
       

IF num_periods gt (N_ELEMENTS(all_periods)) THEN num_periods = (N_ELEMENTS(all_periods))
new_datacube = FLTARR(xsize,ysize,num_periods)

print,''
print,'Power cube will be : '+ARR2STR(xsize,/trim)+'x by '+ARR2STR(ysize,/trim)+'y by '+ARR2STR(num_periods,/trim)+' periods'
print,''

image_ave = TOTAL(datacube,3) / zsize
LOADCT,0,/silent

window,0,title='Wavelet analysis',xsize=800,ysize=400
mult,2,1

lowest_power = 1.e22
FOR a = 0,(xsize-1) DO BEGIN
    FOR b= 0, (ysize-1) DO BEGIN
        lightcurve = REFORM(datacube[a,b,*])
	IF MAX(lightcurve) ne MIN(lightcurve) THEN BEGIN 
	;lightcurve = (lightcurve - AVERAGE(lightcurve)) 
        time = FINDGEN(zsize) * dt
	xrange = [0,timelength]
	recon_lightcurve = wave_recon(lightcurve,cadence)
        wave = WAVELET(recon_lightcurve,cadence,PERIOD=period,SCALE=scale,COI=coi,/PAD,$
    	       SIGNIF=signif,siglvl=0.99,FFT_THEOR=fft_theor,mother=mother,DJ=dj)
        power = (ABS(wave))^2
        global_ws = TOTAL(power,1)/zsize
        J = N_ELEMENTS(scale) - 1
        pos1 = [0.12,0.55,0.95,0.9]
        ;PLOT,time,recon_lightcurve,xst=1,yst=1,charsize=2, $
        ;     YTITLE='intensity',XTITLE='Time (min)',POSITION=pos1,title='Cadence of data = '+strtrim(cadence,1)+'s'
	yrange = [period(start_period_element),period(end_period_element)]
        p=max(power)
	levels = [0.3*p,0.5*p, 0.7*p, 0.9*p]
	IF levels[3] eq levels[0] THEN levels = [0, 0.25, 0.5, 0.75]
	colors = [200,150,100,50]
        period2 = FIX(ALOG(period)/ALOG(2))
        ytickv = 2.^(period2(UNIQ(period2)))
        pos2 = [pos1(0),0.12,pos1(2),0.45]
        ;CONTOUR,power^1.0,time,period,/NOERASE,POSITION=pos2,charsize=2, $
        ;        YRANGE=yrange,xst=1, $ 
        ;       LEVELS=levels,C_COLORS=colors,/FILL, $
        ;        XTITLE='Time (min)',YTITLE='Period (sec)'
    	;PLOTS,time,coi,NOCLIP=0
	power_period = REFORM(power[*,start_period_element:end_period_element]) ; power in form POWER[lightcurve element , period]
    	new_datacube[a,b,*] = REFORM(TOTAL(power_period,1)/zsize)
	IF MIN(new_datacube[a,b,*]) lt lowest_power THEN lowest_power = MIN(new_datacube[a,b,*])
        ENDIF
	countup = countup + 1.
	systm = SYSTIME(1)
    	loop_time = (systm - t1)
    	t1 = systm
	minutes = (FIX(((loop_time)*((xsize*ysize) - countup))/60.))
	hours = FIX(minutes/60.)
	hourminutes = hours*60.
	remain_minutes = minutes - hourminutes
	IF b MOD 200 eq 0 THEN tvim,image_ave>0.,title='Average Image',pcharsize=2
	IF b MOD 200 eq 0 THEN tvim,(TOTAL(new_datacube,3)/num_periods)>lowest_power,title='Total Power',pcharsize=2
	IF b MOD 100 eq 0 THEN print,'I am '+ARR2STR(((countup/(xsize*ysize))*100.),/trim) + $
	                             ' percent complete with approximately '+ARR2STR(hours,/trim)+' hours and '+$
				     ARR2STR(FIX(remain_minutes),/trim)+' minutes left at current CPU speed'
    ENDFOR
ENDFOR

mult,1,1

print,'The periods (in seconds) which have been computed correspond to the values of:'
FOR i = 0,(num_periods-1) DO print,i,': ',save_period[i]
SAVE,FILENAME='wavelet_periods_'+ARR2STR(FIX(first_period_in_sec),/trim)+'s_'+ARR2STR(FIX(last_period_in_sec),/trim)+'s.sav',save_period
print,''
print,'These periods have been saved into the file wavelet_periods_'+ARR2STR(FIX(first_period_in_sec),/trim)+'s_'+ARR2STR(FIX(last_period_in_sec),/trim)+'s.sav'

RETURN, new_datacube

END
