;+
;NAME:
;    	SDBWAVE
;
;PURPOSE:
;   	Calculate and display the wavelet (power) transform of a time series
;   	(CURVE) with the time series above the plot for comparison with the
;   	transform, and also show the global wavelet spectrum to the right,
;   	since it represents the summation of the power transform over time,
;   	and is therefore analogous to the Fourier transform of the series.
;
;CALLING SEQUENCE:
;   	sdbwave,curve,print=print,pc=pc,title=title,nocon=nocon,delt=delt,$
;    	    mother=mother,offset=offset
;
;
;INPUTS:
;   	CURVE - a 1-D time series, e.g. a light curve
;
;INTERACTIVE INPUTS DURING EXECUTION:
;
;   	The user is first asked for a maximum value for the low-pass filter
;   	called by WAVE_EXAMPLE. For the SECIS data, an unfiltered time series
;   	will be analysed by specifying half the data sampling rate of 44.2 Hz
;   	(i.e. 22.1 Hz).
;
;   	The second parameter requested is the significance (confidence)
;   	level which the WAVELET function will calculate as a set of contours
;   	in the power transform. Concentrations of power enclosed within these
;   	contours will be significant to that level. E.g., if the user specifies
;   	a value of 0.99 (i.e. 99%) then there is a 1% chance that the power
;   	within the contours is due to noise.
;
;
;KEYWORDS:
;   	PRINT - if set, then set parameters applicable to printing to a
;   	    	postscript file. NOTE: Does not actually set the display
;   	    	to the PostScript device - this must be done manually.
;
;   	PC - the plot colour for the contours and cone-of-influence
;   	        boundary and cross-hatching. This is a value from 0 to
;   	    	255. In a monochromatic colour table (e.g. loadct,3 for
;   	    	red temperature) 0 will be black, 255 will be white, and
;   	    	any other colours will be somewhere in between (e.g.
;   	    	scale linearly for colour tables 0,1,3, & 8;
;   	    	logarithmically for colour table 9).
;
;   	TITLE - a string, specifying the title over the wavelet transform
;   	    	part of the display output.
;
;   	NOCON - if set, then do not display either contours or the
;   	    	cone-of-influence calculated by the WAVELET routine
;
;OUTPUTS:
;   	None, except for the display of the wavelet transform etc. to
;   	either X windows (the default) or (optionally) the PostScript device.
;
;
;EXAMPLES:
;   	IDL> wave4,series,title='Wavelet Transform for our series thing'
;
;   	 to write the output to an Encapsulated PostScript file:
;   	IDL> set_plot,'PS
;   	IDL> DEVICE,/encapsul,/color,filename='test1.eps',xsize=5,ysize=5
;   	IDL> sdbwave,series,title='Another Title',/print,pc=0
;   	IDL> DEVICE,/close
;   	IDL> SET_PLOT,'X'
;
;   	 or, if you have TOGGLE.PRO in your path:
;   	IDL> TOGGLE,filename='wavelet.ps',/landscape,/color
;   	IDL> sdbwave,series,title='Another Title',/print,pc=0
;   	IDL> TOGGLE
;
;CALLS:
;   	ECLIPSFILTER.PRO, WAVELET.PRO, WAVE_SIGNIF.PRO
;
;HISTORY:
;   	Based on ECLIPSE.PRO by Eoghan O'Shea, QUB, 2000
;   	Modified in stages by D. Williams for the SECIS eclipse data
;   	    "    "    "    "  S. Bloomfield for Sac. Peak 2002 data
;-

PRO sdbwave,curve,print=print,pc=pc,title=title,nocon=nocon,delt=delt,$
    mother=mother,offset=offset,cone=cone,fast=fast,outscale=outscale,$
    outper=outper,outfreq=outfreq,outtime=outtime,color=color,fig=fig,$
    short=short,long=long,ylog=ylog,lthick=lthick,sig_fudge=sig_fudge,$
    maxper=maxper;xticksep=xticksep

;check the number of elements in the input time series
elem=N_ELEMENTS(curve)
print,'No. of data-points',elem
;set the ever-useful angstrom plotting symbol in IDL character specification terms
ang = '!3!sA!r!u!9 %!3!n'

;set the yrange for plotting the timeseries in the first window and the upper panel
min1=MIN(curve)
max1=MAX(curve)

IF NOT KEYWORD_SET(pc) then pc=255

;change the titles from here
IF NOT KEYWORD_SET(title) THEN title = 'a) Time series'
mainTitle='Wavelet analysis and global power spectrum'	;overall title for the main wavelet window

;time interval betweend data points (in seconds), i.e. the cadence
IF NOT KEYWORD_SET(delt) THEN delt=(1*10.)
IF NOT KEYWORD_SET(sig_fudge) THEN sig_fudge=1.

;set the maximum frequency to be plotted, NOT the cut-off frequency for the low-pass filter:
topF=1./(2*delt)     	; the Nyquist frequency

;need to be able to work out the maximum of the COI for the
;purposes of determining the maximum trustable period
middle=elem/2
middle=FIX(middle)

;set up the time array. This is another example of customising in the actual routine,
;rather than being generic enough to be set at the command line.
IF NOT KEYWORD_SET(offset) THEN OFFSET=0.
time=(findgen(elem)*delt)+offset;+(39.46)
print,'Min. time',min(time)
print,'Max. time',max(time)

as=curve  	;avoid negative values by taking the modulus of the values only.
series=as  	;a duplicate of the filtered curve
;why Eoghan ever called it series is anyone's guess...

;for control over the size of the plot text. Should be 1 by default, but another routine
;might have set it to a different value.
!p.charsize=1.0
IF NOT KEYWORD_SET(lthick) THEN lthick=1.

;IF NOT KEYWORD_SET(xticksep) THEN xticksep=20

;the interactive part of the routine, pretty self-explanatory
IF NOT KEYWORD_SET(nocon) THEN $
;    READ,'Input level of significance, 0 -> 1 (i.e. 0.95 for 95%): ',sigg
sigg=0.99

IF NOT KEYWORD_SET(color) THEN color=1.5
LOADCT,color
;GAMMA_CT=0.4
;BLUE/WHITE=1

;
;========================= PLOT TIME SERIES ===============================
;
;again, dodge a window call that'd only confuse the PostScript writer.
IF KEYWORD_SET(print) THEN GOTO, SKIPPING

WINDOW,30,title='(30) Wavelet Plot',xsize=860,ysize=900
;WINDOW,30,title='(30) Wavelet Plot',xsize=640,ysize=710
;WINDOW,30,title='(30) Wavelet Plot',xsize=800,ysize=640

SKIPPING:

;set up the device co-ordinates for the top (time series) panel
pos1 = [0.12,0.80,0.7,0.95]

IF KEYWORD_SET(fig) THEN GOTO, FIG
;plot the filtered time series for comparison with the wavelet power transform
;which will be plotted below it once it's been calculated

;s=SIN(2. * !PI * (time+125) / 220.)

PLOT,time,series,yrange=[min(series),max(series)],$
       xtitle='Time [s]', ytitle='DN', charsize=2.3,$
       title=title,xthick=lthick,ythick=lthick,$
       position=pos1,/NOERASE,xstyle=1,ystyle=1,$
       xtickinterval=xticksep,xtickformat='(I9)'

;PLOT,time,series,yrange=[-10,12],$
;       xtitle='Time [s]', ytitle='DN', charsize=2,$
;       title=title,xthick=lthick,ythick=lthick,$
;       position=pos1,/NOERASE,xstyle=1,ystyle=1;,$
;       xtickinterval=xticksep,xtickformat='(I9)'

;OPLOT,time,0.09*s,lines=2,color=100

FIG:

;
;========================== WAVELET TRANSFORM ==============================
;
;This is the guts of the operation. The inputs are series (the copy of the time series
;made above), DELT (the time step between the data points) and SIGG (which specifies the
;percentage level you want WAVELET.PRO to calculate the significance contours for. The
;keyword PAD is set to pad the start and the end of the time series with values to
;minimise the edge effects inherent in a finite time series
;
;The outputs are as follows:
;   PERIOD is the list of the periods at which the wavelet transform was calculated
;   SCALE is the list of the corresponding wavelet scales
;   COI is the list of Cone-Of-Influence values so that it can be plotted over the
;   	wavelet transform
;   SIGGLVL is set to SIGG, the significance level
;   SIGNIF is the array of significance contours calculated for SIGG percent
;   FFT_THEOR is a theoretical FFT background spectrum for the time series.


IF NOT KEYWORD_SET(mother) THEN BEGIN
    mother='morlet'
    param = 6
END


    wave = WAVELET(series,delt,SIGNIF=signif,PERIOD=period,SCALE=scale,COI=coi,siglvl=sigg,FFT_THEOR=fft_theor,mother=mother,/PAD)



;get the number of wavelet scales actually used.
nscale = N_ELEMENTS(period)

print,'Min. period is',min(period)
print,'Max. period is',max(period)

IF NOT KEYWORD_SET(short) THEN short=min(period)
IF NOT KEYWORD_SET(long) THEN long=max(period)



;we're interested in the POWER transform, not the amplitude or phase,
;so we square the values in the transform (WAVE), and call the power
;transform power. Eoghan's name, not mine
power=abs(wave)^2   	;the power plot from the complex image
;power=power*elem/(moment(curve))(1)
power=power^0.5 	    ;decrease the power so you can actually see what's there

;
;============================ PLOT THE WAVELET POWER TRANSFORM =====================
;
;set the number of filled contour steps from the minimum of the power transform
;to the maximum. Sometimes need to reduce this to around 8 for older postscript
;printers, which may not be able to handle a large number of vertices in the plot
;description in the PostScript language.

	levels=25
	step=(max(power)-min(power))/levels
	userlevels=indgen(levels)*step + min(abs(wave)^2)

;set the device plotting co-ordinates for the wavelet power spectrum
IF NOT KEYWORD_SET(fig) THEN begin
    pos2 = [pos1(0),0.33,pos1(2),0.68]
    wavetitle='b) Morlet wavelet power transform'
ENDIF ELSE begin
    pos2 = [pos1(0),0.33,pos1(2),0.68]
    wavetitle=''
    XYOUTS,0.4,0.960,'a) Morlet wavelet power transform',/norm,align=0.5,charsize=2
ENDELSE

POLYFILL,[pos2(0),pos2(2),pos2(2),pos2(0)],$
    [pos2(1),pos2(1),pos2(3),pos2(3)],/norm,color=0

pickpeakperiod,power,maxper,period,perpower

IF NOT KEYWORD_SET(ylog) THEN begin

    CONTOUR,power,time,period,position=pos2, $
    XSTYLE=1,XTITLE='Time [s]',YTITLE='Period [s]',$
    TITLE=waveTitle,charsize=2.3,nLEVELS=25,/FILL,/NOERASE,$
    ystyle=1,YRANGE=[short,long],xthick=lthick,ythick=lthick,$
    ;xtickinterval=xticksep,xtickformat='(I9)',$
    background=0;,ytickinterval=0.002,/over,/normal
    horline,(elem-1)*delt/(3*sqrt(2)),line=2,color=255

ENDIF ELSE begin

    CONTOUR,power,time,period,position=pos2, $
    XSTYLE=1,XTITLE='Time [s]',YTITLE='Period [s]',$
    TITLE=waveTitle,charsize=2.3,nLEVELS=25,/FILL,/NOERASE,$
    ystyle=1,YRANGE=[short,long],/ytype,xthick=lthick,ythick=lthick,$
    ;xtickinterval=xticksep,xtickformat='(I9)',$
    background=0;,ytickinterval=0.002,/over,/normal
    horline,(elem-1)*delt/(3*sqrt(2)),line=2,color=255

ENDELSE


; mark what frequency the maximum power occurs at at each time
for i=0,N_ELEMENTS(maxper)-1 do plots,$
    [time(i),time(i)],[maxper(i),maxper(i)],$
       psym=-3,color=0

n=n_elements(series)

;power has been dropped to the 1/5th power, so we're saying that POWER is the
;TRUE power transform, so that we can sum over time to create the global
;wavelet spectrum.
power1 = (ABS(wave))^2
global_ws = TOTAL(power1,1)/n    ; global wavelet spectrum (GWS)

;
;================= CALCULATE THE GLOBAL SIGNIFICANCE FOR THE GLOBAL WAVELET SPECTRUM ==============
;
;we calculate this as a function of frequency, to be plotted over GLOBAL_WS on the right
;hand side as a dotted line. It's not the main thing to base your analysis on,
;but it can be helpful extra guide.
;
;this ensures that DOF is a VECTOR rather than a single number.
;this is necessary for WAVE_SIGNIF with SIGTEST=1 (what we
;need)
dof = n - scale     	;the -scale has the effect of correcting for padding at edges

;calcluate the global significance level:
    global_signif = WAVE_SIGNIF(series,delt,scale,1, $
    LAG1=0.0,DOF=dof,siglvl=sigg,MOTHER=mother,param=param)




;grab the sizes of PERIOD (output from WAVELET) and
;TIME (set at the beginning of the routine).
s_scale=n_elements(period)
s_time=n_elements(time)

;SIGNIF is actually just a 1-D vector of significance values for each
;associated period looked at by WAVELET. In order to divide the power
;transform by the significance levels, the following use of REBIN copies
;the vector S_TIME number of times in the Y-direction so that it's a 2-D
;array, and then TRANSPOSE flips it about the -45 degree axis so that
;the X values (in period space) become Y values and the Y values (in
;time) become the X values, for use below. power will be divided through
;by this array, which now has compatible dimensions
signif = REBIN(TRANSPOSE(signif),s_time,s_scale)

;tailor it to match the solid contours we've already plotted
;for the wavelet power plot.
;signif = signif*elem/(moment(curve))(1)
signif=signif*sig_fudge
signif=signif^0.5

;this bit is a mess. I don't know why it's even printed to the terminal
;but perhaps there used to be a good reason for it.
sigg2a=sigg*100.
sigg2a=string(sigg2a)
sigg2a=strtrim(sigg2a,2)

    sigg2a=strtrim(string(sigg*100.),2)
    percent=strcompress(strmid(sigg2a,0,2)+'%',/remove_all)

;print,'sigg2a is '+sigg2a
sigg3a=strmid(sigg2a,0,5)
;print,'sigg3a is '+sigg3a
sigg4a=strcompress(sigg3a+'%',/remove_all)
;print,'sigg4a is ',sigg4a

;
;====================== OVERPLOT THE SIGNIFICANCE CONTOURS ===================
;
;this first IF statement is to avoid overplotting contours and the COI with
;its associated cross-hatching to indicate the 'forbidden' area which the COI
;demarcates.
IF KEYWORD_SET(nocon) then GOTO,NO_CONTOURS

;this IF statement takes care of my own preferences for printing and
;displaying in an X-window respectively. I prefer white contours for
;terminal display, and like to be able to set the contour colours for
;PostScript device output.
;
;If the /PRINT option is used, then the PC variable can be set to anything
;from 0 (black) to 255 (white assuming you're using one of the monochromatic
;colour tables like 0,1,3,8 or 9). I choose to use colour table 1 (blue),
;which is set above, but that's easily changeable. The colour set by PC
;controls the colour of the significance contours, the COI boundary and the
;COI cross-hatching. This is easily changed, too, though.

CONTOUR,power/signif,time,period,$
/OVERPLOT,LEVEL=1,C_ANNOT=percent,color=pc,thick=3,C_CHARTHICK=lthick,$
xtitle='Time [s]',yrange=[short,long],ystyle=1

NO_CONTOURS:

IF KEYWORD_SET(cone) THEN BEGIN

    x = [MIN(time),time,MAX(time)]
    y = [MAX(period),coi,MAX(period)]
    PLOTS,time,coi,COLOR=pc,NOCLIP=0,THICK=lthick
    POLYFILL,x,y,ORIEN=+45,SPACING=0.5,COLOR=pc,NOCLIP=0,THICK=lthick
    POLYFILL,x,y,ORIEN=-45,SPACING=0.5,COLOR=pc,NOCLIP=0,THICK=lthick
	    
END


;print the value of the period at the top of the COI, which is
;therefore the highest 'reliable' period.
;print,'Highest credible period is: '+$
;    ARR2STR(max(coi(middle-20:middle+20)),/trim)+' seconds'

;
;============= PLOT THE GLOBAL WAVELET SPECTRUM AND GLOBAL SIGNIFICANCE LEVEL ==============
;
pos3 = [0.75,pos2(1),0.95,pos2(3)]

IF NOT KEYWORD_SET(ylog) THEN begin

    PLOT,global_ws,period,/NOERASE,POSITION=pos3, $
        XTITLE='Power [DN]',charsize=0.,$;,TITLE='c) Global'
        yrange=[short,long],ystyle=9,thick=lthick,$
        xrange=[MIN(global_ws),MAX(global_ws)],xstyle=11,$
        xtickinterval=5e-4,xthick=lthick,ythick=lthick;,xtickformat='(E9.2)'

ENDIF ELSE begin

    PLOT,global_ws,period,/NOERASE,POSITION=pos3, $
        XTITLE='',charsize=0.,$;,TITLE='c) Global'
        yrange=[short,long],ystyle=9,/ytype,thick=lthick,$
        xrange=[MIN(global_ws),MAX(global_ws)],xstyle=11,$
        xtickinterval=5e-4,xthick=lthick,ythick=lthick;,xtickformat='(E9.2)'
    XYOUTS,0.85,0.30,'Power [DN]',/norm,align=0.5,charsize=2.3
ENDELSE

IF NOT KEYWORD_SET(fig) THEN begin
    XYOUTS,0.85,0.730,'c) Fourier and ',/norm,align=0.5,charsize=2.3
    XYOUTS,0.85,0.705,'   Global Spectra',/norm,align=0.5,charsize=2.3
ENDIF ELSE begin
;    XYOUTS,0.855,0.980,'b) Fourier and ',/norm,align=0.5,charsize=1.2
;    XYOUTS,0.855,0.955,'   Global Spectra',/norm,align=0.5,charsize=1.2
    XYOUTS,0.85,0.960,'b) Fourier and Global Spectra',/norm,align=0.5,charsize=2.3
ENDELSE

IF KEYWORD_SET(cone) THEN BEGIN

;;;;;    OPLOT,global_signif,period,LINES=3   ;;;;;;;;; THIS IS THE FOURIER SIG LEVELS PLOTTED IN c)
;	plots,[MIN(global_ws),MAX(global_ws)],$
;	[max(coi(middle-20:middle+20)),$
;	max(coi(middle-20:middle+20))],line=1
	horline,(elem-1)*delt/(3*sqrt(2)),line=2
	
END

fser = FOURIER2(curve,delt,/norm)
pser = 1./(fser(0,*))
pser = pser
var_curve=(MOMENT(curve))(1)

IF NOT KEYWORD_SET(ylog) THEN begin

    AXIS,/xaxis,xrange=[min(fser(1,*)),max(fser(1,*))],$
    	xticks=2,xstyle=2,/save,xtickinterval=10,xthick=lthick
    AXIS,/yaxis,yrange=[short,long],/save,/yst,ythick=lthick

    oplot,fser(1,*),pser,psym=-1,color=175,thick=lthick

;;;;; DETERMINE THE PREDICTED POWER AT A GIVEN CONFIDENCE LEVEL
    level=signif_conf(curve,sigg) * elem / var_curve
    verline,level,lines=1

ENDIF ELSE begin

    AXIS,/xaxis,xrange=[min(fser(1,*)),max(fser(1,*))],$
    	xticks=2,xstyle=2,/save,xtickinterval=10,xthick=lthick
    AXIS,/yaxis,yrange=[short,long],/ylog,/save,/yst,ythick=lthick

    oplot,fser(1,*),pser,psym=-1,color=175,thick=lthick

    level=signif_conf(curve,sigg) * elem / var_curve
;;;;;;     verline,level,/ylog,lines=1    ;;;;; THIS IS THE GLOBAL WAVELET LINE IN c)

ENDELSE

;
;=================================== PRINT THE MAIN TITLE ==================================
;prev=!p.charsize
;!P.charsize=1.5
;XYOUTS,0.5,0.01,mainTitle,/nor,align=0.5
;!p.charsize=prev

IF KEYWORD_SET(fast) THEN GOTO,nearly

; find the confidence levels by the Banerjee randomisation test
randomwave,curve,conf,delt=delt

IF (N_ELEMENTS(sigg) eq 0) THEN sigg=0.95
linelevel = sigg*100.

IF NOT KEYWORD_SET(fig) THEN begin
    pos4 = [pos1(0),0.06,pos1(2),0.2]
    randomtitle='d) Confidence level for peak power'
ENDIF ELSE begin
    pos4 = [pos1(0),0.09,pos1(2),0.25]
    randomtitle=''
    XYOUTS,0.4,0.3,'c) Confidence level for peak power',/norm,align=0.5,charsize=1.3    
ENDELSE

;plot,time,perpower,yst=3,position=pos4,/noerase,/xst
plot,time,conf,yst=1,position=pos4,charsize=2,/noerase,/xst,$
    ytitle='Confidence %',title=randomtitle,xthick=lthick,ythick=lthick,$
    xtitle='Time [s]',yrange=[min(conf) < ((100*sigg)-2),(max(conf) > ((100*sigg)+2))];,$
    ;xtickinterval=xticksep,xtickformat='(I9)'

plots,[min(time),max(time)],[linelevel,linelevel],lines=3

nearly:

outscale=scale
IF (N_ELEMENTS(period) GT 0) THEN outper=period ELSE outper=0
IF (N_ELEMENTS(f_inver) GT 0) THEN outfreq=f_inver ELSE outfreq=0
IF (N_ELEMENTS(time) GT 0) THEN outtime=time ELSE outtime=0
;end procedure
;STOP
END
