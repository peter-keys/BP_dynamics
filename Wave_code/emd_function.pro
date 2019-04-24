;+
; NAME:
;     EMD_FUNCTION ; slight modifiaction of EMDECOMP.pro
;
; PURPOSE:
;     The purpose of EMD is to perform the Empirical Mode Decomposition,
;     method created by Huang, N. E. and others. EMD decomposes a signal
;     in intrinsic mode functions (imf), which can be used to build the
;     power-spectra (using hilbertspec.pro) like with Wavelets.
;
; AUTHOR:
;       JAUME TERRADAS CALAFELL
;       Departament de Fisica
; 	Universitat de les Illes Balears
;       E-mail: jaume.terradas@uib.es
;
; CATEGORY:
;     Signal Processing.
;
; CALLING SEQUENCE:
;
;     emdecomp,signal,sda,imf,DT=dt,NWIN1=nwin1,NWIN2=nwin2,PSOUT=psout,$
;              PSFILE=psfile
;
; INPUTS:
;     signal:    A 1D time series
;     sda:       Standard deviation to be achieved before accepting an IMF
;                (recommended value between 0.2 and 0.3; perhaps even smaller)
;
; KEYWORD PARAMETERS:
;     dt:        Time step between successive data values
;                (default dt=1). This parameter is only used for ;printing
;                purposes
;     nwin1:     Window used to plot the IMF (default nwin1=8)
;     nwin2:     Window used to plot the signal and residue (default nwin2=9)
;     psout:     Set this string variable to 'y' to produce ps output
;     psfile:    Filename used in to output plots in ps format
;
; OUTPUTS:
;     The variable imf contains the decomposed intrinsic mode functions
;     First index: number of component
;     Second index: the imf component
;
; SIDE EFFECTS:
;     Headache.....
;
; RESTRICTIONS:
;     The signal can be decomposed in a maximum number of nimf imf
;
; EXAMPLE:
;     EMD, signal, std_dev,DT=2.112,PSFILE='My_purdy_picture.ps'
;
; MODIFICATION HISTORY:
;      Written by:     Jaume Terradas, 29 DIC 2001.
;   Jaume Terradas, 20 JAN 2002. option to add characteristic waves at ends
;   Ramon Oliver, 2004, routines written in nice idl format
;   Jaume Terradas 22-6-2009, minor changes
;   Nabil Freij, 11/11/11, Changed the plotting routines
;   Peter Keys, Making it nicer/easier to understand and changing the plots
;   Peter Keys, Changed to function to output the plots as various arrays
;   Peter Keys, Chnaged the cap on number of IMFs (was 12 originally now 20)


FUNCTION EMD_FUNCTION,signal,sda,imf,DT=dt,PSFILE=psfile

;--- MAKE SURE SIGNAL IS DOUBLE PRECISION ---
ss=SIZE(signal)
IF (ss(0) LE 1) THEN signal = DOUBLE(signal)
IF (ss(0) GE 3) THEN BEGIN 
  STOP
ENDIF

;--- SET UP TIME STEP ---
IF (not keyword_set(dt)) THEN dt = 1.

;-- TO EXTEND AT THE ENDS OF EXTREMA ---
waveextensionm='y'

;--- TO CONTROL LARGE SWINGS AT THE ENDS ---
controlext='y'

nelem = N_ELEMENTS(signal)
time = DINDGEN(nelem)

;--- SETTING UP SOME IMF/EMD PARAMETERS ---
nimf = 20.					;Maximum number of IMFs that can be created
imf = DBLARR(nimf,nelem)
d = DBLARR(nelem-1)
h = DBLARR(nelem)
copy = DOUBLE(signal)
x = time
time1 = time
ncomp = 1

FOR t=0, nimf-1 DO BEGIN
  h = copy
  sd = 1.
  control = 0
  WHILE sd GT sda DO BEGIN
    FOR k = 0, nelem-2 DO BEGIN
    	d(k) = h(k+1)-h(k)
    ENDFOR

    maxmin=[0]
    
;--- DEAL WITH THE EXTREMA ---

    FOR i=0, nelem-3 DO BEGIN

      IF ((d(i) EQ 0.d0) AND (i NE 0)) THEN BEGIN
        IF (d(i-1)*d(i+1) LT 0.) THEN maxmin=[maxmin,i]
      ENDIF ELSE BEGIN
      IF (d(i)*d(i+1) lt 0.) THEN maxmin=[maxmin,i+1]
      ENDELSE
    ENDFOR
    IF (n_elements(maxmin) GT 1) THEN maxmin = maxmin(1:*)
    smaxmin = N_ELEMENTS(maxmin)
    IF smaxmin LE 2 THEN BEGIN
      control = -1
      GOTO,jump2
    ENDIF
    maxes=[0]
    mins=[0]
    IF h(maxmin(0)) GT h(maxmin(1)) THEN BEGIN
      FOR j=0,smaxmin-1,2 DO BEGIN
	      maxes = [maxes,maxmin(j)]
	      IF (j+1) LE (smaxmin-1) THEN mins = [mins,maxmin(j+1)]
      ENDFOR
    ENDIF ELSE BEGIN
      FOR j=0,smaxmin-1,2 DO BEGIN
	      mins = [mins,maxmin(j)]
	      IF (j+1) LE (smaxmin-1) THEN maxes=[maxes,maxmin(j+1)]
      ENDFOR

    ENDELSE

    maxes = maxes(1:*)
    mins = mins(1:*)
    nmax = N_ELEMENTS(maxes)
    nmin = N_ELEMENTS(mins)

;--- BEGIN EXTENDING AT THE ENDS OF THE EXTREMA ---

  IF waveextensionm EQ 'y' THEN BEGIN

    maxes1 = [0,maxes(0),maxes+maxes(0),-maxes(nmax-1)+2*(nelem-1)+maxes(0),2*((nelem-1)-maxes(nmax-1))+(nelem-1)+maxes(0)]
    hmodmax = [h(maxes(0)),h(maxes(0)),h(maxes),h(maxes(nmax-1)),h(maxes(nmax-1))]
    y2 = SPL_INIT(maxes1,hmodmax)
    maxenv = SPL_INTERP(maxes1,hmodmax,y2,time1+maxes(0))

    mins1 = [0,mins(0),mins+mins(0),-mins(nmin-1)+2*(nelem-1)+mins(0),2*((nelem-1)-mins(nmin-1))+(nelem-1)+mins(0)]
    hmodmin = [h(mins(0)),h(mins(0)),h(mins),h(mins(nmin-1)),h(mins(nmin-1))]
    y2 = SPL_INIT(mins1,hmodmin)
    minenv = SPL_INTERP(mins1,hmodmin,y2,time1+mins(0))

  ENDIF

    increm = 0.
    timeenvmax = time1
    timeenvmin = time1

    IF controlext EQ 'y' THEN BEGIN
      IF (h(0) GT maxenv(0)) AND (h(nelem-1) GT maxenv(nelem-1)) THEN BEGIN
        maxes1 = [0,maxes,nelem-1]
        hmodmax = [h(0)+increm,h(maxes),h(nelem-1)+increm]
        timeenvmax = time1
        y2 = SPL_INIT(maxes1,hmodmax)
        maxenv = SPL_INTERP(maxes1,hmodmax,y2,time1)
      ENDIF ELSE BEGIN
        IF h(0) GT maxenv(0) THEN BEGIN
          maxes1 = [0,maxes,-maxes(nmax-1)+2*(nelem-1),2*((nelem-1)-maxes(nmax-1))+(nelem-1)]
          hmodmax = [h(0)+increm,h(maxes),h(maxes(nmax-1)),h(maxes(nmax-1))]
          timeenvmax = time1
          y2 = SPL_INIT(maxes1,hmodmax)
          maxenv = SPL_INTERP(maxes1,hmodmax,y2,time1)
        ENDIF

        IF h(nelem-1) GT maxenv(nelem-1) THEN BEGIN
          maxes1 = [0,maxes(0),maxes+maxes(0),maxes(0)+nelem-1]
          hmodmax = [h(maxes(0)),h(maxes(0)),h(maxes),h(nelem-1)+increm]
          timeenvmax = time1+maxes(0)
          y2 = SPL_INIT(maxes1,hmodmax)
          maxenv = SPL_INTERP(maxes1,hmodmax,y2,time1+maxes(0))
        ENDIF
      ENDELSE

    IF (h(0) LT minenv(0)) AND (h(nelem-1) LT minenv(nelem-1)) THEN BEGIN
      mins1 = [0,mins,nelem-1]
      hmodmin = [h(0)-increm,h(mins),h(nelem-1)-increm]
      timeenvmin = time1
      y2 = SPL_INIT(mins1,hmodmin)
      minenv = SPL_INTERP(mins1,hmodmin,y2,time1)
    ENDIF ELSE BEGIN
      IF h(0) LT minenv(0) THEN BEGIN
        mins1 = [0,mins,-mins(nmin-1)+2*(nelem-1),2*((nelem-1)-mins(nmin-1))+(nelem-1)]
        hmodmin = [h(0)-increm,h(mins),h(mins(nmin-1)),h(mins(nmin-1))]
        timeenvmin = time1
        y2 = SPL_INIT(mins1,hmodmin)
        minenv = SPL_INTERP(mins1,hmodmin,y2,time1)
      ENDIF

      IF h(nelem-1) LT minenv(nelem-1) THEN BEGIN
        mins1 = [0,mins(0),mins+mins(0),mins(0)+nelem-1]
        hmodmin = [h(mins(0)),h(mins(0)),h(mins),h(nelem-1)-increm]
        timeenvmin = time1+mins(0)
        y2 = SPL_INIT(mins1,hmodmin)
        minenv = SPL_INTERP(mins1,hmodmin,y2,time1+mins(0))
      ENDIF
    ENDELSE

    ENDIF

    IF (ABS(minenv(0)) EQ !values.f_NaN) OR (ABS(maxenv(0)) EQ !values.f_NaN) THEN BEGIN
      STOP
    ENDIF

    m = (maxenv+minenv)/2.
    prevh = h
    h = h-m

    sd = total((prevh-h)^2/(prevh^2))
  ENDWHILE
  
jump2:
  ss = SIZE(maxmin)
  imf(t,*) = h(*)

  IF (ss(1) le 2) OR (control eq -1) THEN BEGIN
    ncomp = ncomp+1
    GOTO,jump1
  ENDIF
  ncomp = ncomp+1
  copy = copy-h
ENDFOR


jump1:
ncomp = ncomp-1

; --- REMOVE UNNEEDED SPACE IN THE ARRAY IMF ---
imf = imf(0:ncomp-1,*)

; --- RESULTS ARE NOW PLOTTED ---
WINDOW,0,xsize=900,ysize=900,title=' EMD - IMFs '

xtitle = 'Time (min)'
ytitle = strarr(ncomp)
ytitle(ncomp-1) ='r(t)'
FOR i = 0,ncomp-2 DO ytitle(i) = 'c!d'+STRCOMPRESS(strtrim(i+1),/REMOVE)+'!n Amplitude'

ncompv = ncomp

;SET_PLOT,'ps'
;DEVICE,/ENCAPSULATED,/PORT,/INCH,XSIZE=6,YSIZE=5,/COLOR,BITS=8,/helvetica,decomposed=0
;DEVICE,FILENAME=psfile
;DEVICE,/ENCAPSULATED,/PORT,/INCH,xsize=xsize,ysize=ysize,/helvetica,FILENAME=psfile,decomposed=0

;LOADCT,13,/SILENT		; RAINBOW COLOURED
LOADCT,0,/SILENT
;!P.FONT = 0  
;!P.CHARSIZE = 1
!p.background=255.
!p.color=0.
!p.multi = [0,2,(ncompv)/2]

FOR i=0,ncompv-1 DO BEGIN
  cgplot,time*dt/60,imf(i,*),charsize=1,xtitle=xtitle, ytitle=ytitle(i)
ENDFOR

; --- RECONSTRUCTION OF THE ORIGINAL SIGNAL ---
;sum = 0.d0
;FOR i = 0,ncomp-1 DO BEGIN
;  sum = sum+imf(i,*)
;ENDFOR

;cgplot,time*dt/60,signal,xstyle=2,ystyle=2,xtitle=xtitle,charsize=1
;cgplots,time*dt/60,imf(ncomp-1,*),color=200,thick=2

;DEVICE,/CLOSE

;SET_PLOT,'x'
!p.multi=0
LOADCT,0,/SILENT

RETURN, imf
END
