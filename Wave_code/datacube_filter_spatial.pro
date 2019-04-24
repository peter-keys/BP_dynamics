
PRO datacube_filter_spatial, datacube, cadence, arcsecpx, period_filters, spatial_sizes, save_directory=save_directory, save_memory=save_memory, super_save_memory=super_save_memory, make_plot=make_plot, no_komega=no_komega

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; CHECK FOR DEPENDENCIES FIRST:
;       
;
; datacube = INPUT DATACUBE, NORMALLY IN THE FORM OF [x, y, t] - BUT SHOULD BE A FILENAME IN THE CASE OF USING THE KEYWORD super_save_memory (BUT SAVE FILES MUST HAVE THE CUBE SAVED UNDER THE VARIABLE datacube)
; cadence = DELTA TIME BETWEEN SUCESSIVE FRAMES - GIVEN IN SECONDS
; arcsecpx = SPATIAL SAMPLING OF THE INPUT DATACUBE - GIVEN IN ARCSECONDS PER PIXEL
; period_filters = AN ARRAY OF PERIODS THAT YOU WISH TO FILTER - GIVEN IN SECONDS
; spatial_sizes = AN ARRAY OF SPATIAL SIZES OVER WHICH YOU WISH TO FILTER - GIVEN IN ARCSECONDS
; save_directory = THE LOCATION WHERE YOU WANT THE OUTPUT FILES TO BE SAVED - DEFAULT IS './'
; save_memory = OPTIONAL KEYWORD THAT HELPS TO SAVE MEMORY (RAM) - SEE BELOW
; super_save_memory = OPTIONAL KEYWORD THAT HELPS TO SAVE A LOT OF MEMORY (RAM) - DO NOT USE WITH save_memory!!!
;            WITH NO save_memory OR super_save_memory KEYWORDS THE CREATED ARRAYS ARE:
;                    datacube = [x, y, t] - ORIGINAL DATACUBE
;                    threedft = [x, y, t] - FOURIER DATACUBE
;                    filtered_FFT_cube = [x, y, t] - FOURIER FILTERED DATACUBE
;                    filtered_cube = [x, y, t] - FINAL FILTERED DATACUBE
;            WITH save_memory KEYWORD THE CREATED ARRAYS ARE:
;                    datacube = [x, y, t] - ORIGINAL DATACUBE (OVERWRITTEN BY FINAL FILTERED DATACUBE)
;                    threedft = [x, y, t] - FOURIER DATACUBE
;                    filtered_FFT_cube = [x, y, t] - FOURIER FILTERED DATACUBE
;            WITH super_save_memory KEYWORD THE CREATED ARRAYS ARE:
;                    datacube = [x, y, t] - ORIGINAL DATACUBE (OVERWRITTEN BY FINAL FILTERED DATACUBE)
;                    threedft = [x, y, t] - FOURIER DATACUBE (OVERWRITTEN BY FOURIER FILTERED DATACUBE)
; NOTE THAT USING THE super_save_memory KEYWORD WILL RESULT IN MUCH LONGER PROCESSING TIMES AS THE 3D FFT OF THE ORIGINAL DATACUBE MUST BE TAKEN FOR EACH PERIOD/SIZE LOOP SINCE WE CANNOT SAVE ITS CONTENTS		       		  
; make_plot = OPTIONAL KEYWORD THAT SAVES THE RESULTING K-OMEGA DIAGRAM AS AN AESTHETICALLY PLEASING EPS FILE IN THE CHOSEN save_directory
; no_komega = OPTIONAL KEYWORD THAT WILL NOT PRE-COMPUTE A K-OMEGA DIAGRAM, THUS SAVING SOME TIME AND MEMORY, BUT MISSING OUT ON 'COOL' FIGURES
;
; USAGE EXAMPLE: datacube_filter_spatial, datacube, 19., 0.18, [45,60,120,180,240,300,360,420,480,540,600], [0.4,0.6,1,3,5,10,15,20,30,40,50], save_directory='./', /save_memory, /make_plot
;                datacube_filter_spatial, 'Halpha_5s_cadence_140x1_1969x2_10y1_1899y2_integer.fits', 5., 0.09, [30,60,120,240,360,600], [0.2,0.5,1,5,10,20,50], save_directory='./', /super_save_memory, /no_komega
;                datacube_filter_spatial, 'Halpha_5s_cadence_1hr_140x1_1969x2_10y1_1899y2_integer.fits', 5., 0.09, [30,60,120,240,360,600], [0.2,0.5,1,5,10,20,50], save_directory='./', /super_save_memory, /make_plot
;                datacube_filter_spatial, datacube, 1.78, 0.138, [30,120,240,360,600,1200], [0.3,1.,5.,10.,30.,50.], save_directory='./', /save_memory, /make_plot
;
; THIS WILL RESULT IN FREQUENCY BINS OF 45:60, 60:120, 120:180 seconds, etc.
; THIS WILL RESULT IN SPATIAL BINS OF 1:3, 3:5, 5:10, 10:15, 15:20 arcsec, etc.
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


; IF A DATACUBE IS PASSED IN THEN WE DON'T NEED TO LOAD ANYTHING
IF ISA(datacube,/string) eq 0 THEN loading_method = 0

; IF A FILENAME IS PASSED IN AS datacube THEN THIS WILL NEED TO BE LOADED, EITHER THROUGH 'RESTORE' OR 'READFITS' BUT THIS DEPENDS ON THE FILETYPE
IF ISA(datacube,/string) eq 1 THEN BEGIN
    final_dot = STRPOS(datacube, '.', /reverse_search)
    file_ending = STRMID(datacube, final_dot)
    IF (file_ending eq '.sav') OR (file_ending eq '.save') THEN loading_method = 1
    IF (file_ending eq '.fit') OR (file_ending eq '.fits') THEN loading_method = 2
    filename = datacube
ENDIF    

; DEPENDING ON WHETHER A FITS/SAVE FILE IS PASSED INTO THE PROGRAMME, WE NEED TO LOAD THE APPROPRIATE DATA INTO THE datacube VARIABLE
IF loading_method eq 1 THEN restore,filename
IF loading_method eq 2 THEN datacube = READFITS(filename,/silent)

; CHOOSE THE SMALLEST/LARGEST STRUCTURE SIZE TO PLOT IN ARCSEC
smallest_size = arcsecpx * 2. ; ARCSEC
IF N_ELEMENTS(datacube[*,0,0]) le N_ELEMENTS(datacube[0,*,0]) THEN largest_size = arcsecpx * N_ELEMENTS(datacube[*,0,0]) ; ARCSEC
IF N_ELEMENTS(datacube[*,0,0]) gt N_ELEMENTS(datacube[0,*,0]) THEN largest_size = arcsecpx * N_ELEMENTS(datacube[0,*,0]) ; ARCSEC
; CHOOSE THE SMALLEST/LARGEST PERIODS TO PLOT
smallest_period = cadence * 2 ; SECONDS
largest_period = cadence * (N_ELEMENTS(datacube[0,0,*])/2.) ; SECONDS

; DEFINE WHERE ALL THE SAVED OUTPUTS WILL GO (k-omega DIAGRAM AND FILTERED TIME SERIES)
IF NOT KEYWORD_SET(save_directory) THEN save_directory = './'

range_max = MIN(period_filters)
new_range_max = (1./range_max) 

xsize_cube = N_ELEMENTS(datacube[*,0,0])
ysize_cube = N_ELEMENTS(datacube[0,*,0])
zsize_cube = N_ELEMENTS(datacube[0,0,*])
period_bins = N_ELEMENTS(period_filters) - 1
spatial_bins = N_ELEMENTS(spatial_sizes) - 1

print,''
print,''
print,''
IF loading_method eq 0 THEN print,'You have supplied an input datacube (i.e., not a file name!)'
IF loading_method eq 1 THEN print,'The input datacube is an IDL SAVE file with the file name: ' + filename
IF loading_method eq 2 THEN print,'The input datacube is a FITS file with the file name: ' + filename
print,'The input datacube is of size: ['+ARR2STR(xsize_cube,/trim)+', '+ARR2STR(ysize_cube,/trim)+', '+ARR2STR(zsize_cube,/trim)+']'
print,''
print,'The 2-pixel size is '+ARR2STR(smallest_size,/trim)+' arcsec, which corresponds to a spatial frequency of '+ARR2STR(((2.*!pi)/FLOAT(smallest_size)),/trim)+' arcsec^-1'
print,'The duration of the time series is '+ARR2STR(cadence*zsize_cube,/trim)+' seconds with a Nyquist period of '+ARR2STR(smallest_period,/trim)+' seconds (or '+ARR2STR((1./FLOAT(smallest_period))*1000.,/trim)+' mHz)'
print,''

IF NOT KEYWORD_SET(no_komega) THEN BEGIN
    sp_out = DBLARR(xsize_cube/2.,zsize_cube/2.)
    print,''
    print,'Constructing a k-omega diagram of the input datacube..........'
    ; MAKE THE k-omega DIAGRAM USING THE PROVEN METHOD OF ROB RUTTEN
    kopower = plotkopower_funct(datacube, sp_out, arcsecpx, cadence, apod=0.,  kmax=1., fmax=1.)

    ; MAKE THE FULL RANGE OF WAVENUMBERS
    ; CALCULATED FOR THE PLOTTING OF THE K-OMEGA DIAGRAM
    full_xrange = 2.*!pi/(arcsecpx*2.)*FINDGEN((xsize_cube/2.))/((xsize_cube/2.)-1)
    wavenumber_elements = WHERE(full_xrange ge ((2.*!pi)/largest_size) AND full_xrange le ((2.*!pi)/smallest_size))
    wavemin = MIN(wavenumber_elements)
    wavemax = MAX(wavenumber_elements)
    new_xrange = REFORM(full_xrange[wavenumber_elements])

    ; MAKE THE FULL RANGE OF FREQUENCIES
    ; CALCULATED FOR THE PLOTTING OF THE K-OMEGA DIAGRAM
    full_yrange = (1./(cadence*2.)) * (FINDGEN(zsize_cube/2.)/((zsize_cube/2.)-1)) * 1e3
    frequency_elements = WHERE(full_yrange ge ((1./largest_period)*1000.) AND full_yrange le ((1./smallest_period)*1000.))
    freqmin = MIN(frequency_elements)
    freqmax = MAX(frequency_elements)
    new_yrange = REFORM(full_yrange[frequency_elements])

    smoothkopower = CONVOL(REBIN(kopower[wavemin:wavemax,freqmin:freqmax], (wavemax-wavemin+1)*5, (freqmax-freqmin+1)*5), FLTARR(5,5)+1/(5.*5.), /edge_truncate)
    kopower_map = MAKE_MAP(ALOG10(smoothkopower)-MIN(ALOG10(smoothkopower)), dx=(MAX(new_xrange)-MIN(new_xrange))/((wavemax-wavemin+1.)*5.), dy=(MAX(new_yrange)-MIN(new_yrange))/((freqmax-freqmin+1.)*5.), xc=((MAX(new_xrange)-MIN(new_xrange))/2.)+MIN(new_xrange), yc=((MAX(new_yrange)-MIN(new_yrange))/2.)+MIN(new_yrange), time='', units='arcsecs')
    kopower_max = MAX(kopower_map.data)
    kopower_min = MIN(kopower_map.data)
ENDIF

print,''
print,''
print,'Making a 3D Fourier transform of the input datacube..........'
print,''
threedft = FFT(datacube,-1)

; WORK OUT THINGS FOR DISPLAYING PURPOSES
; TEMPORAL FREQUENCY
; CALCULATED FOR THE TIME FILTERING OF WAVES
mid_point_time = (N_ELEMENTS(threedft[0,0,*])/2) + 1 
freq = INDGEN(N_ELEMENTS(threedft[0,0,*])) 
IF (N_ELEMENTS(threedft[0,0,*]) MOD 2) ne 0 THEN freq[mid_point_time] = mid_point_time - N_ELEMENTS(threedft[0,0,*]) + FINDGEN(mid_point_time-1)
IF (N_ELEMENTS(threedft[0,0,*]) MOD 2) eq 0 THEN freq[mid_point_time] = mid_point_time - N_ELEMENTS(threedft[0,0,*]) + FINDGEN(mid_point_time-2)
time_freq = freq / (N_ELEMENTS(threedft[0,0,*])*cadence) 
; WORK OUT THINGS FOR DISPLAYING PURPOSES
; SPATIAL_x FREQUENCY
; CALCULATED FOR THE SPATIAL-X FILTERING OF WAVES
mid_point = (N_ELEMENTS(threedft[*,0,0])/2) + 1 
freq = INDGEN(N_ELEMENTS(threedft[*,0,0])) 
IF (N_ELEMENTS(threedft[*,0,0]) MOD 2) ne 0 THEN freq[mid_point] = mid_point - N_ELEMENTS(threedft[*,0,0]) + FINDGEN(mid_point-1)
IF (N_ELEMENTS(threedft[*,0,0]) MOD 2) eq 0 THEN freq[mid_point] = mid_point - N_ELEMENTS(threedft[*,0,0]) + FINDGEN(mid_point-2)
spatial_x_freq = freq / (N_ELEMENTS(threedft[*,0,0])*(arcsecpx/(2.*!pi))) 
; WORK OUT THINGS FOR DISPLAYING PURPOSES
; SPATIAL_y FREQUENCY
; CALCULATED FOR THE SPATIAL-Y FILTERING OF WAVES
mid_point = (N_ELEMENTS(threedft[0,*,0])/2) + 1 
freq = INDGEN(N_ELEMENTS(threedft[0,*,0])) 
IF (N_ELEMENTS(threedft[0,*,0]) MOD 2) ne 0 THEN freq[mid_point] = mid_point - N_ELEMENTS(threedft[0,*,0]) + FINDGEN(mid_point-1)
IF (N_ELEMENTS(threedft[0,*,0]) MOD 2) eq 0 THEN freq[mid_point] = mid_point - N_ELEMENTS(threedft[0,*,0]) + FINDGEN(mid_point-2)
spatial_y_freq = freq / (N_ELEMENTS(threedft[0,*,0])*(arcsecpx/(2.*!pi))) 

; MAKE A SPATIALLY AVERAGED POWER SPECTRUM
threedft_spatial_sum = REFORM(TOTAL(TOTAL(threedft,2),1)) / (FLOAT(N_ELEMENTS(threedft[*,0,0])) * FLOAT(N_ELEMENTS(threedft[0,*,0])))
; MAKE A TEMPORALLY AVERAGED POWER SPECTRA AND CONVERT TO A MAP
threedft_time_sum = REFORM(TOTAL(threedft,3)) / (FLOAT(N_ELEMENTS(threedft[0,0,*])))
;IF KEYWORD_SET(no_komega) THEN 
threedft_max = MAX(REAL_PART(ALOG10(SHIFT(threedft_time_sum,(N_ELEMENTS(threedft[*,0,0])/2),(N_ELEMENTS(threedft[0,*,0])/2)))-MIN(ALOG10(SHIFT(threedft_time_sum,(N_ELEMENTS(threedft[*,0,0])/2),(N_ELEMENTS(threedft[0,*,0])/2))))))
;IF KEYWORD_SET(no_komega) THEN 
threedft_min = MIN(REAL_PART(ALOG10(SHIFT(threedft_time_sum,(N_ELEMENTS(threedft[*,0,0])/2),(N_ELEMENTS(threedft[0,*,0])/2)))-MIN(ALOG10(SHIFT(threedft_time_sum,(N_ELEMENTS(threedft[*,0,0])/2),(N_ELEMENTS(threedft[0,*,0])/2))))))
IF (xsize_cube MOD 2 eq 0) AND (ysize_cube MOD 2 eq 0) THEN threedft_time_sum_map = MAKE_MAP(ALOG10(SHIFT(threedft_time_sum,((N_ELEMENTS(threedft[*,0,0])/2)-1),((N_ELEMENTS(threedft[0,*,0])/2)-1)))-MIN(ALOG10(SHIFT(threedft_time_sum,((N_ELEMENTS(threedft[*,0,0])/2)-1),((N_ELEMENTS(threedft[0,*,0])/2)-1)))), dx=(MAX(spatial_x_freq)-MIN(spatial_x_freq))/(FLOAT(N_ELEMENTS(spatial_x_freq))), dy=(MAX(spatial_y_freq)-MIN(spatial_y_freq))/(FLOAT(N_ELEMENTS(spatial_y_freq))), xc=((MAX(spatial_x_freq)-MIN(spatial_x_freq))/2.)+MIN(spatial_x_freq), yc=((MAX(spatial_y_freq)-MIN(spatial_y_freq))/2.)+MIN(spatial_y_freq), time='', units='arcsecs')
IF (xsize_cube MOD 2 eq 0) AND (ysize_cube MOD 2 ne 0) THEN threedft_time_sum_map = MAKE_MAP(ALOG10(SHIFT(threedft_time_sum,((N_ELEMENTS(threedft[*,0,0])/2)-1),((N_ELEMENTS(threedft[0,*,0])/2))))-MIN(ALOG10(SHIFT(threedft_time_sum,((N_ELEMENTS(threedft[*,0,0])/2)-1),((N_ELEMENTS(threedft[0,*,0])/2))))), dx=(MAX(spatial_x_freq)-MIN(spatial_x_freq))/(FLOAT(N_ELEMENTS(spatial_x_freq))), dy=(MAX(spatial_y_freq)-MIN(spatial_y_freq))/(FLOAT(N_ELEMENTS(spatial_y_freq))), xc=((MAX(spatial_x_freq)-MIN(spatial_x_freq))/2.)+MIN(spatial_x_freq), yc=((MAX(spatial_y_freq)-MIN(spatial_y_freq))/2.)+MIN(spatial_y_freq), time='', units='arcsecs')
IF (xsize_cube MOD 2 ne 0) AND (ysize_cube MOD 2 eq 0) THEN threedft_time_sum_map = MAKE_MAP(ALOG10(SHIFT(threedft_time_sum,((N_ELEMENTS(threedft[*,0,0])/2)),((N_ELEMENTS(threedft[0,*,0])/2)-1)))-MIN(ALOG10(SHIFT(threedft_time_sum,((N_ELEMENTS(threedft[*,0,0])/2)),((N_ELEMENTS(threedft[0,*,0])/2)-1)))), dx=(MAX(spatial_x_freq)-MIN(spatial_x_freq))/(FLOAT(N_ELEMENTS(spatial_x_freq))), dy=(MAX(spatial_y_freq)-MIN(spatial_y_freq))/(FLOAT(N_ELEMENTS(spatial_y_freq))), xc=((MAX(spatial_x_freq)-MIN(spatial_x_freq))/2.)+MIN(spatial_x_freq), yc=((MAX(spatial_y_freq)-MIN(spatial_y_freq))/2.)+MIN(spatial_y_freq), time='', units='arcsecs')
IF (xsize_cube MOD 2 ne 0) AND (ysize_cube MOD 2 ne 0) THEN threedft_time_sum_map = MAKE_MAP(ALOG10(SHIFT(threedft_time_sum,((N_ELEMENTS(threedft[*,0,0])/2)),((N_ELEMENTS(threedft[0,*,0])/2))))-MIN(ALOG10(SHIFT(threedft_time_sum,((N_ELEMENTS(threedft[*,0,0])/2)),((N_ELEMENTS(threedft[0,*,0])/2))))), dx=(MAX(spatial_x_freq)-MIN(spatial_x_freq))/(FLOAT(N_ELEMENTS(spatial_x_freq))), dy=(MAX(spatial_y_freq)-MIN(spatial_y_freq))/(FLOAT(N_ELEMENTS(spatial_y_freq))), xc=((MAX(spatial_x_freq)-MIN(spatial_x_freq))/2.)+MIN(spatial_x_freq), yc=((MAX(spatial_y_freq)-MIN(spatial_y_freq))/2.)+MIN(spatial_y_freq), time='', units='arcsecs')


print,''
print,'Displaying spectra..........'
print,''
IF (NOT KEYWORD_SET(no_komega)) AND (KEYWORD_SET(make_plot)) THEN BEGIN
    set_plot, 'ps'
    xsize=24 ; xsize in cm
    ysize=20
    device, filename = save_directory + 'Filtered_datacube_k_omega_diagram.eps', /encapsulated, xsize=xsize, ysize=ysize, /tt_font, set_font='Times', font_size=4, /color, bits_per_pixel=8
    LOADCT,0,/silent
    !p.background = 255.
    !p.color = 0.
    mult,1,1
    LOADCT,39,/silent
    IF (MAX(new_yrange) le 50) THEN plot_map, kopower_map, charsize=3, xticklen=-.025, yticklen=-.025, xtickinterval=5, ytickinterval=5, dmin=kopower_min, dmax=(kopower_max-1), xtitle='Wavenumber (arcsec!U-1!N)', ytitle='Frequency (mHz)', /notitle, xst=9, yst=9, position=[0.1,0.1,0.9,0.83]
    IF (MAX(new_yrange) gt 50) AND (MAX(new_yrange) le 100) THEN plot_map, kopower_map, charsize=3, xticklen=-.025, yticklen=-.025, xtickinterval=5, ytickinterval=10, dmin=kopower_min, dmax=(kopower_max-1), xtitle='Wavenumber (arcsec!U-1!N)', ytitle='Frequency (mHz)', /notitle, xst=9, yst=9, position=[0.1,0.1,0.9,0.83]
    IF (MAX(new_yrange) gt 100) AND (MAX(new_yrange) le 200) THEN plot_map, kopower_map, charsize=3, xticklen=-.025, yticklen=-.025, xtickinterval=5, ytickinterval=20, dmin=kopower_min, dmax=(kopower_max-1), xtitle='Wavenumber (arcsec!U-1!N)', ytitle='Frequency (mHz)', /notitle, xst=9, yst=9, position=[0.1,0.1,0.9,0.83]
    IF (MAX(new_yrange) gt 200) AND (MAX(new_yrange) le 500) THEN plot_map, kopower_map, charsize=3, xticklen=-.025, yticklen=-.025, xtickinterval=5, ytickinterval=50, dmin=kopower_min, dmax=(kopower_max-1), xtitle='Wavenumber (arcsec!U-1!N)', ytitle='Frequency (mHz)', /notitle, xst=9, yst=9, position=[0.1,0.1,0.9,0.83]
    IF (MAX(new_yrange) gt 500) AND (MAX(new_yrange) le 1000) THEN plot_map, kopower_map, charsize=3, xticklen=-.025, yticklen=-.025, xtickinterval=5, ytickinterval=100, dmin=kopower_min, dmax=(kopower_max-1), xtitle='Wavenumber (arcsec!U-1!N)', ytitle='Frequency (mHz)', /notitle, xst=9, yst=9, position=[0.1,0.1,0.9,0.83]
    IF (MAX(new_yrange) gt 1000) THEN plot_map, kopower_map, charsize=3, xticklen=-.025, yticklen=-.025, xtickinterval=5, ytickinterval=250, dmin=kopower_min, dmax=(kopower_max-1), xtitle='Wavenumber (arcsec!U-1!N)', ytitle='Frequency (mHz)', /notitle, xst=9, yst=9, position=[0.1,0.1,0.9,0.83]
    AXIS, XSTYLE=1, XAXIS=1, XTICKFORMAT='convert_wavenumber_to_distance', CHARSIZE=3, XTITLE='Spatial size (arcsec)', xticklen=-.025
    AXIS, YSTYLE=1, YAXIS=1, YTICKFORMAT='convert_mHz_to_period', CHARSIZE=3, YTITLE='Period (s)', yticklen=-.025
    tickmarknames = STRARR(FIX(kopower_max))
    FOR tick = 0, (FIX(kopower_max)-1) DO tickmarknames[tick] = tick
    colorbar, bottom=0, ncolors=255, charsize=1.5, color=0, divisions=FIX(kopower_max)-1, minrange=kopower_min, maxrange=(kopower_max-1), position=[0.1, 0.92, 0.9, 0.95], /horizontal, /top, ticknames=tickmarknames, title='Oscillation Power (log)', font=-1
    device,/close
    set_plot,'x'
ENDIF
IF (NOT KEYWORD_SET(no_komega)) THEN BEGIN
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    ;k-omega diagram
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    window,1,xsize=1200,ysize=800,title='k-omega diagram'
    loadct,0,/silent
    !p.color=0.
    !p.background = 255.
    mult,1,1
    loadct,39,/silent
    IF (MAX(new_yrange) le 50) THEN plot_map, kopower_map, charsize=3, xticklen=-.025, yticklen=-.025, xtickinterval=5, ytickinterval=5, dmin=kopower_min, dmax=(kopower_max-1), xtitle='Wavenumber (arcsec!U-1!N)', ytitle='Frequency (mHz)', /notitle, xst=9, yst=9, position=[0.1,0.1,0.9,0.83]
    IF (MAX(new_yrange) gt 50) AND (MAX(new_yrange) le 100) THEN plot_map, kopower_map, charsize=3, xticklen=-.025, yticklen=-.025, xtickinterval=5, ytickinterval=10, dmin=kopower_min, dmax=(kopower_max-1), xtitle='Wavenumber (arcsec!U-1!N)', ytitle='Frequency (mHz)', /notitle, xst=9, yst=9, position=[0.1,0.1,0.9,0.83]
    IF (MAX(new_yrange) gt 100) AND (MAX(new_yrange) le 200) THEN plot_map, kopower_map, charsize=3, xticklen=-.025, yticklen=-.025, xtickinterval=5, ytickinterval=20, dmin=kopower_min, dmax=(kopower_max-1), xtitle='Wavenumber (arcsec!U-1!N)', ytitle='Frequency (mHz)', /notitle, xst=9, yst=9, position=[0.1,0.1,0.9,0.83]
    IF (MAX(new_yrange) gt 200) AND (MAX(new_yrange) le 500) THEN plot_map, kopower_map, charsize=3, xticklen=-.025, yticklen=-.025, xtickinterval=5, ytickinterval=50, dmin=kopower_min, dmax=(kopower_max-1), xtitle='Wavenumber (arcsec!U-1!N)', ytitle='Frequency (mHz)', /notitle, xst=9, yst=9, position=[0.1,0.1,0.9,0.83]
    IF (MAX(new_yrange) gt 500) AND (MAX(new_yrange) le 1000) THEN plot_map, kopower_map, charsize=3, xticklen=-.025, yticklen=-.025, xtickinterval=5, ytickinterval=100, dmin=kopower_min, dmax=(kopower_max-1), xtitle='Wavenumber (arcsec!U-1!N)', ytitle='Frequency (mHz)', /notitle, xst=9, yst=9, position=[0.1,0.1,0.9,0.83]
    IF (MAX(new_yrange) gt 1000) THEN plot_map, kopower_map, charsize=3, xticklen=-.025, yticklen=-.025, xtickinterval=5, ytickinterval=250, dmin=kopower_min, dmax=(kopower_max-1), xtitle='Wavenumber (arcsec!U-1!N)', ytitle='Frequency (mHz)', /notitle, xst=9, yst=9, position=[0.1,0.1,0.9,0.83]
    AXIS, XSTYLE=1, XAXIS=1, XTICKFORMAT='convert_wavenumber_to_distance', CHARSIZE=3, XTITLE='Spatial size (arcsec)', xticklen=-.025
    AXIS, YSTYLE=1, YAXIS=1, YTICKFORMAT='convert_mHz_to_period', CHARSIZE=3, YTITLE='Period (s)', yticklen=-.025
    tickmarknames = STRARR(FIX(kopower_max))
    FOR tick = 0, (FIX(kopower_max)-1) DO tickmarknames[tick] = tick
    ;colorbar, bottom=0, ncolors=255, charsize=2, color=0, divisions=5, minrange=0, maxrange=5, position=[0.95, 0.1, 0.97, 0.85], /vertical, /right, ticknames=tickmarknames, title='Oscillation Power (log)', font=-1
    colorbar, bottom=0, ncolors=255, charsize=1.5, color=0, divisions=FIX(kopower_max)-1, minrange=kopower_min, maxrange=(kopower_max-1), position=[0.1, 0.92, 0.9, 0.95], /horizontal, /top, ticknames=tickmarknames, title='Oscillation Power (log)', font=-1
ENDIF
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;Fourier spectra
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
window,0,xsize=1200,ysize=800,title='Fourier filtering spectrum'
loadct,0,/silent
!p.color=0.
!p.background = 255.
mult,2,1
loadct,5,/silent
plot_map, threedft_time_sum_map, charsize=2, xticklen=-.025, yticklen=-.025, xtickinterval=10, ytickinterval=10, dmin=(threedft_min+2), dmax=(threedft_max-2), xtitle='Wavenumber x (arcsec!U-1!N)', ytitle='Wavenumber y (arcsec!U-1!N)', title='Temporally averaged Fourier spectrum'
PLOTS, [MIN(spatial_x_freq),MAX(spatial_x_freq)], [0,0], line=2, thick=1
PLOTS, [0,0], [MIN(spatial_y_freq),MAX(spatial_y_freq)], line=2, thick=1
plot,SHIFT(time_freq, -mid_point_time), SHIFT(ABS(threedft_spatial_sum), -mid_point_time),/ylog,thick=3,$
     xtitle='Frequency (Hz)',ytitle='Power',title='Spatially averaged Fourier power',charsize=2,xst=1,yst=1,xr=[(-new_range_max),new_range_max]
;tvim,ALOG10(SHIFT(threedft_time_sum,(N_ELEMENTS(threedft[*,0,0])/2),(N_ELEMENTS(threedft[0,*,0])/2))),/sc,range=[-4,-2],xr=[MIN(spatial_x_freq),MAX(spatial_x_freq)],yr=[MIN(spatial_y_freq),MAX(spatial_y_freq)],pcharsize=2,title='Temporally averaged Fourier power',ytitle='Wavenumber (y)',xtitle='Wavenumber (x)'
;window,0,xsize=1800,ysize=1000,title='Fourier filtering spectrum'
;mult,2,1
;loadct,5,/silent
;plot_map, threedft_time_sum_map, charsize=2, xticklen=-.025, yticklen=-.025, xtickinterval=10, ytickinterval=10, dmin=(kopower_min+2), dmax=(kopower_max-2), xtitle='Wavenumber x (arcsec!U-1!N)', ytitle='Wavenumber y (arcsec!U-1!N)', title='Temporally averaged Fourier spectrum'
;plot,SHIFT(time_freq, -mid_point_time), SHIFT(ABS(threedft_spatial_sum), -mid_point_time),/ylog,thick=3,$
;     xtitle='Frequency (Hz)',ytitle='Power',title='Spatially averaged Fourier power',charsize=2,xst=1,yst=1,xr=[(-new_range_max),new_range_max]
;tvim,ALOG10(SHIFT(threedft_time_sum,(N_ELEMENTS(threedft[*,0,0])/2),(N_ELEMENTS(threedft[0,*,0])/2))),/sc,range=[-4,-2],xr=[MIN(spatial_x_freq),MAX(spatial_x_freq)],yr=[MIN(spatial_y_freq),MAX(spatial_y_freq)],pcharsize=2,title='Temporally averaged Fourier power',ytitle='Wavenumber (y)',xtitle='Wavenumber (x)'
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

WAIT, 0.5

countup = 0.
t0 = SYSTIME(1) 
t1 = t0

; SET UP AN INTEGER-ROUNDED PYTHAGORAS DISTANCE ARRAY THAT IS DIRECTLY MATCHED TO THE Kx and Ky FREQUENCY DOMAINS 
distances = ko_dist(xsize_cube,ysize_cube)      ;INTEGER-ROUNDED PYTHAGORAS DISTANCE ARRAY
maxdist = (MIN([xsize_cube,ysize_cube])/2)-1    ;LARGEST QUARTER CIRCLE

; CREATE A BINARY-TYPE MASKING MAP THAT WILL ONLY INCLUDE SPATIAL FREQUENCIES OF INTEREST
filter_map = FLTARR(xsize_cube,ysize_cube)

IF (NOT KEYWORD_SET(save_memory)) AND (NOT KEYWORD_SET(super_save_memory)) THEN filtered_cube = datacube
IF NOT KEYWORD_SET(super_save_memory) THEN filtered_fft_cube = threedft

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; START TEMPORAL AND SPATIAL FILTERING
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

FOR z = 0,(period_bins-1) DO BEGIN 
    WSET,0
    mult,2,1
    min_freq = 1./(period_filters[z+1])
    max_freq = 1./(period_filters[z]) 
    ; SET FREQUENCY FILTERING RANGES REMEMBERING TO INCLUDE NEGATIVE FREQUENCIES!!!
    filter_range1 = WHERE(time_freq gt (-min_freq) AND time_freq lt (min_freq)) 
    filter_range2 = WHERE(time_freq lt (-max_freq)) 
    filter_range3 = WHERE(time_freq gt (max_freq)) 
    temp_sum = threedft_spatial_sum
    temp_sum[filter_range1] = 0. 
    temp_sum[filter_range2] = 0. 
    temp_sum[filter_range3] = 0. 
    ;;;;;;;;;;;;;;;;;;;;;;;;;;
    ; START SPATIAL FILTERING
    ;;;;;;;;;;;;;;;;;;;;;;;;;;
    FOR i=0, (spatial_bins-1) DO BEGIN 
        IF KEYWORD_SET(super_save_memory) THEN BEGIN
	   datacube = 1. ; THIS JUST CLEARS THE MEMORY FROM THE ORIGINAL datacube ARRAY, WHICH IF NOT DONE CAN CAUSE MEMORY CRASHES!
	   IF loading_method eq 1 THEN restore,filename
           IF loading_method eq 2 THEN datacube = READFITS(filename,/silent)
	ENDIF
	IF KEYWORD_SET(super_save_memory) THEN threedft = 1. ; THIS JUST CLEARS THE MEMORY FROM THE ORIGINAL threedft ARRAY, WHICH IF NOT DONE CAN CAUSE MEMORY CRASHES!
	IF KEYWORD_SET(super_save_memory) THEN threedft = FFT(datacube,-1)
        ; TAKE THE USER-DEFINED SPATIAL SIZES AND CONVERT TO SPATIAL WAVENUMBERS
	min_spatial_size = spatial_sizes[i]
	max_spatial_size = spatial_sizes[i+1]
	min_spatial_size_wavenumber = (2.*!pi) / FLOAT(min_spatial_size)
	max_spatial_size_wavenumber = (2.*!pi) / FLOAT(max_spatial_size)
	; values_spatial_frequency ARE THE LOCATIONS WHERE THE PRE-COMPUTED SPATIAL FREQUENCIES MATCH THE DESIRED SPATIAL SIZES
	values_spatial_frequency = WHERE(((spatial_y_freq ge max_spatial_size_wavenumber) AND (spatial_y_freq le min_spatial_size_wavenumber)) OR ((spatial_x_freq ge max_spatial_size_wavenumber) AND (spatial_x_freq le min_spatial_size_wavenumber)))
	distance_min = distances[MIN(values_spatial_frequency),0]
	distance_max = distances[MAX(values_spatial_frequency),0]
	; FIND THE LOCATIONS WHERE THE distances ARRAY MATCHES THE SPATIAL-FREQUENCY PIXEL VALUES (NOW CONVERTED TO DISTANCES WITH MAX/MIN EQUAL TO distance_max/distance_min)
	; NOTE THAT IF THE x/y INPUT ARRAY SIZES ARE NOT IDENTICAL THIS WILL RESULT IN SLIGHTLY OVER-ESTIMATED kx/ky FREQUENCIES WITH: 
	;                           min_spatial_size_wavenumber > MIN(spatial_y_freq[values_spatial_frequency])   AND/OR
	;                           max_spatial_size_wavenumber<MAX(spatial_y_freq[values_spatial_frequency])     AND/OR
	;                           min_spatial_size_wavenumber > MIN(spatial_x_freq[values_spatial_frequency])   AND/OR
	;                           max_spatial_size_wavenumber<MAX(spatial_x_freq[values_spatial_frequency])
	values = WHERE((distances ge distance_min) AND (distances le distance_max))
        filter_map[*] = 0.
        filter_map[values] = 1.
        used_mask = filter_map
	used_mask[*] = 0.
	used_mask[values] = threedft_time_sum_map.data[values]
	; PLOT THE TEMPORALLY-AVERAGED SPATIAL FFT AND THE SPATIALLY-AVERAGED FREQUENCY FFT
	LOADCT,5,/silent
	plot_map, threedft_time_sum_map, charsize=2, xticklen=-.025, yticklen=-.025, xtickinterval=10, ytickinterval=10, dmin=(threedft_min+2), dmax=(threedft_max-2), xtitle='Wavenumber x (arcsec!U-1!N)', ytitle='Wavenumber y (arcsec!U-1!N)', title='Temporally averaged Fourier spectrum'
	PLOTS, [MIN(spatial_x_freq),MAX(spatial_x_freq)], [0,0], line=2, thick=1
	PLOTS, [0,0], [MIN(spatial_y_freq),MAX(spatial_y_freq)], line=2, thick=1
	filter_map_map = MAKE_MAP(SHIFT(filter_map, xsize_cube/2, ysize_cube/2), dx=(MAX(spatial_x_freq)-MIN(spatial_x_freq))/(FLOAT(N_ELEMENTS(spatial_x_freq))), dy=(MAX(spatial_y_freq)-MIN(spatial_y_freq))/(FLOAT(N_ELEMENTS(spatial_y_freq))), xc=((MAX(spatial_x_freq)-MIN(spatial_x_freq))/2.)+MIN(spatial_x_freq), yc=((MAX(spatial_y_freq)-MIN(spatial_y_freq))/2.)+MIN(spatial_y_freq), time='', units='arcsecs')
        used_mask_map = MAKE_MAP(SHIFT(used_mask, xsize_cube/2, ysize_cube/2), dx=(MAX(spatial_x_freq)-MIN(spatial_x_freq))/(FLOAT(N_ELEMENTS(spatial_x_freq))), dy=(MAX(spatial_y_freq)-MIN(spatial_y_freq))/(FLOAT(N_ELEMENTS(spatial_y_freq))), xc=((MAX(spatial_x_freq)-MIN(spatial_x_freq))/2.)+MIN(spatial_x_freq), yc=((MAX(spatial_y_freq)-MIN(spatial_y_freq))/2.)+MIN(spatial_y_freq), time='', units='arcsecs')
        LOADCT,39,/silent
	;plot_map, used_mask_map, /blend, /over, alpha=0.3, dmin=2, dmax=5
	plot_map, filter_map_map, /contour, /overlay, levels=[0.5], c_color=155, c_thick=3
	plot,SHIFT(time_freq, -mid_point_time), SHIFT(ABS(threedft_spatial_sum), -mid_point_time),/ylog,thick=3,$
             xtitle='Frequency (Hz)',ytitle='Power',title='Spatially averaged Fourier power',charsize=2,xst=1,yst=1,xr=[(-new_range_max),new_range_max]
        oplot,SHIFT(time_freq, -mid_point_time), SHIFT(ABS(temp_sum), -mid_point_time),thick=3,color=235
        oplot,SHIFT(time_freq, -mid_point_time), SHIFT(ABS(temp_sum), -mid_point_time),thick=3,color=235
	WAIT, 1.
	IF (NOT KEYWORD_SET(save_memory)) AND (NOT KEYWORD_SET(super_save_memory)) THEN filtered_cube = 1 ELSE datacube = 1 ; THIS JUST CLEARS THE MEMORY FROM THE ORIGINAL datacube ARRAY, WHICH IF NOT DONE CAN CAUSE MEMORY CRASHES!
        IF NOT KEYWORD_SET(super_save_memory) THEN filtered_fft_cube[*] = 0.
        ; USE THE USER-DEFINED SPATIAL FILTER TO ISOLATE SPATIAL POWER AT THESE SPATIAL FREQUENCIES
	IF NOT KEYWORD_SET(super_save_memory) THEN FOR t = 0, (N_ELEMENTS(threedft[0,0,*])-1) DO filtered_fft_cube[*,*,t] = REFORM(threedft[*,*,t]) * filter_map
	IF KEYWORD_SET(super_save_memory) THEN FOR t = 0, (N_ELEMENTS(threedft[0,0,*])-1) DO threedft[*,*,t] = REFORM(threedft[*,*,t]) * filter_map
	; NOW PERFORM FREQUENCY FILTERING USING THE USER-DEFINED PERIODS FOR THOSE ELEMENTS WITH NON-ZERO POWER
	FOR positive_values = 0, (N_ELEMENTS(values)-1) DO BEGIN
	    pixels = XYPOSITION(filter_map, values[positive_values])
	    pixel_x = pixels[0] & pixel_y = pixels[1]
	    IF NOT KEYWORD_SET(super_save_memory) THEN lc_fft = REFORM(filtered_fft_cube[pixel_x, pixel_y, *])
	    IF KEYWORD_SET(super_save_memory) THEN lc_fft = REFORM(threedft[pixel_x, pixel_y, *])
	    lc_fft[filter_range1] = 0. 
            lc_fft[filter_range2] = 0. 
            lc_fft[filter_range3] = 0. 
            IF NOT KEYWORD_SET(super_save_memory) THEN filtered_fft_cube[pixel_x, pixel_y, *] = lc_fft
	    IF KEYWORD_SET(super_save_memory) THEN threedft[pixel_x, pixel_y, *] = lc_fft
        ENDFOR 
	; FOR THE PURPOSES OF FILE NAMING, CALCULATE THE CHOSEN PERIOD MIN/MAX AND SPATIAL-SIZE MIN/MAX
	Pmin = ARR2STR(period_filters[z],/trim) 
        Pmax = ARR2STR(period_filters[z+1],/trim)
	remainder_of_spatial_size_min =  spatial_sizes[i] - FIX(spatial_sizes[i])
	addon_spatial_size_min = remainder_of_spatial_size_min * 100
	IF FIX(spatial_sizes[i]) ne spatial_sizes[i] THEN BEGIN
	    IF FIX(spatial_sizes[i]) le 9 THEN spatial_size_min = '00'+ARR2STR(FIX(spatial_sizes[i]),/trim)+'pt'+ARR2STR(FIX(addon_spatial_size_min),/trim)
	    IF (FIX(spatial_sizes[i]) gt 9) AND (FIX(spatial_sizes[i]) le 99) THEN spatial_size_min = '0'+ARR2STR(FIX(spatial_sizes[i]),/trim)+'pt'+ARR2STR(FIX(addon_spatial_size_min),/trim)
	    IF (FIX(spatial_sizes[i]) gt 99) AND (FIX(spatial_sizes[i]) le 999) THEN spatial_size_min = ARR2STR(FIX(spatial_sizes[i]),/trim)+'pt'+ARR2STR(FIX(addon_spatial_size_min),/trim)
	ENDIF
	IF FIX(spatial_sizes[i]) eq spatial_sizes[i] THEN BEGIN
	    IF FIX(spatial_sizes[i]) le 9 THEN spatial_size_min = '00'+ARR2STR(FIX(spatial_sizes[i]),/trim)+'pt00'
	    IF (FIX(spatial_sizes[i]) gt 9) AND (FIX(spatial_sizes[i]) le 99) THEN spatial_size_min = '0'+ARR2STR(FIX(spatial_sizes[i]),/trim)+'pt00'
	    IF (FIX(spatial_sizes[i]) gt 99) AND (FIX(spatial_sizes[i]) le 999) THEN spatial_size_min = ARR2STR(FIX(spatial_sizes[i]),/trim)+'pt00'
	ENDIF
	remainder_of_spatial_size_max =  spatial_sizes[i+1] - FIX(spatial_sizes[i+1])
	addon_spatial_size_max = remainder_of_spatial_size_max * 100
	IF FIX(spatial_sizes[i+1]) ne spatial_sizes[i+1] THEN BEGIN
	    IF FIX(spatial_sizes[i+1]) le 9 THEN spatial_size_max = '00'+ARR2STR(FIX(spatial_sizes[i+1]),/trim)+'pt'+ARR2STR(FIX(addon_spatial_size_max),/trim)
	    IF (FIX(spatial_sizes[i+1]) gt 9) AND (FIX(spatial_sizes[i+1]) le 99) THEN spatial_size_max = '0'+ARR2STR(FIX(spatial_sizes[i+1]),/trim)+'pt'+ARR2STR(FIX(addon_spatial_size_max),/trim)
	    IF (FIX(spatial_sizes[i+1]) gt 99) AND (FIX(spatial_sizes[i+1]) le 999) THEN spatial_size_max = ARR2STR(FIX(spatial_sizes[i+1]),/trim)+'pt'+ARR2STR(FIX(addon_spatial_size_max),/trim)
	ENDIF
	IF FIX(spatial_sizes[i+1]) eq spatial_sizes[i+1] THEN BEGIN
	    IF FIX(spatial_sizes[i+1]) le 9 THEN spatial_size_max = '00'+ARR2STR(FIX(spatial_sizes[i+1]),/trim)+'pt00'
	    IF (FIX(spatial_sizes[i+1]) gt 9) AND (FIX(spatial_sizes[i]) le 99) THEN spatial_size_max = '0'+ARR2STR(FIX(spatial_sizes[i+1]),/trim)+'pt00'
	    IF (FIX(spatial_sizes[i+1]) gt 99) AND (FIX(spatial_sizes[i]) le 999) THEN spatial_size_max = ARR2STR(FIX(spatial_sizes[i+1]),/trim)+'pt00'
	ENDIF
	; TAKE THE FILTERED FFT CUBE AND REFORM BACK INTO [x, y, t]
	IF (NOT KEYWORD_SET(save_memory)) AND (NOT KEYWORD_SET(super_save_memory)) THEN filtered_cube = REAL_PART(FFT(filtered_fft_cube, 1)) 
	IF KEYWORD_SET(save_memory) THEN datacube = REAL_PART(FFT(filtered_fft_cube, 1))
	; TO SAVE ON MEMORY, THE THREEDFT VARIABLE IS REUSED ALONGSIDE THE temporary COMMAND TO ASSIST WITH MEMORY REUSE
	IF KEYWORD_SET(super_save_memory) THEN threedft = REAL_PART(FFT(TEMPORARY(threedft), 1))
        ; SAVE THE FILTERED ARRAYS
	IF (NOT KEYWORD_SET(save_memory)) AND (NOT KEYWORD_SET(super_save_memory)) THEN WRITEFITS,save_directory + 'Filtered_datacube_' + Pmin + 'Pmin_' + Pmax + 'Pmax_' + spatial_size_min + 'arcsecMin_' + spatial_size_max + 'arcsecMax.fits',filtered_cube 
        IF KEYWORD_SET(save_memory) THEN WRITEFITS,save_directory + 'Filtered_datacube_' + Pmin + 'Pmin_' + Pmax + 'Pmax_' + spatial_size_min + 'arcsecMin_' + spatial_size_max + 'arcsecMax.fits',datacube 
        IF KEYWORD_SET(super_save_memory) THEN WRITEFITS,save_directory + 'Filtered_datacube_' + Pmin + 'Pmin_' + Pmax + 'Pmax_' + spatial_size_min + 'arcsecMin_' + spatial_size_max + 'arcsecMax.fits',threedft 
    ENDFOR
    ; ESTIMATE THE REMAINING PROCESSING TIME
    countup = countup + 1. 
    systm = SYSTIME(1) 
    loop_time = (systm - t1) 
    t1 = systm 
    minutes = (FIX(((loop_time)*(period_bins - countup))/60.)) 
    hours = FIX(minutes/60.) 
    hourminutes = hours*60. 
    remain_minutes = minutes - hourminutes 
    print,'I am '+ARR2STR(((countup/(period_bins))*100.),/trim) + $
          ' percent complete with approximately '+ARR2STR(hours,/trim)+' hours and '+$
	  ARR2STR(FIX(remain_minutes),/trim)+' minutes left at current CPU speed' 
ENDFOR

mult,1,1
loadct,0,/silent

END
