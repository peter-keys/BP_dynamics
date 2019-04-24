;	Working out the Energy/Wavenumbers/Theoretically whether they're surface/body
;---------------------------------------------------------------------------------------
; Following approaches laid out by Moreels et al 2013, Grant et al. 2015 etc.
;
; Need to work out the following in order:
;	Va	-	Alfven Speed
;	Cs	-	Sound Speed
;	A	-	Ratio of Intensities & areas
;	Vph	-	Phase Speed
;	k	-	Wavenumber
;	Ct	-	Tube speed
;	Kappa_squared	-	From Moreels for S/B clarification
;	=> S/B mode
;
; NB G-band formation Height is taken as 100km (Jess et al. 2012)
;CHANGE TO VAL A
;--------------------------------------------------------------------------------------
;Alfven speed:

mu_naught = 1.25663706E-6		;Permeability of free space (constant)
rho = 1.604E-4				;Typical density (VAL-A Model)
Bz = 1100 * 0.0001			;Typical Line of sight B component (HMI data) in Tesla

Va = Bz / (SQRT(mu_naught * rho))	;In m/s
;--------------------------------------------------------------------------------------
;Sound speed:

gamma = 5/3.				;Ratio of specific Heats
R = 8.3144621				;Gas Constant
T = 5450				;Gas temp. in K (from VAL-A)
mu_hat = (1.609E-7) / 6.86E16		;Av. Mass per particle (from VAL-A)
mu = mu_hat / 1.67262158E-27		;Mean molecular Weight

Cs = SQRT((gamma * R * T) / mu)		;Might need to check that at the min 
Cs = Cs * 1000.				;In m/s
;--------------------------------------------------------------------------------------
;Tube Speed:

Ct = (Cs * Va)/ (SQRT((Cs*Cs) + (Va*Va)))
;--------------------------------------------------------------------------------------
;Phase Speed:

;10Dec2011 (1200G):			09Dec2011 (1200G):		30Sept2012(820G):		06Mar2013 (850G):
;	  +	460s = 2.17mHz		+	480s = 2.08mHz		+	360s = 2.78mHz		+	280s = 3.57mHZ
;	  +	216s = 4.90mHz		+	260s = 3.85mHz		+	220s = 4.55mHz		+	145s = 6.9mHz
;	  +	145s = 6.90mHz		+	135s = 7.41mHz		+	90s = 11.11mHz		+	85s = 11.76mHz

;17Aug2013 (990G):			15Apr2014 (1100G):		11Jul2011 (1300):
;	  +	400s = 2.5mHz		+	350s = 2.86mHz		+	500s = 2.00mHz
;	  +	245s = 4.08mHz		+	190s = 5.26mHz		+	230s = 4.35mHz
;	  +	140s = 7.14mHz		+	100s = 10mHz		+	105s = 9.5mHz

freq = 3E-8/(430.55*1E-9)				;Freq. of Observation
h = 6.62606957E-34					;Planck Constant
kb = 1.3806488E-23					;Boltzmann Constant
;LOAD IN THE INTENSITY VALUES HERE PLEASE!
;restore,'/data/rosa3/oldrosa1/phk/data/wavelet_datasets/11Jul2012/Gband/key_parameters_gband_juldata.sav',/ver		;11Jul2011
;restore,'/data/rosa3/oldrosa1/phk/data/wavelet_datasets/09Dec2011/Gband/key_parameters_gband_09decdata.sav',/ver	;09Dec2011
;restore,'/data/rosa3/oldrosa1/phk/data/wavelet_datasets/10Dec2011/Gband/key_parameters_gband_10decdata.sav',/ver	;10Dec2011
;restore,'/data/rosa3/oldrosa1/phk/data/wavelet_datasets/30Sept2012/key_parameters_septdata.sav',/ver			;30Sept2012
;restore,'/data/rosa3/oldrosa1/phk/data/wavelet_datasets/06Mar2013/key_parameters_mardata.sav',/ver			;06Mar2013
;restore,'/data/rosa3/oldrosa1/phk/data/wavelet_datasets/17Aug2013/Gband/key_parameters_gband_augdata.sav',/ver		;17Aug2013
restore,'/data/rosa3/oldrosa1/phk/data/wavelet_datasets/Apr2014/key_parameters_gband_aprdata.sav',/ver			;15Apr2014

delta_I = (MAX(marporeinten)-MIN(marporeinten))/2.	;Intensity Amplitude
I_naught = MEAN(marporeinten)				;Intensity Mean
delta_S = (MAX(marporeareas)-MIN(marporeareas))/2.	;Area Amplitude
S_naught = MEAN(marporeareas)				;Area Mean

A = (delta_I / I_naught) / (delta_s / S_naught)

Vph_top = A - 1
Vph_bottom_1 = (gamma - 1) * ((h * freq) / (kb*T))
Vph_bottom = (A - 1) + Vph_bottom_1
Vph_half = SQRT(Vph_top / Vph_bottom)
Vph = Cs * Vph_half

Vph = Ct
;Vph = Cs * (SQRT((A - 1) / ((A - 1) + ((gamma - 1)*((h * freq) / (kb * T))))))
;--------------------------------------------------------------------------------------
;Longitudinal Wavenumber:

;10Dec2011 (1200G):			09Dec2011 (1200G):		30Sept2012(820G):		06Mar2013 (850G):
;	  +	460s = 2.17mHz		+	480s = 2.08mHz		+	360s = 2.78mHz		+	280s = 3.57mHZ
;	  +	216s = 4.90mHz		+	260s = 3.85mHz		+	220s = 4.55mHz		+	145s = 6.9mHz
;	  +	145s = 6.90mHz		+	135s = 7.41mHz		+	90s = 11.11mHz		+	85s = 11.76mHz

;17Aug2013 (990G):			15Apr2014 (1100G):		11Jul2011 (1300):
;	  +	400s = 2.5mHz		+	350s = 2.86mHz		+	500s = 2.00mHz
;	  +	245s = 4.08mHz		+	190s = 5.26mHz		+	230s = 4.35mHz
;	  +	140s = 7.14mHz		+	100s = 10mHz		+	105s = 9.5mHz

P = 100.				;Period of observed oscillation
k = (2 * !PI) / (P * Vph)		;in m^-1
k_Mm = k/(1./1000000.)			;in Mm^-1

print,k_Mm
;--------------------------------------------------------------------------------------
;Wavelength:

lambda = (2 * !PI)/k		;In m
lambda_km = lambda / 1000.	;In km
;--------------------------------------------------------------------------------------
;Tube Speed:

Ct = (Cs * Va)/ (SQRT((Cs*Cs) + (Va*Va)))
;--------------------------------------------------------------------------------------
;Kappa_squared:

omega = (2 * !PI) / P		;The angular freq. of the oscillation
kappa_top = (k*k*Cs*Cs - omega*omega) * (k*k*Va*Va - omega*omega)
kappa_bottom = (Cs*Cs + Va*Va) * (k*k*Ct*Ct - omega*omega)

kappa_squared = kappa_top / kappa_bottom

print,k,lambda_km,omega,kappa_squared

print,va/1000.,vph/1000.,cs/1000., A,Ct/1000.
;--------------------------------------------------------------------------------------
;Surface/Body?

print,'Kappa^2 = '+arr2str(kappa_squared,/trim)
print,' '
print,'+VE is Surface   ;   -VE is Body'
print,' '
;--------------------------------------------------------------------------------------

