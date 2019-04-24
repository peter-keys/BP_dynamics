;-------------------------------------------------------------------------------------------------------------
;		ENERGY CALC USING MOREELS et al. 2015
;-------------------------------------------------------------------------------------------------------------
; Really convoluted....
;---------------------------------------------------------------------------------------
; Following approaches laid out by Moreels et al 2013/2015, Grant et al. 2015 etc.
;
; Need to work out the following:
;	Va	-	Alfven Speed
;	Cs	-	Sound Speed
;	A	-	Ratio of Intensities & areas
;	Vph	-	Phase Speed
;	k	-	Wavenumber
;	Ct	-	Tube speed
;
; NB G-band formation Height is taken as 100km (Jess et al. 2012)
;--------------------------------------------------------------------------------------
;First State Some Values for your wave (from obs):
;Reminder of the Freq's.... (Chromological order)
;11Jul2011 (1300):			09Dec2011 (1200G):		10Dec2011 (1200G):
;	+	500s = 2.00mHz		+	480s = 2.08mHz	  	+	460s = 2.17mHz
;	+	230s = 4.35mHz		+	260s = 3.85mHz	  	+	216s = 4.90mHz
;	+	105s = 9.5mHz		+	135s = 7.41mHz	  	+	145s = 6.90mHz

;30Sept2012(820G):			06Mar2013 (850G):		17Aug2013 (990G):
;	  +	360s = 2.78mHz		+	280s = 3.57mHZ		  +	400s = 2.5mHz	
;	  +	220s = 4.55mHz		+	145s = 6.9mHz		  +	245s = 4.08mHz
;	  +	90s = 11.11mHz		+	85s = 11.76mHz	 	  +	140s = 7.14mHz

;15Apr2014 (1100G):		 
;	  +	350s = 2.86mHz		
;	  +	190s = 5.26mHz		
;	  +	100s = 10mHz	
;--------------------------------------------------------------------------------------
;State the period and the LOS B component for the wave (from above)

P = 105.				;in secs
Bz = 1300.				;in Gauss (HMI) 
Bz = Bz * 0.0001			;Change to Tesla	
omega = (2 * !PI) / P			;The angular freq. of the oscillation	
;--------------------------------------------------------------------------------------
;Alfven speed:

mu_naught = 1.25663706E-6		;Permeability of free space (constant)
;rho = 7.5E-5				;Typical density (VAL-A Model)
rho = 4.E-5
freq = 3E-8/(430.55*1E-9)		;Freq. of Observation
h = 6.62606957E-34			;Planck Constant
kb = 1.3806488E-23			;Boltzmann Constant

Va_alt = Bz / (SQRT(mu_naught * rho))	;In m/s
Va = 12000.
;--------------------------------------------------------------------------------------
;Sound speed:

gamma = 5/3.				;Ratio of specific Heats
R = 8.3144621				;Gas Constant
T = 5450				;Gas temp. in K (from VAL-A)
mu_hat = (1.609E-7) / 6.86E16		;Av. Mass per particle (from VAL-A)
mu = mu_hat / 1.67262158E-27		;Mean molecular Weight

Cs = SQRT((gamma * R * T) / mu)		;Might need to check that at the min 
Cs = Cs * 1000.				;In m/s
;cs = 6400.
;--------------------------------------------------------------------------------------
;Tube Speed:

Ct = (Cs * Va)/ (SQRT((Cs*Cs) + (Va*Va)))
;ct=5000.
;--------------------------------------------------------------------------------------
;Phase Speed:

Vph = Ct

;NB Don't have info. at multiple heights, therfore, Vph is estimated as the tube speed.
;--------------------------------------------------------------------------------------
;Wavenumber (k):

k = (2 * !PI) / (P * Vph)		;in m^-1
k_Mm = k/(1./1000000.)			;in Mm^-1 (useful)
;--------------------------------------------------------------------------------------
;Calculate the Dimensionless variable A:
;LOAD IN THE INTENSITY VALUES HERE PLEASE!
;restore,'/data/rosa3/oldrosa1/phk/data/wavelet_datasets/11Jul2012/Gband/key_parameters_gband_juldata.sav',/ver		;11Jul2011
;restore,'/data/rosa3/oldrosa1/phk/data/wavelet_datasets/09Dec2011/Gband/key_parameters_gband_09decdata.sav',/ver	;09Dec2011
;restore,'/data/rosa3/oldrosa1/phk/data/wavelet_datasets/10Dec2011/Gband/key_parameters_gband_10decdata.sav',/ver	;10Dec2011
;restore,'/data/rosa3/oldrosa1/phk/data/wavelet_datasets/30Sept2012/key_parameters_septdata.sav',/ver			;30Sept2012
;restore,'/data/rosa3/oldrosa1/phk/data/wavelet_datasets/06Mar2013/key_parameters_mardata.sav',/ver			;06Mar2013
;restore,'/data/rosa3/oldrosa1/phk/data/wavelet_datasets/17Aug2013/Gband/key_parameters_gband_augdata.sav',/ver		;17Aug2013
;restore,'/data/rosa3/oldrosa1/phk/data/wavelet_datasets/Apr2014/key_parameters_gband_aprdata.sav',/ver			;15Apr2014

;delta_I = (MAX(marporeinten)-MIN(marporeinten))/2.	;Intensity Amplitude
;I_naught = MEAN(marporeinten)				;Intensity Mean
;delta_S = (MAX(marporeareas)-MIN(marporeareas))/2.	;Area Amplitude
;S_naught = MEAN(marporeareas)				;Area Mean

;A = (delta_I / I_naught) / (delta_s / S_naught)
;--------------------------------------------------------------------------------------
;Take deltaI/I and deltaa/A values from EMD results (comment out as appropriate):
;------------
;11Jul2011:
;------------
;1st:
;percent_I = 1.55/100.
;percent_A = 1.06/100.
;2nd
;percent_I = 1.23/100.
;percent_A = 1.71/100.
;3rd
;percent_I = 1.63/100.
;percent_A = 2.045/100.
;------------
;avspercent_I = 1.47/100.
;avspercent_A = 1.61/100.
;------------
;09Dec2011
;------------
;1st:
;percent_I = 1.899/100.
;percent_A = 2.048/100.
;2nd
;percent_I = 2.014/100.
;percent_A = 1.612/100.
;3rd
;percent_I = 2.571/100.
;percent_A = 0.990/100.
;------------
;avspercent_I = 2.16/100.
;avspercent_A = 1.72/100.
;------------
;10Dec2011
;------------
;1st:
;percent_I = 2.649/100.
;percent_A = 4.274/100.
;2nd
;percent_I = 1.821/100.
;percent_A = 4.757/100.
;3rd
;percent_I = 1.821/100.
;percent_A = 2.992/100.
;------------
;avspercent_I = 2.1/100.
;avspercent_A = 4.01/100.
;------------
;30Sept2012
;------------
;1st:
;percent_I = 2.916/100.
;percent_A = 2.869/100.
;2nd
;percent_I = 3.296/100.
;percent_A = 2.032/100.
;3rd
;percent_I = 3.249/100.
;percent_A = 2.04/100.
;------------
;avspercent_I = 3.15/100.
;avspercent_A = 2.45/100.
;------------
;06Mar2013
;------------
;1st:
;percent_I = 2.132/100.
;percent_A = 2.727/100.
;2nd
;percent_I = 1.379/100.
;percent_A = 2.727/100.
;3rd
;percent_I = 2.2/100.
;percent_A = 2.934/100.
;------------
;avspercent_I = 1.9/100.
;avspercent_A = 2.8/100.
;------------
;17Aug2013
;------------
;1st:
;percent_I = 1.495/100.
;percent_A = 1.281/100.
;2nd
;percent_I = 1.055/100.
;percent_A = 2.218/100.
;3rd
;percent_I = 0.719/100.
;percent_A = 2.571/100.
;------------
;avspercent_I = 1.09/100.
;avspercent_A = 2.02/100.
;------------
;15Apr2014
;------------
;1st:
percent_I = 2.412/100.
percent_A = 1.805/100.
;2nd
;percent_I = 2.795/100.
;percent_A = 2.370/100.
;3rd
;percent_I = 2.412/100.
;percent_A = 2.476/100.
;------------
;avspercent_I = 2.54/100.
;avspercent_A = 2.22/100.
;------------

A = percent_I/percent_A

;pixel oscillation conversion to area osc:
;print,SQRT((184.6*((0.16/2.)*725000.)^2)/!PI)
;vph =2200.
;Vph = Cs * (SQRT((-1)*(A - 1) / ((-1)*(A - 1) + ((gamma - 1)*((h * freq) / (kb * T))))))
;--------------------------------------------------------------------------------------
;Calculate the m variable

;Moreels:
m_top = 	((k*k*Cs*Cs) - (omega*omega)) * ((k*k*Va*Va) - (omega*omega))
m_bottom = 	    ((Cs*Cs) + (Va*Va)) *  ((k*k*Ct*Ct) - (omega*omega))
;m_top = 		k*k*((k*k*Cs*Cs) - (omega*omega)) * ((k*k*Va*Va) - (omega*omega))
;m_bottom = 			((k*k*Cs*Cs) + (k*k*Va*Va)) * ((k*k*Ct*Ct) - (omega*omega))

;Morton:
m_top_morton = 	((omega*omega)-((Va*Va)*(k*k))) * ((omega*omega)-((Cs*Cs)*(k*k)))
m_bottom_morton = 		((Cs*Cs)+(Va*Va)) * ((omega*omega) - ((Ct*Ct)*(k*k)))

m_morton = m_top_morton/m_bottom_morton

m = m_top/m_bottom
;--------------------------------------------------------------------------------------
;Work out the Lagrangian displacment variable 'funk'

;funk = (SQRT(ABS(m))) / (rho * ((omega*omega) - ((Va*Va)*(k*k))))
;funk = funk * A
;funk = (1.72/100)*((5.6*1e6)/2.)

;Lagrangian displacement worked out with the %change and average diameters:
;The Average diameters are:
;11Jul2011 (1300):			09Dec2011 (1200G):		10Dec2011 (1200G):
;	+	500s =		+	480s = 	  	+	460s = 
;	+	230s = 		+	260s = 	  	+	216s = 
;	+	105s = 		+	135s = 	  	+	145s = 

;30Sept2012(820G):			06Mar2013 (850G):		17Aug2013 (990G):
;	  +	360s = 		+	280s = 		  +	400s = 	
;	  +	220s = 		+	145s = 		  +	245s = 
;	  +	90s = 		+	85s = 	 	  +	140s = 

;15Apr2014 (1100G):		 
;	  +	350s = 		
;	  +	190s = 		
;	  +	100s = 
;	
;11Jul2011: d = 8.6Mm  
;r = 4.3E6

;09Dec2011: d = 5.6 Mm 
;r = 2.8E6

;10Dec2011: d = 7.7Mm 
;r = 3.85E6

;30Sept2012: d = 1.5Mm 
;r = 7.5E5

;06Mar2013: d = 1.3Mm 
;r = 6.5E5

;17Aug2013: d = 3.0Mm 
;r = 1.5E6

;15Apr2014: d = 6.3Mm 
r = 3.15E6

deltar = r*percent_A

;--------------------------------------------------------------------------------------
;Work out the energy using these....

f = 0.1				;Filling factor...

E = f * (1 + ALOG(1/f))*(rho/2.) * (omega*omega) * (deltar^2) * Va
E_kW = E/1000.			;in kW/m2

print,'--------------------------------------------------------------------------------------'
print,'omega,Va,cs,ct,A,Vph,k,k_Mm, Va_alt'
print,omega,Va,cs,ct,A,Vph,k,k_Mm, Va_alt
print,'--------------------------------------------------------------------------------------'
print,'m,m_morton,deltar,E_kW'
print,m,m_morton,deltar,E_kW
print,'--------------------------------------------------------------------------------------'

;--------------------------------------------------------------------------------------
vz = 1000.
E_alt = f * (rho/2.) * (Vz^2) * Ct

;Morton Approach:

;Work out Vr:

;Vr = -1* A * ((omega^2 - (k^2 * Cs^2))/(m * omega^2))

;Vz = -1 * A * ((k * cs^2)/(omega^2))

;Br = -1 * (k / omega) * Bz * Vr

;Bbz = A * ((omega^2 - (k^2 * Cs^2))/(omega^3))*Bz


;E = 0.25 * Vph * ((Vr^2 + Vz^2)+ (Br^2 + Bbz^2))
 
;E = 0.25 * Vph * ((rho * vr^2)+(Br^2)/mu_naught)
