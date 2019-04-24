;==================================================================
;	STITCHING BACK THE SUB-FIELDS GENERATED TO SPEED UP 
;			INVERSIONS
;==================================================================
;23/3/2018
; This is the process for taking the two arrays 
; of an MBP sub-field over several frames, 
; then sticthing them back together to make it easier to analyse

; Output Directories of inverted data (1 is 1st batch chronologically 2 is 2nd...)
diro1 = '/data/solarstore3/phk/nicole_2018/NICOLE/run/out_MBP35/r3/outputs/'
diro2 = '/data/solarstore3/phk/nicole_2018/test_environment/NICOLE2/run/out_MBP35/r3/outputs/'
; Input directories
diri1 = '/data/solarstore3/phk/nicole_2018/NICOLE/run/out_MBP35/r3/inputs/'
diri2 = '/data/solarstore3/phk/nicole_2018/test_environment/NICOLE2/run/out_MBP35/r3/inputs/'

; Best Model directories
dirm1 = '/data/solarstore3/phk/nicole_2018/NICOLE/run/out_MBP35/r1/outputs/'
dirm2 = '/data/solarstore3/phk/nicole_2018/test_environment/NICOLE2/run/out_MBP35/r1/outputs/'

s1i = READ_PROFILE(diro1+'synthetic.pro',s1q,s1u,s1v)
s2i = READ_PROFILE(diro2+'synthetic.pro',s2q,s2u,s2v)

o1i = READ_PROFILE(diri1+'observed_MBP35.nic',o1q,o1u,o1v)
o2i = READ_PROFILE(diri2+'observed_MBP35.nic',o2q,o2u,o2v)

; Good Models
mg1 = READ_MODEL(dirm1+'modelout.mod')
mg2 = READ_MODEL(dirm2+'modelout.mod')

;Final Models
mf1 = READ_MODEL(diro1+'modelout.mod')
mf2 = READ_MODEL(diro2+'modelout.mod')

;Dimensions of 1st sub-field 
RESTORE,diri1+'MBP35_xyf_info.sav',/ver
xsize1 = x1 - x0 + 1
nims1 = nx1/xsize1


;Dimensions of 2nd sub-field 
RESTORE,diri2+'MBP35_xyf_info.sav',/ver
xsize2 = x1 - x0 + 1
nims2 = nx1/xsize2

IF (nims1+nims2 NE nf1) THEN print,'Number of images dont match up...'


; Make arrays of the new I and V Observed and Synhtetic profiles for the sub-field... (don't need Q and U)
oi_new = FLTARR(xsize1,ny1,nf1,92)
ov_new = FLTARR(xsize1,ny1,nf1,92)

FOR i = 0, nims1-1 DO BEGIN &$
	oi_new[*,*,i,*] = o1i[xsize1*i:(xsize1*(i+1))-1,*,*] &$
	ov_new[*,*,i,*] = o1v[xsize1*i:(xsize1*(i+1))-1,*,*] &$
ENDFOR
FOR i = nims1, nf1-1 DO BEGIN &$
	imnum = i - nims1 &$
	oi_new[*,*,i,*] = o2i[xsize2*imnum:(xsize2*(imnum+1))-1,*,*] &$
	ov_new[*,*,i,*] = o2v[xsize2*imnum:(xsize2*(imnum+1))-1,*,*] &$
ENDFOR
	

si_new = FLTARR(xsize1,ny1,nf1,92)
sv_new = FLTARR(xsize1,ny1,nf1,92)

FOR i = 0, nims1-1 DO BEGIN &$
	si_new[*,*,i,*] = s1i[xsize1*i:(xsize1*(i+1))-1,*,*] &$
	sv_new[*,*,i,*] = s1v[xsize1*i:(xsize1*(i+1))-1,*,*] &$
ENDFOR
FOR i = nims1, nf1-1 DO BEGIN &$
	imnum = i - nims1 &$
	si_new[*,*,i,*] = s2i[xsize2*imnum:(xsize2*(imnum+1))-1,*,*] &$
	sv_new[*,*,i,*] = s2v[xsize2*imnum:(xsize2*(imnum+1))-1,*,*] &$
ENDFOR

SAVE,FILENAME=diro1+'Obs_final_MBP35_series.sav',oi_new,ov_new
SAVE,FILENAME=diro1+'Syn_final_MBP35_series.sav',si_new,sv_new

; Make arrays for the Good Models...
; Take specific Tau Values (-2, -1 & 0...)
; For MBP35 tau = 0.025 occurs at 85
;	    tau = -0.994 occurs at 73
; 	    tau = -2.014 occurs at 61

; Make arrays for the interesting things...
B_tau0 = FLTARR(xsize1,ny1,nf1)		
B_tau1 = FLTARR(xsize1,ny1,nf1)		;B_LOS_Z at tau = -1
B_tau2 = FLTARR(xsize1,ny1,nf1)

T_tau0 = FLTARR(xsize1,ny1,nf1)		
T_tau1 = FLTARR(xsize1,ny1,nf1)		;Temp at tau = -1
T_tau2 = FLTARR(xsize1,ny1,nf1)

V_tau0 = FLTARR(xsize1,ny1,nf1)		
V_tau1 = FLTARR(xsize1,ny1,nf1)		;Vel_LOS at tau = -1
V_tau2 = FLTARR(xsize1,ny1,nf1)

ElP_tau0 = FLTARR(xsize1,ny1,nf1)		
ElP_tau1 = FLTARR(xsize1,ny1,nf1)	;El_P at tau = -1
ElP_tau2 = FLTARR(xsize1,ny1,nf1)

GasP_tau0 = FLTARR(xsize1,ny1,nf1)		
GasP_tau1 = FLTARR(xsize1,ny1,nf1)	;Gas_P at tau = -1
GasP_tau2 = FLTARR(xsize1,ny1,nf1)

FOR i = 0, nims1-1 DO BEGIN &$
	B_tau0[*,*,i] = mg1.B_LOS_Z[xsize1*i:(xsize1*(i+1))-1,*,85] &$
	B_tau1[*,*,i] = mg1.B_LOS_Z[xsize1*i:(xsize1*(i+1))-1,*,73] &$
	B_tau2[*,*,i] = mg1.B_LOS_Z[xsize1*i:(xsize1*(i+1))-1,*,61] &$
	
	T_tau0[*,*,i] = mg1.T[xsize1*i:(xsize1*(i+1))-1,*,85] &$
	T_tau1[*,*,i] = mg1.T[xsize1*i:(xsize1*(i+1))-1,*,73] &$
	T_tau2[*,*,i] = mg1.T[xsize1*i:(xsize1*(i+1))-1,*,61] &$
	
	V_tau0[*,*,i] = mg1.v_LOS[xsize1*i:(xsize1*(i+1))-1,*,85] &$
	V_tau1[*,*,i] = mg1.v_LOS[xsize1*i:(xsize1*(i+1))-1,*,73] &$
	V_tau2[*,*,i] = mg1.v_LOS[xsize1*i:(xsize1*(i+1))-1,*,61] &$

	ELP_tau0[*,*,i] = mg1.EL_P[xsize1*i:(xsize1*(i+1))-1,*,85] &$
	ELP_tau1[*,*,i] = mg1.EL_P[xsize1*i:(xsize1*(i+1))-1,*,73] &$
	ELP_tau2[*,*,i] = mg1.EL_P[xsize1*i:(xsize1*(i+1))-1,*,61] &$

	GasP_tau0[*,*,i] = mg1.Gas_P[xsize1*i:(xsize1*(i+1))-1,*,85] &$
	GasP_tau1[*,*,i] = mg1.Gas_P[xsize1*i:(xsize1*(i+1))-1,*,73] &$
	GasP_tau2[*,*,i] = mg1.Gas_P[xsize1*i:(xsize1*(i+1))-1,*,61] &$
ENDFOR
FOR i = nims1, nf1-1 DO BEGIN &$
	imnum = i - nims1 &$
	B_tau0[*,*,i] = mg2.B_LOS_Z[xsize1*imnum:(xsize1*(imnum+1))-1,*,85] &$
	B_tau1[*,*,i] = mg2.B_LOS_Z[xsize1*imnum:(xsize1*(imnum+1))-1,*,73] &$
	B_tau2[*,*,i] = mg2.B_LOS_Z[xsize1*imnum:(xsize1*(imnum+1))-1,*,61] &$
	
	T_tau0[*,*,i] = mg2.T[xsize1*imnum:(xsize1*(imnum+1))-1,*,85] &$
	T_tau1[*,*,i] = mg2.T[xsize1*imnum:(xsize1*(imnum+1))-1,*,73] &$
	T_tau2[*,*,i] = mg2.T[xsize1*imnum:(xsize1*(imnum+1))-1,*,61] &$
	
	V_tau0[*,*,i] = mg2.v_LOS[xsize1*imnum:(xsize1*(imnum+1))-1,*,85] &$
	V_tau1[*,*,i] = mg2.v_LOS[xsize1*imnum:(xsize1*(imnum+1))-1,*,73] &$
	V_tau2[*,*,i] = mg2.v_LOS[xsize1*imnum:(xsize1*(imnum+1))-1,*,61] &$

	ELP_tau0[*,*,i] = mg2.EL_P[xsize1*imnum:(xsize1*(imnum+1))-1,*,85] &$
	ELP_tau1[*,*,i] = mg2.EL_P[xsize1*imnum:(xsize1*(imnum+1))-1,*,73] &$
	ELP_tau2[*,*,i] = mg2.EL_P[xsize1*imnum:(xsize1*(imnum+1))-1,*,61] &$

	GasP_tau0[*,*,i] = mg2.Gas_P[xsize1*imnum:(xsize1*(imnum+1))-1,*,85] &$
	GasP_tau1[*,*,i] = mg2.Gas_P[xsize1*imnum:(xsize1*(imnum+1))-1,*,73] &$
	GasP_tau2[*,*,i] = mg2.Gas_P[xsize1*imnum:(xsize1*(imnum+1))-1,*,61] &$
ENDFOR
SAVE,FILENAME=diro1+'GModel_final_MBP35_series_tau0.sav',B_tau0,T_tau0,V_tau0,ELP_tau0,GasP_tau0
SAVE,FILENAME=diro1+'GModel_final_MBP35_series_taum1.sav',B_tau1,T_tau1,V_tau1,ELP_tau1,GasP_tau1
SAVE,FILENAME=diro1+'GModel_final_MBP35_series_taum2.sav',B_tau2,T_tau2,V_tau2,ELP_tau2,GasP_tau2


; Make arrays for the Final Models... (not necessarily good...)
; Take specific Tau Values (-2, -1 & 0...)
; For MBP35 tau = 0.025 occurs at 85
;	    tau = -0.994 occurs at 73
; 	    tau = -2.014 occurs at 61

; Make arrays for the interesting things...
B_tau0 = FLTARR(xsize1,ny1,nf1)		
B_tau1 = FLTARR(xsize1,ny1,nf1)		;B_LOS_Z at tau = -1
B_tau2 = FLTARR(xsize1,ny1,nf1)

T_tau0 = FLTARR(xsize1,ny1,nf1)		
T_tau1 = FLTARR(xsize1,ny1,nf1)		;Temp at tau = -1
T_tau2 = FLTARR(xsize1,ny1,nf1)

V_tau0 = FLTARR(xsize1,ny1,nf1)		
V_tau1 = FLTARR(xsize1,ny1,nf1)		;Vel_LOS at tau = -1
V_tau2 = FLTARR(xsize1,ny1,nf1)

ElP_tau0 = FLTARR(xsize1,ny1,nf1)		
ElP_tau1 = FLTARR(xsize1,ny1,nf1)	;El_P at tau = -1
ElP_tau2 = FLTARR(xsize1,ny1,nf1)

GasP_tau0 = FLTARR(xsize1,ny1,nf1)		
GasP_tau1 = FLTARR(xsize1,ny1,nf1)	;Gas_P at tau = -1
GasP_tau2 = FLTARR(xsize1,ny1,nf1)

FOR i = 0, nims1-1 DO BEGIN &$
	B_tau0[*,*,i] = mf1.B_LOS_Z[xsize1*i:(xsize1*(i+1))-1,*,85] &$
	B_tau1[*,*,i] = mf1.B_LOS_Z[xsize1*i:(xsize1*(i+1))-1,*,73] &$
	B_tau2[*,*,i] = mf1.B_LOS_Z[xsize1*i:(xsize1*(i+1))-1,*,61] &$
	
	T_tau0[*,*,i] = mf1.T[xsize1*i:(xsize1*(i+1))-1,*,85] &$
	T_tau1[*,*,i] = mf1.T[xsize1*i:(xsize1*(i+1))-1,*,73] &$
	T_tau2[*,*,i] = mf1.T[xsize1*i:(xsize1*(i+1))-1,*,61] &$
	
	V_tau0[*,*,i] = mf1.v_LOS[xsize1*i:(xsize1*(i+1))-1,*,85] &$
	V_tau1[*,*,i] = mf1.v_LOS[xsize1*i:(xsize1*(i+1))-1,*,73] &$
	V_tau2[*,*,i] = mf1.v_LOS[xsize1*i:(xsize1*(i+1))-1,*,61] &$

	ELP_tau0[*,*,i] = mf1.EL_P[xsize1*i:(xsize1*(i+1))-1,*,85] &$
	ELP_tau1[*,*,i] = mf1.EL_P[xsize1*i:(xsize1*(i+1))-1,*,73] &$
	ELP_tau2[*,*,i] = mf1.EL_P[xsize1*i:(xsize1*(i+1))-1,*,61] &$

	GasP_tau0[*,*,i] = mf1.Gas_P[xsize1*i:(xsize1*(i+1))-1,*,85] &$
	GasP_tau1[*,*,i] = mf1.Gas_P[xsize1*i:(xsize1*(i+1))-1,*,73] &$
	GasP_tau2[*,*,i] = mf1.Gas_P[xsize1*i:(xsize1*(i+1))-1,*,61] &$
ENDFOR
FOR i = nims1, nf1-1 DO BEGIN &$
	imnum = i - nims1 &$
	B_tau0[*,*,i] = mf2.B_LOS_Z[xsize1*imnum:(xsize1*(imnum+1))-1,*,85] &$
	B_tau1[*,*,i] = mf2.B_LOS_Z[xsize1*imnum:(xsize1*(imnum+1))-1,*,73] &$
	B_tau2[*,*,i] = mf2.B_LOS_Z[xsize1*imnum:(xsize1*(imnum+1))-1,*,61] &$
	
	T_tau0[*,*,i] = mf2.T[xsize1*imnum:(xsize1*(imnum+1))-1,*,85] &$
	T_tau1[*,*,i] = mf2.T[xsize1*imnum:(xsize1*(imnum+1))-1,*,73] &$
	T_tau2[*,*,i] = mf2.T[xsize1*imnum:(xsize1*(imnum+1))-1,*,61] &$
	
	V_tau0[*,*,i] = mf2.v_LOS[xsize1*imnum:(xsize1*(imnum+1))-1,*,85] &$
	V_tau1[*,*,i] = mf2.v_LOS[xsize1*imnum:(xsize1*(imnum+1))-1,*,73] &$
	V_tau2[*,*,i] = mf2.v_LOS[xsize1*imnum:(xsize1*(imnum+1))-1,*,61] &$

	ELP_tau0[*,*,i] = mf2.EL_P[xsize1*imnum:(xsize1*(imnum+1))-1,*,85] &$
	ELP_tau1[*,*,i] = mf2.EL_P[xsize1*imnum:(xsize1*(imnum+1))-1,*,73] &$
	ELP_tau2[*,*,i] = mf2.EL_P[xsize1*imnum:(xsize1*(imnum+1))-1,*,61] &$

	GasP_tau0[*,*,i] = mf2.Gas_P[xsize1*imnum:(xsize1*(imnum+1))-1,*,85] &$
	GasP_tau1[*,*,i] = mf2.Gas_P[xsize1*imnum:(xsize1*(imnum+1))-1,*,73] &$
	GasP_tau2[*,*,i] = mf2.Gas_P[xsize1*imnum:(xsize1*(imnum+1))-1,*,61] &$
ENDFOR
SAVE,FILENAME=diro1+'FModel_final_MBP35_series_tau0.sav',B_tau0,T_tau0,V_tau0,ELP_tau0,GasP_tau0
SAVE,FILENAME=diro1+'FModel_final_MBP35_series_taum1.sav',B_tau1,T_tau1,V_tau1,ELP_tau1,GasP_tau1
SAVE,FILENAME=diro1+'FModel_final_MBP35_series_taum2.sav',B_tau2,T_tau2,V_tau2,ELP_tau2,GasP_tau2

