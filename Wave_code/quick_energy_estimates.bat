;Energy Estimates


;F = f rho/2 Vz^2 ct

rho = 4E-5	;kg/m3
ct = 5300.	;m/s
f = 0.2		;generic coeificient
;vz = 3500.	;LOS vel m/s
vz = 674.

e = f*(rho/2.)*(Vz*Vz)*ct

print,e,' W/m2'
print,e/1000.,'k W/m2'

;--------------------------------------
;Lagrangian Displacement

radius = 4.3E6			;(6.3E6)/2.;av radius
P = 105.
percent_a = 0.0106			;Area perturbation
vph = 6260.

Area = 2*!PI*radius^2			;av area assuming circle
Area_n = area + (area*percent_a)	;New area given radial change (perturbations in area signal gives percent_a)
delta_r = SQRT(area_n/(2*!PI)) - radius


;--------------------------------------

;Other method

f = 0.1
rho = 4E-5
omega = (2 * !PI) / P
;vph = 5300.


E = (f*(1 + (ALOG(1/f)))) * (rho/2) * (omega*omega) * (delta_r*delta_r) * vph
print,e/1000.

print,area_n, delta_r, E/1000.
