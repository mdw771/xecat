function elsigel,z,v,help=help

;+
; function elsigel(z,v)
;
; given Z of an atom and electron voltage V in keV, calculates
; the total elastic scattering cross section sig_el in
; m^2 according to formula 1 in
;	J.P. Langmore and M.F. Smith, ``Quantitatitve energy-filtered
;	electron microscopy of biological molecules in ice'',
;	{\it Ultramicroscopy \bf 46}, 349--373 (1992).
;
; CJJ Feb. 93
;-

if (keyword_set(help) or (n_elements(z) eq 0)) then begin
	print,'result=elsigel(z,v)'
	print,'  returns sig_el in m^2 with V in keV'
	return,0.
endif

electron_mass=511.003414
c=2.99792458e8
gamma=double(1.+v/electron_mass)
beta=sqrt(1.-1./(gamma*gamma))

z=float(z)
; the 1.0e-18 is to return in m^2 rather than nm^2
result=1.0e-18*1.4e-6*(z^1.5)*(1.-(0.26*z/(137.*beta)))/(beta*beta)

return,result
end
