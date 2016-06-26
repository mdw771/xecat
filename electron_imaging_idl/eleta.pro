function eleta,v,s0,alpha=alpha, $
; polynomial=polynomial,
	help=help

;+
; function eleta,v,s0
;
; given voltage v in keV and either:
;	spatial frequency cutoff of objective aperture s0 in
;		inverse nm (typically 2-5 nm^{-1})
;	if /alpha switch is used, objective aperture half-angle
;		in radians (typically 0.008-0.018)
;;	/polynomial
; calculates the fraction of elastically scattered electrons
; that are scattered outside the objective aperture according to
; formula 4 in
;	J.P. Langmore and M.F. Smith, ``Quantitatitve energy-filtered
;	electron microscopy of biological molecules in ice'',
;	{\it Ultramicroscopy \bf 46}, 349--373 (1992).
;
; CJJ Feb. 93
;-

if (keyword_set(help) or (n_elements(v) eq 0)) then begin
	print,'result=eleta(v,s0,[/alpha])'
	print,'  v is in keV, returns sigma_el in m^2 with'
	print,'  objective aperture half-angle s0 in nm^-1 or alpha in radians'
	return,0.
endif

electron_mass=511.003414 ; keV/c^2
c=2.99792458e8 ; m/sec
h=6.62617636e-34 ; eV sec
joule_ev=1./1.602189246e-19
gamma=double(1.+v/electron_mass)
beta=sqrt(1.-1./(gamma*gamma))
lambda=h*joule_ev*c/(gamma*1000.*electron_mass*beta)

; Langmore eq 4 uses s0 in nm^{-1}
if keyword_set(alpha) then begin
	s0=1.0e-9*2.*sin(alpha0*0.5)/lambda
endif
eta=1-0.1*s0

; if keyword_set(polynomial) then begin
;	polyfit=[0.93460987,-126.66198,251.65296,-202.59160, $
;		84.715841,-19.433256,2.3220602,-0.11307297]
;	pn=n_elements(polyfit)
;
;	lnz=alog(z)
;
;	correction=polyfit(0)
;	for pni=1,(pn-1) do begin
;		correction=correction+polyfit(pni)*(lnz^float(pni))
;	endfor
;	eta=eta*correction
; endif

return,eta
end
