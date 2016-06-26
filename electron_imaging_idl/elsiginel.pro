function elsiginel,z,v,energy_loss,help=help

;+
; function elsiginel,z,v,energy_loss
;
; given Z of an atom, electron voltage V in keV, and
; average inelastic scatter energy loss energy_loss in eV, calculates
; the total inelastic scattering cross section sig_inel in
; m^2 according to formula 5 in
;	J.P. Langmore and M.F. Smith, ``Quantitatitve energy-filtered
;	electron microscopy of biological molecules in ice'',
;	{\it Ultramicroscopy \bf 46}, 349--373 (1992).
;
; CJJ Feb. 93
;-

if (keyword_set(help) or (n_elements(z) eq 0)) then begin
	print,'result=elsiginelg(z,v,energy_loss)'
	print,'  returns sig_inel in m^2 with V in keV'
	print,'  and energy_loss in eV'
	return,0.
endif

electron_mass=511.003414 ; keV/c^2
c=2.99792458e8 ; m/sec
h=6.62617636e-34 ; eV sec
joule_ev=1./1.602189246e-19
gamma=double(1.+v/electron_mass)
beta=sqrt(1.-1./(gamma*gamma))
lambda=h*joule_ev*c/(gamma*1000.*electron_mass*beta)

theta_e=1.0e-3*energy_loss/(beta*beta*(v+electron_mass))

if (n_elements(z) eq 1) then begin
  local_z=1.
  local_z=z(0)
  if (abs(local_z-1.) lt 0.1) then begin
	gamma100=double(1.+100./electron_mass)
	beta100sq=1.-1./(gamma100*gamma100)
	const=8.8/(1.5*alog(2./theta_e))*beta100sq
	local_z=const*const
  endif
endif else begin
  local_z=z
  hydrogen_indices=where((abs(local_z-1.) lt 0.1),h_count)
  if (h_count gt 0) then begin
	gamma100=double(1.+100./electron_mass)
	beta100sq=1.-1./(gamma100*gamma100)
	const=8.8/(1.5*alog(2./theta_e))*beta100sq
	local_z(hydrogen_indices)=const*const
  endif
endelse

; the 1.0e-18 is to return in m^2 rather than nm^2
result=1.0e-18*1.5e-6*sqrt(local_z)*alog(2./theta_e)/(beta*beta)

return,result
end
