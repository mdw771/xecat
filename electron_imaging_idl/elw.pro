pro elw,energy,cs,defocus,f_scherzer,freqs,w,normdefocus=normdefocus, $
	help=help

if (keyword_set(help) or (n_elements(energy) eq 0)) then begin
  print,'elctf,energy,cs,defocus,f_scherzer,freqs,w'
  print,'                        ^^^^^^^^^^^^^^^^^^ Returned to you'
  print,'  energy in keV, cs (spherical aberration) in meters,'
  print,'    defocus as fraction of Scherzer, frequencies in 1/meters'
  print,'  [,/normdefocus: give defocus as fraction of Scherzer, '
  print,'     return in meters]
  return
endif

electron_mass=511.003414 ; keV/c^2
c=2.99792458e8 ; m/sec
h=6.62617636e-34 ; eV sec
joule_ev=1./1.602189246e-19
gamma=double(1.+energy/electron_mass)
beta=sqrt(1.-1./(gamma*gamma))
lambda=h*joule_ev*c/(gamma*1000.*electron_mass*beta)

f_scherzer=sqrt(2.)/(cs*lambda*lambda*lambda)^0.25

if (n_elements(freqs) lt 2) then begin
	freqs=findgen(100)
endif
w=freqs
freqs=(freqs/max(freqs))*2.*f_scherzer

w=0.5*!pi*(cs*lambda^3.0*freqs^4.0 - $
	2*(defocus*sqrt(cs*lambda))*lambda*freqs*freqs)

if keyword_set(normdefocus) then begin
	defocus=defocus*sqrt(cs*lambda)
endif

return
end



