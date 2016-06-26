pro elwarr,energy,cs,defocus,f_scherzer,array,pixel_size, $
	normdefocus=normdefocus,help=help

if (keyword_set(help) or (n_elements(energy) eq 0)) then begin
  print,'elwarr,energy,cs,defocus,f_scherzer,array,pixel_size'
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

f_max=1./(2.*pixel_size)
f_max_sq=f_max*f_max

print,'  Maximum frequencies: Scherzer=' + $
	strtrim(string(f_scherzer*1.e-6,format='(f10.1)'),2)+ $
	', your array='+strtrim(string(f_max*1.e-6,format='(f10.1)'),2)+$
	' inverse microns.'

nx=n_elements(array(*,0))
ny=n_elements(array(0,*))

xarr=(findgen(nx)-0.5*nx)/(0.5*nx)
xarr=xarr*xarr

; Make the array be the frequency squared
for yindex=0,(ny-1) do begin
	array(0:(nx-1),yindex)= $
		xarr(0:(nx-1))+xarr(yindex)
endfor
array=sqrt(array)*f_max

; this used to have lambda^3*array^4 in the first term, and
; lambda*array^2 in the second term
array=array*lambda
array=0.5*!pi*(cs*array^4.0/lambda - $
	2*(defocus*sqrt(cs*lambda))*array*array/lambda)

arrmax=max(array,min=arrmin)
if (abs(arrmin) gt abs(arrmax)) then arrmax=arrmin
print,'  Maximum aberration is '+strtrim(string(arrmax,format='(f10.2)'),2)+$
	' radians.'

if keyword_set(normdefocus) then begin
	defocus=defocus*sqrt(cs*lambda)
endif

return
end



