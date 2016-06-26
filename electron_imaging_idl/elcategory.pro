pro elcategory,compound_name,density,thicknesses, $
	energy,f_cutoff,inel_loss, $
	i_noscat,i_one_el,i_mult_el,i_elout,i_inel, $
	debug=debug

if (keyword_set(help) or (n_elements(energy) eq 0)) then begin
  print,'elcategory,compound,density,thickness,energy,f_cutoff,inel_loss,'
  print,'  i_noscat,i_one_el,i_multel,i_elout,i_inel'
  return
endif

compound,compound_name,density,atwt,z_array

el=0.
el_eta=0.
inel=0.

maxz=n_elements(z_array)
for i=0,(maxz-1) do begin
    if (z_array(i) gt 0.) then begin
	this_z=i+1
	this_el=elsigel(this_z,energy)
	this_eta=eleta(energy,f_cutoff)
	this_inel=elsiginel(this_z,energy,inel_loss)

	el=el+z_array(i)*this_el
	el_eta=el_eta+z_array(i)*this_el*this_eta
	inel=inel+z_array(i)*this_inel
    endif
endfor

; this is number density times thickness in 1/nm^2
a_over_t=density*602./atwt
kinel=double(1.e18*inel*a_over_t)
kelout=double(1.e18*el_eta*a_over_t)
kelin=double(1.e18*(el-el_eta)*a_over_t)
kel=kelin+kelout
kappa=kinel+kelout+kelin

if keyword_set(debug) then begin
  print,'Compound: '+compound_name+'  Density: ' + $
	strtrim(string(density,format='(f6.3)'),2) + $
	' g/cm^3.  Atomic weight: ' + $
	strtrim(string(atwt,format='(f10.1)'),2)+' g/mol.  Compositon:'
  zshow,z_array
  print,'          pm^2/atwt   lambda (nm)'
  print,'el        ' + $
	strtrim(string((1.0e24*el/atwt),format='(f9.2)'),2) + $
	'    '+strtrim(string((1./kel),format='(f7.2)'),2)
  print,'el,in     ' + $
	strtrim(string((1.0e24*(el-el_eta)/atwt),format='(f9.2)'),2) + $
	'    '+strtrim(string((1./kelin),format='(f7.2)'),2)
  print,'el,out    ' + $
	strtrim(string((1.0e24*el_eta/atwt),format='(f9.2)'),2) + $
	'    '+strtrim(string((1./kelout),format='(f7.2)'),2)
  print,'inel      ' + $
	strtrim(string((1.0e24*inel/atwt),format='(f9.2)'),2) + $
	'    '+strtrim(string((1./kinel),format='(f7.2)'),2)
  print,'tot       ' + $
	strtrim(string((1.0e24*(el+inel)/atwt),format='(f9.2)'),2) + $
	'    '+strtrim(string((1./kappa),format='(f7.2)'),2)
endif

i_noscat=exp(-kappa*thicknesses)
i_one_el=i_noscat*kelin*thicknesses
i_in=exp(-kelout*thicknesses)
i_elout=1.-i_in
i_in_noinel=exp(-(kelout+kinel)*thicknesses)
i_inel=i_in-i_in_noinel
i_mult_el=i_in_noinel - i_noscat - i_one_el

return
end



