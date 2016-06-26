pro ellambda,compound_name,density,energy,f_cutoff,inel_loss, $
	lambda_el,lambda_elout,lambda_inel,help=help

;+
; pro ellambda,compound_name,density,energy,f_cutoff,inel_loss, $' 
;    lambda_el,lambda_elout,lambda_inel' 
;    (Lambdas in nm returned to you)' 
;  energy in keV can be an array
;  f_cutoff is spatial frequency cutoff in 1/nm of the obj. aperture
;  inel_loss is mean inelastic energy loss in eV
;
; After Langmore 1992
; CJJ Feb. 93
;-

if (keyword_set(help) or (n_elements(compound_name) eq 0)) then begin
  print,'ellambda,compound_name,density,energy,f_cutoff,inel_loss, $'
  print,'       lambda_el,lambda_elout,lambda_inel'
  print,'       (Lambdas in nm returned to you)'
  print,'  energy in keV can be an array'
  print,'  f_cutoff is spatial frequency cutoff in 1/nm of the obj. aperture'
  print,'  inel_loss is mean inelastic energy loss in eV'
  return
endif

compound,compound_name,density,atwt,z_array

el=0.
el_eta=0.
inel=0.

maxz=n_elements(z_array)
for i=0,(maxz-1) do begin
    if (z_array(i) gt 0.) then begin
	this_el=elsigel((i+1),energy)
	this_eleta=this_el*eleta(energy,f_cutoff)
	this_inel=elsiginel((i+1),energy,inel_loss)

	el=el+z_array(i)*this_el
	el_eta=el_eta+z_array(i)*this_eleta
	inel=inel+z_array(i)*this_inel
    endif
endfor

; this is number density times thickness in 1/nm^2
a_over_t=density*602./atwt
lambda_inel=1./(1.e18*inel*a_over_t)
lambda_el=1./(1.e18*el*a_over_t)
lambda_elout=1./(1.e18*el_eta*a_over_t)

return
end

