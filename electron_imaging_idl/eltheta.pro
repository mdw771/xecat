pro eltheta,feature_name,rho_f,t_f, $
	background_name,rho_b,thicknesses, $
	energy,f_cutoff,inel_loss, $
	theta_b,theta_bphi,theta_bf,theta_bfphi,inel_in_f, $
	epsilon=epsilon,help=help

;+
; PRO eltheta
;
; eltheta,feature_name,rho_f,t_f, $' 
;     background_name,rho_b,thicknesses, $' 
;     energy,f_cutoff,inel_loss,epsilon, $' 
;     theta_b,theta_bphi,theta_bf,theta_bfphi,inel_in_f
;
;   rho in g/cm^3, thickness in nm, energy in keV,' 
;     f_cutoff in 1/nm, inel_loss in eV.' 
;   theta: contrast as follows...
;   _b: bright field, no phase contrast, no energy filter' 
;   _bphi: bright field, phase contrrast, no energy filter' 
;   _bf: bright field, no phase contrast, energy filter' 
;   _bfphi: bright field, phase contrast, energy filter' 
;   inel_in_f: fraction of incident electrons which are' 
;    are inelastically scattered by the feature' 
;   /epsilon=(leakage of inelastic electrons in filter)
;-

if (keyword_set(help) or (n_elements(feature_name) eq 0)) then begin
  print,'eltheta,feature_name,rho_f,t_f, $'
  print,'    background_name,rho_b,thicknesses, $'
  print,'    energy,f_cutoff,inel_loss, $'
  print,'    theta_b,theta_bphi,theta_bf,theta_bfphi,inel_in_f'
  print,''
  print,'  rho in g/cm^3, thickness in nm, energy in keV,'
  print,'    f_cutoff in 1/nm, inel_loss in eV.'
  print,'  theta: contrast as follows...
  print,'  _b: bright field, no phase contrast, no energy filter'
  print,'  _bphi: bright field, phase contrrast, no energy filter'
  print,'  _bf: bright field, no phase contrast, energy filter'
  print,'  _bfphi: bright field, phase contrast, energy filter'
  print,'  inel_in_f: fraction of incident electrons which are'
  print,'    are inelastically scattered by the feature'
  return
endif

if keyword_set(epsilon) then begin
	local_epsilon=epsilon
endif else begin
	local_epsilon=0.0
endelse

; First calculate constants on the feature
compound,feature_name,rho_f,atwt,z_array

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
a_over_t=rho_f*602./atwt
kinel_f=double(1.e18*inel*a_over_t)
kelout_f=double(1.e18*el_eta*a_over_t)
kelin_f=double(1.e18*(el-el_eta)*a_over_t)
ktot_f=kinel_f+kelout_f+kelin_f

; Then calculate constants for the background material
compound,background_name,rho_b,atwt,z_array

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
a_over_t=rho_b*602./atwt
kinel_b=double(1.e18*inel*a_over_t)
kelout_b=double(1.e18*el_eta*a_over_t)
kelin_b=double(1.e18*(el-el_eta)*a_over_t)
ktot_b=kinel_b+kelout_b+kelin_b
tbs=double(thicknesses-t_f)

; If the contrast is used to calculate exposure and dose, we
; will find out how many electrons are needed on the whole
; thickness of background plus feature in order to see the
; feature.  However, we also want to find out what fraction
; of those electrons produce inelastic scattering in the
; feature.  We'll assume the feature is halfway within the background,
; so that the number of electrons getting to the feature is
; given by the number of electrons remaining at that point.
; Then find out the number which are inelastically scattered
; by the feature.
inel_in_f=exp(-kelout_b*0.5*tbs) * $
	(1.-exp(-kinel_f*t_f))

i_noscat_b=exp(-ktot_b*thicknesses)
i_oneel_b=kelin_b*thicknesses*i_noscat_b
i_innoinel_b=exp(-(kinel_b+kelout_b)*thicknesses)
i_multel_b=i_innoinel_b-i_noscat_b-i_oneel_b
i_in_b=exp(-kelout_b*thicknesses)
i_elout_b=1.-i_in_b
i_inel_b=i_in_b-i_innoinel_b

i_noscat_f=exp(-ktot_b*tbs)*exp(-ktot_f*t_f)
i_oneel_f=(kelin_b*tbs+kelin_f*t_f)*i_noscat_f
i_oneelf_f=kelin_f*t_f*i_noscat_f
i_innoinel_f=exp(-(kinel_b+kelout_b)*tbs) * $
	exp(-(kinel_f+kelout_f)*t_f)
i_multel_f=i_innoinel_f-i_noscat_f-i_oneel_f
i_in_f=exp(-kelout_b*tbs)*exp(-kelout_f*t_f)
i_elout_f=1.-i_in_f
i_inel_f=i_in_f-i_innoinel_f

; no phase contrast, no filter
theta_b=abs(i_innoinel_b-i_innoinel_f)/sqrt(i_in_f+i_in_b)
; no phase contrast, with filter
theta_bf=abs(i_innoinel_b-i_innoinel_f) / $
	sqrt(i_innoinel_f+i_innoinel_b+local_epsilon*i_inel_f)
; phase contrast, no filter
theta_bphi=(abs(i_innoinel_b-i_innoinel_f) + $
		2.*sqrt(i_noscat_f*i_oneelf_f)) / $
	sqrt(i_in_f+i_in_b)
; phase contrast, with filter
theta_bfphi=(abs(i_innoinel_b-i_innoinel_f) + $
		2.*sqrt(i_noscat_f*i_oneelf_f)) / $
	sqrt(i_innoinel_f+i_innoinel_b+local_epsilon*i_inel_f)

return
end

