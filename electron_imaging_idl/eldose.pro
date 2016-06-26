pro eldose,feature_name,rho_feature,t_feature, $
	background_name,rho_background,thicknesses_background, $
	energy,f_cutoff,inel_loss,pixel_size, $
	num_filt,num_nofilt,dose_filt,dose_nofilt, $
	nophase=nophase,epsilon=epsilon,help=help

if (keyword_set(help) or (n_elements(feature_name) eq 0)) then begin
  print,'eldose,feature_name,rho_feature,t_feature, $'
  print,'   background_name,rho_background,thicknesses_background, $'
  print,'   energy,f_cutoff,inel_loss,pixel_size, $'
  print,'   num_filt,num_nofilt,dose_filt,dose_nofilt'
  print,'   ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ Returned to you'
  print,' rho in g/cm^3, thickness in nm, energy in keV, f_cutoff'
  print,'   in 1/nm, inel_loss in eV, pixel_size in nm, dose in Gy.'
  print,' /nophase: no phase contrast'
  print,' epsilon= (0.01 default; fraction of inel sneaking through)'
  return
endif

if (not keyword_set(epsilon)) then epsilon=0.01

eltheta,feature_name,rho_feature,t_feature, $
	background_name,rho_background,thicknesses_background, $
	energy,f_cutoff,inel_loss, $
	theta_b,theta_bphi,theta_bf,theta_bfphi,ief, $
	epsilon=epsilon

; Feature mass in kg: 1.e-3 for kg from g, 1.e-21 for nm^3 -> cm^3
inverse_feature_mass=1./(1.e-24*rho_feature*t_feature*pixel_size*pixel_size)

if keyword_set(nophase) then begin
  num_filt=25./(theta_bf*theta_bf)
  fixers=where((num_filt lt 25.),fixcount)
  if (fixcount gt 0) then num_filt(fixers)=25.
  num_nofilt=25./(theta_b*theta_b)
  fixers=where((num_nofilt lt 25.),fixcount)
  if (fixcount gt 0) then num_nofilt(fixers)=25.
endif else begin
  num_filt=25./(theta_bfphi*theta_bfphi)
  fixers=where((num_filt lt 25.),fixcount)
  if (fixcount gt 0) then num_filt(fixers)=25.
  num_nofilt=25./(theta_bphi*theta_bphi)
  fixers=where((num_nofilt lt 25.),fixcount)
  if (fixcount gt 0) then num_nofilt(fixers)=25.
endelse

dose_filt=1.6e-19*inel_loss*num_filt*ief*inverse_feature_mass
dose_nofilt=1.6e-19*inel_loss*num_nofilt*ief*inverse_feature_mass

return
end

