pro eldarr,feature_name,rho_feature,t_feature, $
	background_name,rho_background,thicknesses_background, $
	energies,f_cutoff,inel_loss,pixel_size, $
	nums_filt,nums_nofilt,doses_filt,doses_nofilt, $
	nophase=nophase,epsilon=epsilon,help=help

if (keyword_set(help) or (n_elements(feature_name) eq 0)) then begin
  print,'eldarr,feature_name,rho_feature,t_feature, $'
  print,'   background_name,rho_background,thicknesses_background, $'
  print,'   energies,f_cutoff,inel_loss,pixel_size, $'
  print,'   nums_filt,nums_nofilt,doses_filt,doses_nofilt'
  print,'   ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ Returned to you'
  print,' rho in g/cm^3, thickness in nm, energies in keV, f_cutoff'
  print,'   in 1/nm, inel_loss in eV, pixel_size in nm, dose in Gy.'
  print,' /nophase: no phase contrast'
  print,' epsilon=(fraction of inel sneaking through)'
  return
endif

if keyword_set(epsilon) then begin
endif else begin
	epsilon=0.0
endelse

num_energies=n_elements(energies)
num_thicknesses=n_elements(thicknesses_background)

nums_filt=fltarr(num_energies,num_thicknesses)
nums_nofilt=fltarr(num_energies,num_thicknesses)
doses_filt=fltarr(num_energies,num_thicknesses)
doses_nofilt=fltarr(num_energies,num_thicknesses)

for eni=0,(num_energies-1) do begin
  energy=energies(eni)

  if keyword_set(nophase) then begin
    eldose,feature_name,rho_feature,t_feature, $
	background_name,rho_background,thicknesses_background, $
	energy,f_cutoff,inel_loss,pixel_size, $
	filt_num,nofilt_num,filt_dose,nofilt_dose, $
	epsilon=epsilon,/nophase
  endif else begin
    eldose,feature_name,rho_feature,t_feature, $
	background_name,rho_background,thicknesses_background, $
	energy,f_cutoff,inel_loss,pixel_size, $
	filt_num,nofilt_num,filt_dose,nofilt_dose, $
	epsilon=epsilon
  endelse
  nums_filt(eni,0:(num_thicknesses-1))=filt_num(0:(num_thicknesses-1))
  nums_nofilt(eni,0:(num_thicknesses-1))=nofilt_num(0:(num_thicknesses-1))
  doses_filt(eni,0:(num_thicknesses-1))=filt_dose(0:(num_thicknesses-1))
  doses_nofilt(eni,0:(num_thicknesses-1))=nofilt_dose(0:(num_thicknesses-1))
endfor

return
end
