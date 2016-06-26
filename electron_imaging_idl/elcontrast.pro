pro elcontrast,feature_name,rho_feature,t_feature, $
	background_name,rho_background,thicknesses_background, $
	energy,f_cutoff,inel_loss, $
	cph_filt,cph_nofilt,inel_in_feature, $
	ifnpc_filt,ibnpc_filt,ifnpc_nofilt,ibnpc_nofilt, $
	nophase=nophase,help=help,plot=plot

;+
; PRO elcontrast
;
; elcontrast,feature_name,rho_feature,t_feature, $' 
;     background_name,rho_background,thicknesses_background, $' 
;     energy,f_cutoff,inel_loss, $' 
;     cph_filt,cph_nofilt,inel_in_feature, $' 
;     ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ Returned to you' 
;     ifnpc_filt,ibnpc_filt,ifnpc_nofilt,ibnpc_nofilt' 
;     ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ [Additional]' 
;   rho in g/cm^3, thickness in nm, energy in keV,' 
;     f_cutoff in 1/nm, inel_loss in eV.' 
;   cph_filt: phase contrast, energy filter' 
;   cph_nofilt: phase contrast, no energy filter' 
;   inel_in_feature: fraction of incident electrons which are' 
;     are inelastically scattered by the feature' 
;   ifnpc_filt: intensity from feature column, neg. phase contrast' 
;   /nophase: no phase contrast
;   /plot: show you plots of things on the way' 
;-

if (keyword_set(help) or (n_elements(feature_name) eq 0)) then begin
  print,'elcontrast,feature_name,rho_feature,t_feature, $'
  print,'    background_name,rho_background,thicknesses_background, $'
  print,'    energy,f_cutoff,inel_loss, $'
  print,'    cph_filt,cph_nofilt,inel_in_feature, $'
  print,'    ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ Returned to you'
  print,'    ifnpc_filt,ibnpc_filt,ifnpc_nofilt,ibnpc_nofilt'
  print,'    ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ [Additional]'
  print,'  rho in g/cm^3, thickness in nm, energy in keV,'
  print,'    f_cutoff in 1/nm, inel_loss in eV.'
  print,'  cph_filt: phase contrast, energy filter'
  print,'  cph_nofilt: phase contrast, no energy filter'
  print,'  inel_in_feature: fraction of incident electrons which are'
  print,'    are inelastically scattered by the feature'
  print,'  ifnpc_filt: intensity from feature column, neg. phase contrast'
  print,'  /nophase: no phase contrast
  print,'  /plot: show you plots of things on the way'
  return
endif

; First calculate constants on the feature
compound,feature_name,rho_feature,atwt,z_array

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
a_over_t=rho_feature*602./atwt
kinel_feature=double(1.e18*inel*a_over_t)
kelout_feature=double(1.e18*el_eta*a_over_t)
kelin_feature=double(1.e18*(el-el_eta)*a_over_t)
kappa_feature=kinel_feature+kelout_feature+kelin_feature

; Then calculate constants for the background material
compound,background_name,rho_background,atwt,z_array

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
a_over_t=rho_background*602./atwt
kinel_background=double(1.e18*inel*a_over_t)
kelout_background=double(1.e18*el_eta*a_over_t)
kelin_background=double(1.e18*(el-el_eta)*a_over_t)
kappa_background=kinel_background+kelout_background+kelin_background
tbs=double(thicknesses_background-t_feature)

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
inel_in_feature=exp(-kelout_background*0.5*tbs) * $
	(1.-exp(-kinel_feature*t_feature))

i_in_filt_feature=exp(-(kelout_feature+kinel_feature)*t_feature) * $
	exp(-(kelout_background+kinel_background)*tbs)
i_in_nofilt_feature=exp(-(kelout_feature)*t_feature) * $
	exp(-(kelout_background)*tbs)
i_noscat_feature=exp(-kappa_background*tbs) * $
	exp(-kappa_feature*t_feature)
i_phi_feature=i_noscat_feature*kelin_feature*t_feature

i_in_filt_background=exp(-(kelout_background+kinel_background) * $
	thicknesses_background)
i_in_nofilt_background=exp(-kelout_background*thicknesses_background)
i_noscat_background=exp(-kappa_background*thicknesses_background)
i_phi_background=i_noscat_background*kelin_background*t_feature

if keyword_set(nophase) then begin
  ifnpc_filt=i_in_filt_feature
  ibnpc_filt=i_in_filt_background
endif else begin
  ifnpc_filt=i_in_filt_feature + $
	2.*sqrt(i_phi_feature*i_noscat_feature)
  ibnpc_filt=i_in_filt_background
  negcont=i_in_filt_feature - $
	2.*sqrt(i_phi_feature*i_noscat_feature)
  negindices=where(abs(ifnpc_filt-ibnpc_filt) lt $
	abs(negcont-ibnpc_filt),negcount)
  if (negcount gt 0) then begin
	ifnpc_filt(negindices)=negcont(negindices)
  endif
endelse
cph_filt=(ifnpc_filt-ibnpc_filt)/sqrt(ibnpc_filt)

if keyword_set(nophase) then begin
  ifnpc_nofilt=i_in_nofilt_feature
  ibnpc_nofilt=i_in_nofilt_background
endif else begin
  ifnpc_nofilt=i_in_nofilt_feature + $
	2.*sqrt(i_phi_feature*i_noscat_feature)
  ibnpc_nofilt=i_in_nofilt_background
  negcont=i_in_nofilt_feature - $
	2.*sqrt(i_phi_feature*i_noscat_feature)
  negindices=where(abs(ifnpc_nofilt-ibnpc_nofilt) lt $
	abs(negcont-ibnpc_nofilt),negcount)
  if (negcount gt 0) then begin
	ifnpc_nofilt(negindices)=negcont(negindices)
  endif
endelse
cph_nofilt=(ifnpc_nofilt-ibnpc_nofilt)/sqrt(ibnpc_nofilt)

if keyword_set(plot) then begin
  !x.title='Background thickness in nm'

  !p.title='i_in_filt_feature'
  plot,thicknesses_background,i_in_filt_feature
  wshow,0
  print,'Intensity within aperture with filter, from feature'
  taco=''
  read,'Hit return',taco

  !p.title='i_in_nofilt_feature'
  plot,thicknesses_background,i_in_nofilt_feature
  print,'Intensity within aperture with no filter, from feature'
  taco=''
  read,'Hit return',taco

  !p.title='i_phi_feature'
  plot,thicknesses_background,i_phi_feature
  print,'Single scattering from feature'
  taco=''
  read,'Hit return',taco

  !p.title='inel_in_feature'
  plot,thicknesses_background,inel_in_feature
  print,'Fraction inelastically scattered by feature'
  taco=''
  read,'Hit return',taco

  !p.title='i_in_filt_background'
  plot,thicknesses_background,i_in_filt_background
  print,'Intensity within aperture with filter, from background'
  taco=''
  read,'Hit return',taco

  !p.title='i_in_nofilt_background'
  plot,thicknesses_background,i_in_nofilt_background
  print,'Intensity within aperture with no filter, from background'
  taco=''
  read,'Hit return',taco

  !p.title='i_phi_background'
  plot,thicknesses_background,i_phi_background
  print,'Single scattering from background'
  taco=''
  read,'Hit return',taco

  !p.title='ifnpc_filt, ibnpc_filt, diff'
  a_min=min(ifnpc_filt)
  a_max=max(ifnpc_filt)
  b_min=min(ibnpc_filt)
  b_max=max(ibnpc_filt)
  diff=ifnpc_filt-ibnpc_filt
  d_min=min(diff)
  d_max=max(diff)
  if (a_min gt b_min) then a_min=b_min
  if (a_max lt b_max) then a_max=b_max
  if (a_min gt d_min) then a_min=d_min
  if (a_max lt d_max) then a_max=d_max
  plot,thicknesses_background,ifnpc_filt,yrange=[a_min,a_max]
  oplot,thicknesses_background,ibnpc_filt,linestyle=2
  oplot,thicknesses_background,diff,linestyle=3
  print,'Image intensity with filter, from feature and background'
  taco=''
  read,'Hit return',taco

  !p.title='cph_filt'
  plot,thicknesses_background,cph_filt
  print,'Negative phase contrast with filter'
  taco=''
  read,'Hit return',taco

  !p.title='ifnpc_nofilt, ibnpc_nofilt, diff'
  a_min=min(ifnpc_nofilt)
  a_max=max(ifnpc_nofilt)
  b_min=min(ibnpc_nofilt)
  b_max=max(ibnpc_nofilt)
  diff=ifnpc_nofilt-ibnpc_nofilt
  d_min=min(diff)
  d_max=max(diff)
  if (a_min gt b_min) then a_min=b_min
  if (a_max lt b_max) then a_max=b_max
  if (a_min gt d_min) then a_min=d_min
  if (a_max lt d_max) then a_max=d_max
  plot,thicknesses_background,ifnpc_nofilt,yrange=[a_min,a_max]
  oplot,thicknesses_background,ibnpc_nofilt,linestyle=2
  oplot,thicknesses_background,diff,linestyle=3
  print,'Image intensity with no filter, from feature and background'
  taco=''
  read,'Hit return',taco

  !p.title='cph_nofilt'
  plot,thicknesses_background,cph_nofilt
  print,'Negative phase contrast with no filter'
  taco=''
  read,'Hit return',taco
endif

return
end
