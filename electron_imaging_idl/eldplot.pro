pro eldplot,epsfile=epsfile,xinches=xinches,yinches=yinches,$
	feature=feature,rho_feature=rho_feature, $
	background=background,rho_background=rho_background, $
	pixel_size=pixel_size, max_thickness=max_thickness, $
	s0=s0,energy=energy, $
	nophase=nophase,epsilon=epsilon, $
	eps=eps,encapsulated=encapsulated, $
	titlesize=titlesize,annotsize=annotsize,help=help

if keyword_set(help) then begin
	print,'eldplot,[/nophase]'
	print,'  feature=,rho_feature= (protein default)'
	print,'  background=,rho_background= (ice default)'
	print,'  pixel_size= (20 nm default)'
	print,'  max_thickness= (1500 nm default)'
	print,'  s0= (4.12 nm^{-1} default)'
	print,'  energy= (100 keV default)'
	print,'  epsfile=,xinches=,yinches='
	print,'  epsilon= (0.01 default; fraction of inel sneaking through)'
	return
endif

if (not keyword_set(epsilon)) then epsilon=0.01

if (not keyword_set(feature)) then feature='protein'

if (not keyword_set(rho_feature)) then rho_feature=1.

if (not keyword_set(background)) then background='ice'

if (not keyword_set(rho_background)) then rho_background=1.

if (not keyword_set(pixel_size)) then pixel_size=20.

if (not keyword_set(energy)) then energy=100.

if (not keyword_set(s0)) then s0=4.12

if (not keyword_set(max_thickness)) then max_thickness=1500.

if (not keyword_set(xinches)) then xinches=5.5
if (not keyword_set(yinches)) then yinches=4.

if keyword_set(epsfile) then begin
	print,'Writing PostScript file ',epsfile
	old_plot=!d.name
	old_font=!p.font
	old_charsize=!p.charsize
	set_plot,'ps'
	!p.font=0
	device,file=epsfile,/inch,/eps,xsize=xinches,ysize=yinches
endif

if (not keyword_set(titlesize)) then titlesize=1.0

if (not keyword_set(annotsize)) then annotsize=1.0


thicknesses=pixel_size+(max_thickness-pixel_size)*findgen(200)/199.

if keyword_set(nophase) then begin
  eldose,feature,rho_feature,pixel_size, $
	background,rho_background,thicknesses,energy,s0,25.,pixel_size, $
	num_filt,num_nofilt,dose_filt,dose_nofilt, $
	epsilon=epsilon,/nophase
  !p.title=!p.title+', no phase contrast'
endif else begin
  eldose,feature,rho_feature,pixel_size,$
	background,rho_background,thicknesses,energy,s0,25.,pixel_size, $
	num_filt,num_nofilt,dose_filt,dose_nofilt, $
	epsilon=epsilon
  !p.title=!p.title+', with phase contrast'
endelse

plot_io,thicknesses,dose_filt,yrange=[1.e1,1.e8],charsize=titlesize,$
  xtitle=background+' thickness in nm',$
  ytitle='Dose in Gray', $
  title=strtrim(string(pixel_size,format='(i10)'),2)+ $
    ' nm '+feature+', '+strtrim(string(energy,format='(i10)'),2)+' keV'


xline=[700.,950.]
yline=[7.e1,7.e1]
xtit=1000.
ymult=2.5
ymultt=0.85

filterstring='Filter (!Me!X='+strtrim(string(epsilon, $
	format='(f8.4)'),2)+')'
oplot,xline,yline
xyouts,xtit,yline*ymultt,/data,filterstring,charsize=annotsize

oplot,thicknesses,dose_nofilt,linestyle=2
oplot,xline,yline*ymult,linestyle=2
xyouts,xtit,yline*ymult*ymultt,/data, $
	'No filter',charsize=annotsize

if (n_elements(filename) gt 0) then begin
	device,/close
	set_plot,old_name
	!p.font=old_font
	!p.charsize=old_charsize
endif

return
end
