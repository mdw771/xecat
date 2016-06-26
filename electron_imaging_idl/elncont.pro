pro elncont,epsfile=epsfile,xinches=xinches,yinches=yinches,$
	feature=feature,rho_feature=rho_feature, $
	background=background,rho_background=rho_background, $
	pixel_size=pixel_size,nofilt=nofilt,nophase=nophase, $
	s0=s0,inel_loss=inel_loss,epsilon=epsilon, $
	epsilon=epsilon, $
	eps=eps,encapsulated=encapsulated, $
	titlesize=titlesize,annotsize=annotsize,help=help

if (keyword_set(help) or (n_elements(s0) eq 0)) then begin
	print,'elncont,[/nofilt],[/nophase]'
	print,'  epsfile=, xinches=, yinches='
	print,'  feature=,rho_feature= (protein default)'
	print,'  background=,rho_background= (ice default)'
	print,'  pixel_size= (20 nm default)'
	print,'  s0= (4.12 nm^{-1} default)'
	print,'  inel_loss= (25 eV default)
	print,'  /nofilt: no energy filter'
	print,'  /nophase: no phase contrast'
	print,'  /eps or /encapsulated: Encapsulated PostScript'
	print,'  epsilon= (0.01 default; fraction of inel sneaking through)'
	return
endif

if (not keyword_set(xinches)) then xinches=5.5
if (not keyword_set(yinches)) then yinches=5.5

if keyword_set(epsfile) then begin
	print,'Making EPS file '+epsfile
	old_plot=!d.name
	old_font=!p.font
	old_charsize=!p.charsize
	!p.font=0
	set_plot,'ps'
	device,file=epsfile,/inch,/encap,xsize=xinches,ysize=yinches
endif

if (not keyword_set(titlesize)) then titlesize=1.0

if (not keyword_set(annotsize)) then annotsize=1.0

if (not keyword_set(epsilon)) then epsilon=0.01

if (not keyword_set(feature)) then feature='protein'

if (not keyword_set(rho_feature)) then rho_feature=1.

if (not keyword_set(background)) then background='ice'

if (not keyword_set(rho_background)) then rho_background=1.

if (not keyword_set(pixel_size)) then pixel_size=20.

if (not keyword_set(s0)) then s0=4.12

if (not keyword_set(inel_loss)) then inel_loss=25.


thicknesses=20.+1980.*findgen(100)/99.
energies=100.+900.*findgen(100)/99.

if keyword_set(nophase) then begin
  eldarr,feature,rho_feature,pixel_size, $
	background,rho_background,thicknesses, $
	energies,s0,inel_loss,pixel_size, $
	nums_filt,nums_nofilt,doses_filt,doses_nofilt, $
	epsilon=epsilon,/nophase
endif else begin
  eldarr,feature,rho_feature,pixel_size, $
	background,rho_background,thicknesses, $
	energies,s0,inel_loss,pixel_size, $
	nums_filt,nums_nofilt,doses_filt,doses_nofilt, $
	epsilon=epsilon
endelse

my_subtitle=''

if keyword_set(nofilt) then begin
	to_plot=nums_nofilt/400.
	my_subtitle='No energy filter'
endif else begin
	to_plot=nums_filt/400.
	my_subtitle='Energy filter (!Me!X=' + $
		strtrim(string(epsilon,format='(f8.4)'),2)+')'
endelse

if keyword_set(nophase) then begin
	my_subtitle=my_subtitle+', no phase contrast'
endif else begin
	my_subtitle=my_subtitle+', with phase contrast'
endelse

mylevels=[1.,3.16,1.e1,3.16e1,1.e2,3.16e2,1.e3,3.16e3,1.e4,3.16e4, $
	1.e5,3.16e5,1.e6,3.16e6,1.e7,3.16e7,1.e8,3.16e8,1.e9,3.16e9,1.e10]
myannot=['10!U0','','10!U1','','10!U2','','10!U3','','10!U4','', $
	'10!U5','','10!U6','','10!U7','','10!U8','','10!U9','','10!U10']
mylabels=[1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1]
mylinestyles=[0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0]
contour,to_plot,energies,thicknesses,c_charsize=annotsize, $
	levels=mylevels,c_annotation=myannot,c_labels=mylabels, $
	c_linestyle=mylinestyles,xrange=[100,1000],xtype=1, $
	xticks=5,xtickv=[100.,200.,400.,1000.],$
	xtickname=['100','200','400','1000'],$
	xtitle='Electron energy in keV', $
	ytitle='Ice thickness in nm', $
	title='Electrons/nm!S!U2!R  for ' + $
	  strtrim(string(pixel_size,format='(i10)'),2)+ $
	  ' nm '+feature+', s!S!D0!R ='+ $
	  strtrim(string(s0,format='(f3.1)'),2)+' nm'

if keyword_set(epsfile) then begin
	device,/close
	set_plot,old_plot
	!p.font=old_font
	!p.charsize=old_charsize
endif

return
end

