pro elthetaplot,filename,energy=energy,epsilon=epsilon, $
	maxthick=maxthick, $
	eps=eps,encapsulated=encapsulated, $
	titlesize=titlesize,annotsize=annotsize,help=help

if keyword_set(help) then begin
	print,'elthetaplot,[filename],[/eps] or [/encapsulated]'
	print,'  [energy=(keV)],[maxthick=(nm)],[epsilon=] '
	print,'  /titlesize=, /annotsize='
	return
endif

if (n_elements(filename) gt 0) then begin
  print,'Writing PostScript file ',strtrim(filename,2)
  old_name=!d.name
  old_font=!p.font
  set_plot,'ps'
  !p.font=0
  if (keyword_set(eps) or keyword_set(encapsulated)) then begin
	device,file=filename,/portrait,/encapsulated
  endif else begin
	device,file=filename,/inch,/landscape, $
		xsize=7.,ysize=5.
  endelse
endif

if keyword_set(titlesize) then begin
endif else begin
	titlesize=1.
endelse

if keyword_set(annotsize) then begin
endif else begin
	annotsize=1.0
endelse

if keyword_set(epsilon) then begin
endif else begin
	epsilon=0.
endelse

if keyword_set(energy) then begin
endif else begin
	energy=100.
endelse

if keyword_set(maxthick) then begin
endif else begin
	maxthick=2000.
endelse

f_cutoff=4.12
t_f=20.
inel_loss=25.
!x.title='Ice thickness in nm'
!y.title='Contrast !MQ!X'
!p.title=strtrim(string(energy,format='(i5)'),2) + $
	' keV, s!S!D0!R ='+strtrim(string(f_cutoff,format='(f4.2)'),2) + $
	'!S!U-1!R,  '+strtrim(string(t_f,format='(i5)'),2)+' nm protein feature'
!p.subtitle='!Me!X='+strtrim(string(epsilon,format='(f8.4)'),2)

thicknesses=20.+(maxthick-20.)*findgen(200)/199.
eltheta,'protein',dp,t_f,'ice',di,thicknesses, $
	energy,f_cutoff,inel_loss, $
	theta_b,theta_bphi,theta_bf,theta_bfphi,inel_in_f, $
	epsilon=epsilon

plot_io,thicknesses,theta_bfphi,yrange=[1.e-6,1.],charsize=titlesize
oplot,thicknesses,theta_bf,linestyle=2
oplot,thicknesses,theta_bphi,linestyle=3
oplot,thicknesses,theta_b,linestyle=1

lthick=[100.,500.]
ltheta=[6.e-6,6.e-6]
xtext=550.
ymultt=0.85
ymult=2.5

oplot,lthick,ltheta*(ymult^3.)
xyouts,xtext,ltheta*(ymult^3.)*ymultt,'BF!Mj!X',charsize=annotsize

oplot,lthick,ltheta*(ymult^2.),linestyle=2
xyouts,xtext,ltheta*(ymult^2.)*ymultt,'BF',charsize=annotsize

oplot,lthick,ltheta*ymult,linestyle=3
xyouts,xtext,ltheta*ymult*ymultt,'B!Mj!X',charsize=annotsize

oplot,lthick,ltheta,linestyle=1
xyouts,xtext,ltheta*ymultt,'B',charsize=annotsize

if (n_elements(filename) gt 0) then begin
	device,/close
	set_plot,old_name
	!p.font=old_font
endif

return
end

