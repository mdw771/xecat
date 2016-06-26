pro elthetaplot,epsfile=epsfile,xinches=xinches,yinches=yinches,$
	energy=energy,epsilon=epsilon, $
	maxthick=maxthick, $
	eps=eps,encapsulated=encapsulated, $
	titlesize=titlesize,annotsize=annotsize,help=help

if keyword_set(help) then begin
	print,'elthetaplot,[filename],[/eps] or [/encapsulated]'
	print,'  [energy=(keV)],[maxthick=(nm)],[epsilon=] '
	print,'  /titlesize=, /annotsize='
	return
endif

if (not keyword_set(xinches)) then xinches=5.5
if (not keyword_set(yinches)) then yinches=4.

if keyword_set(epsfile) then beign
  print,'Writing Encapsulated PostScript file '+epsfile
  old_plot=!d.name
  old_font=!p.font
  old_charsize=!p.charsize
  set_plot,'ps'
  !p.font=0
  device,file=epsfile,/encap,/inch,xsize=xinches,ysize=yinches
endif

if (not keyword_set(titlesize)) then titlesize=1.

if (not keyword_set(annotsize)) then annotsize=1.0

if (not keyword_set(epsilon)) then epsilon=0.

if (not keyword_set(energy)) then energy=100.

if (not keyword_set(maxthick)) then maxthick=2000.
x
f_cutoff=4.12
t_f=20.
inel_loss=25.

thicknesses=20.+(maxthick-20.)*findgen(200)/199.
eltheta,'protein',dp,t_f,'ice',di,thicknesses, $
	energy,f_cutoff,inel_loss, $
	theta_b,theta_bphi,theta_bf,theta_bfphi,inel_in_f, $
	epsilon=epsilon

plot_io,thicknesses,theta_bfphi,yrange=[1.e-6,1.],charsize=titlesize,$
  xtitle='Ice thickness in nm', $
  ytitle='Contrast !MQ!X', $
  title=strtrim(string(energy,format='(i5)'),2) + $
    ' keV, s!S!D0!R ='+strtrim(string(f_cutoff,format='(f4.2)'),2) + $
    '!S!U-1!R, '+strtrim(string(t_f,format='(i5)'),2)+' nm protein feature',$
  subtitle='!Me!X='+strtrim(string(epsilon,format='(f8.4)'),2)

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

if keyword_set(epsfile) then begin
	device,/close
	set_plot,old_plot
	!p.font=old_font
	!p.charsize=old_charsize
endif

return
end

