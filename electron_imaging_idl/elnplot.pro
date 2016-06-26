pro elnplot,epsfile=epsfile,xinches=xinches,yinches=yinches,$
	,nophase=nophase,nofilter=nofilter, $
	epsilon=epsilon, $
	eps=eps,encapsulated=encapsulated, $
	titlesize=titlesize,annotsize=annotsize,help=help

if keyword_set(help) then begin
	print,'elnplot,[filename],[/nophase],[/nofilter] for PostScript file'
	print,'  /eps or /encapsulated: Encapsulated PostScript'
	print,' [epsilon=(fraction of inel sneaking through)]'
	return
endif

if keyword_set(epsfile) then begin
	print,'Writing EPS file '+epsfile
	old_plot=!d.name
	old_font=!p.font
	old_charsize=!p.charsize
	set_plot,'ps'
	!p.font=0
	device,file=epsfile,/inch,/encap,xsize=xinches,ysize=yinches
endif

if (not keyword_set(epsilon)) then epsilon=0.0

if (not keyword_set(titlesize)) then titlesize=1.

if (not keyword_set(annotsize)) then annotsize=1.


thicknesses=20.+1480.*findgen(200)/199.

if keyword_set(nophase) then begin
    eldose,'protein',dp,20.,'ice',di,thicknesses,100.,4.12,25.,20., $
      num_filt,num_nofilt,dose_filt,dose_nofilt,/nophase,epsilon=epsilon
endif else begin
    eldose,'protein',dp,20.,'ice',di,thicknesses,100.,4.12,25.,20., $
      num_filt,num_nofilt,dose_filt,dose_nofilt,epsilon=epsilon
endelse

if keyword_set(nofilter) then begin
  num=num_nofilt
  if keyword_set(nophase) then begin
    !p.title=!p.title+', no phase contrast'
  endif else begin
    !p.title=!p.title+', with phase contrast'
  endelse
endif else begin
  num=num_filt
  if keyword_set(nophase) then begin
    !p.title=!p.title+', no phase contrast'
  endif else begin
    !p.title=!p.title+', with phase contrast'
  endelse
endelse

plot_io,thicknesses,num,yrange=[1.e2,1.e10],charsize=titlesize,$
  xtitle='Ice thickness in nm', $
  ytitle='Required electrons/(20 nm)!S!U2!R', $
  title='20 nm protein, 100 keV'


xline=[100.,300]
yline=[1.e8,1.e8]
xtit=350.
ymult=3.
ymultt=0.85

filterstring='Filter (!Me!X=' + $
  strtrim(string(epsilon,format='(f8.4)'),2)+')'
oplot,xline,yline
xyouts,xtit,yline*ymultt,/data,filterstring,charsize=annotsize

oplot,thicknesses,num_nofilt,linestyle=2
oplot,xline,yline*ymult,linestyle=2
xyouts,xtit,yline*ymult*ymultt,/data, $
	'Without filter',charsize=annotsize

if keyword_set(epsfile) then begin
	device,/close
	set_plot,old_plot
	!p.font=old_font
	!p.charsize=old_charsize
endif

return
end
