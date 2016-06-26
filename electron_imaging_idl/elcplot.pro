pro elcplot,nofilt=nofilt,nophase=nophase,log=log, $
	epsfile=epsfile,xinches=xinches,yinches=yinches,$
	titlesize=titlesize,annotsize=annotsize,help=help

if keyword_set(help) then begin
  print,'elcplot,[filename],[/nofilt],[/nophase],[/log]'
  print,'  epsfile=, xinches=, yinches= for EPS file.'
  return
endif

old_plot=!d.name
old_font=!p.font
old_charsize=!p.charsize

if (not keyword_set(xinches)) then xinches=5.5
if (not keyword_set(yinches)) then yinches=4.

if keyword_set(epsfile) then begin
	set_plot,'ps'
	!p.font=0
	device,file=filename,/inch,/portrait, $
	   xsize=xinches,ysize=yinches,/encapsulated
endif

if keyword_set(charsize) then begin
	!p.charsize=charsize
endif

if keyword_set(annotsize) then begin
endif else begin
	annotsize=1.
endelse


thicknesses=20.+780.*findgen(200)/199.

if keyword_set(nofilt) then begin
  !p.title=!p.title+', no energy filter'
endif else begin
  !p.title=!p.title+', with energy filter'
endelse

if keyword_set(nophase) then begin
  elcontrast,'protein',dp,20.,'ice',di,thicknesses,100.,4.12,25., $
	cf,cnf,ief,ifnpc_filt,ibnpc_filt,ifnpc_nofilt,ibnpc_nofilt,/nophase
  !p.title=!p.title+', no phase contrast'
endif else begin
  elcontrast,'protein',dp,20.,'ice',di,thicknesses,100.,4.12,25., $
	cf,cnf,ief,ifnpc_filt,ibnpc_filt,ifnpc_nofilt,ibnpc_nofilt
  !p.title=!p.title+', with phase contrast'
endelse

if keyword_set(nofilt) then begin
	feat=ifnpc_nofilt
	back=ibnpc_nofilt
	diff=ifnpc_nofilt-ibnpc_nofilt
endif else begin
	feat=ifnpc_filt
	back=ibnpc_filt
	diff=ifnpc_filt-ibnpc_filt
endelse

if keyword_set(log) then begin
  plot_io,thicknesses,feat,yrange=[1.e-4,10.],$
    xtitle='Ice thickness in nm', ytitle='Intensity (neg. phase)',$
    title='20 nm protein'


  linethick=[350.,500]
  linepos=[8.e-1,8.e-1]
  dx=550.
  ymult=2.3*2.3
  ymultt=0.9
  oplot,linethick,ymult*linepos
  xyouts,dx,ymult*ymultt*linepos,/data,'Protein present',charsize=annotsize

  ymult=2.3
  oplot,thicknesses,back,linestyle=2
  oplot,linethick,ymult*linepos,linestyle=2
  xyouts,dx,ymult*ymultt*linepos,/data,'Protein absent',charsize=annotsize

  ymult=1.
  oplot,thicknesses,abs(diff),linestyle=3
  oplot,linethick,ymult*linepos,linestyle=3
  xyouts,dx,ymult*ymultt*linepos,/data,'ABS(Difference)',charsize=annotsize
endif else begin
  plot,thicknesses,feat,yrange=[min(-0.4),max(1.2)],$
    xtitle='Ice thickness in nm', ytitle='Intensity (neg. phase)',$
    title='20 nm protein'

  linethick=[300.,450]
  linepos=[0.,0.]
  dx=500.
  dy=1.1
  dyt=-0.02
  oplot,linethick,linepos+dy
  xyouts,dx,dy+dyt,/data,'Protein present',charsize=annotsize

  dy=0.9
  oplot,thicknesses,back,linestyle=2
  oplot,linethick,linepos+dy,linestyle=2
  xyouts,dx,dy+dyt,/data,'Protein absent',charsize=annotsize

  dy=0.7
  oplot,thicknesses,diff,linestyle=3
  oplot,linethick,linepos+dy,linestyle=3
  xyouts,dx,dy+dyt,/data,'Difference',charsize=annotsize
endelse

if keyword_set(epsfile) then begin
	device,/close
	set_plot,old_plot
	!p.font=old_font
	!p.charsize=old_charsize
endif

return
end
