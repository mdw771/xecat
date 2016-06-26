PRO ellamplot,compound_name,density,filename,max=max, $
              s0=s0,inel_loss=inel_loss, $
              eps=eps,encapsulated=encapsulated, $
              titlesize=titlesize,annotsize=annotsize,help=help
   
IF (keyword_set(help) OR (n_elements(compound_name) EQ 0)) $
  THEN BEGIN
    print,'ellamplot,compound_name,density,[filename],[max=val]'
    print,'  s0= (4.12 nm^{-1} default), inel_loss= (35 ev default)'
    print,'  /eps or /encapsulated: Encapsulated PostScript'
    return
ENDIF

IF (n_elements(filename) GT 0) THEN BEGIN
    print,'Writing PostScript file ',filename
    old_name=!d.name
    old_font=!p.font
    set_plot,'ps'
    !p.font=0
    IF (keyword_set(eps) OR keyword_set(encapsulated)) THEN BEGIN
        device,file=filename,/inch,/portrait, $
          xsize=6.,ysize=5.,/encapsulated
    ENDIF ELSE BEGIN
        device,file=filename,/inch,/landscape, $
          xsize=6.,ysize=5.
    ENDELSE
ENDIF

IF keyword_set(titlesize) THEN BEGIN
    !p.charsize=titlesize
ENDIF

IF keyword_set(annotsize) THEN BEGIN
ENDIF ELSE BEGIN
    annotsize=1.0
ENDELSE

IF (NOT keyword_set(s0)) THEN s0=4.12
IF (NOT keyword_set(inel_loss)) THEN inel_loss=35.

!x.title='Energy in keV'
!y.title='!Ml!X in nm'
!p.title=compound_name+', s!S!D0!R ='+strtrim(string(s0,format='(f8.2)'),2)+ $
  ' nm!S!U-1!R,  '+strtrim(string(inel_loss,format='(i3)'),2)+ $
  ' eV inel. loss'
!p.subtitle=''

energies=80.+1480.*findgen(100)/99.

ellambda,compound_name,density,energies,s0,inel_loss, $
  lambda_el,lambda_elout,lambda_inel,help=help

kel=1./lambda_el
kelout=1./lambda_elout
kelin=kel-kelout
lambda_elin=1./kelin

allmin=min([min(lambda_el),min(lambda_inel)])
allmax=max([max(lambda_el),max(lambda_inel)])

IF keyword_set(max) THEN BEGIN
    allmin=0.
    allmax=max
ENDIF

help, energies
help, lambda_inel
help, allmin
help, allmax

plot_oi,energies,lambda_inel,yrange=[allmin,allmax], $
  xrange=[100,1500],xstyle=1

xyouts, 200., !y.crange(0)-0.05*(!y.crange(1)-!y.crange(0)), /data, $
  '200', align=0.5
xyouts, 400., !y.crange(0)-0.05*(!y.crange(1)-!y.crange(0)), /data, $
  '400', align=0.5
xyouts, 1500., !y.crange(0)-0.05*(!y.crange(1)-!y.crange(0)), /data, $
  '1500', align=0.5


oplot,energies,lambda_el,linestyle=1
;oplot,energies,lambda_elout,linestyle=2
;oplot,energies,lambda_elin,linestyle=3

xvec=[120.,170]
xt=190.
ymult=1.
ymultt=1.
yvec=[1.,1.]*(!y.crange(0)+0.7*(!y.crange(1)-!y.crange(0)))
plots,xvec,ymult*yvec
yvec=[1.,1.]*(!y.crange(0)+0.69*(!y.crange(1)-!y.crange(0)))
xyouts,xt,ymult*yvec(0)*ymultt,/data,'inelastic',charsize=annotsize

ymult=1.
yvec=[1.,1.]*(!y.crange(0)+0.78*(!y.crange(1)-!y.crange(0)))
plots,xvec,ymult*yvec,linestyle=1
yvec=[1.,1.]*(!y.crange(0)+0.77*(!y.crange(1)-!y.crange(0)))
xyouts,xt,ymult*yvec(0)*ymultt,/data,'elastic',charsize=annotsize

;ymult=1.
;yvec=[1.,1.]*(!y.crange(0)+0.86*(!y.crange(1)-!y.crange(0)))
;plots,xvec,ymult*yvec,linestyle=2
;yvec=[1.,1.]*(!y.crange(0)+0.85*(!y.crange(1)-!y.crange(0)))
;xyouts,xt,ymult*yvec(0)*ymultt,/data,'elastic, out',charsize=annotsize

;ymult=1.
;yvec=[1.,1.]*(!y.crange(0)+0.92*(!y.crange(1)-!y.crange(0)))
;plots,xvec,ymult*yvec,linestyle=3
;yvec=[1.,1.]*(!y.crange(0)+0.91*(!y.crange(1)-!y.crange(0)))
;xyouts,xt,ymult*yvec(0)*ymultt,/data,'elastic, in',charsize=annotsize

IF (n_elements(filename) GT 0) THEN BEGIN
    device,/close
    set_plot,old_name
    !p.font=old_font
ENDIF

return
END




