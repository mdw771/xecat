PRO Elxdcomp, epsfile=epsfile, xinches=xinches, yinches=yinches, $
	      xenergy=xenergy, elenergy=elenergy, $
	      feature=feature, rho_feature=rho_feature, $
	      background=background, rho_background=rho_background, $
	      s0=s0, inel_loss=inel_loss, $
	      pixelsize=pixelsize, epsilon=epsilon, $
	      phasecontrast=phasecontrast, nofilter=nofilter, $
	      nophaseabsorption=nophaseabsorption, $
	      thickrange_nm=thickrange_nm, doserange=doserange, $
	      xfrac=xfrac, yfrac=yfrac, $
              xlog=xlog, $
	      titlesize=titlesize, annotsize=annotsize, help=help

IF keyword_set(help) THEN BEGIN
    print, 'elxdcomp,[/nophase]'
    print, ' epsfile=, xinches=, yinches='
    print, ' [feature=,rho_feature= (default protein)]'
    print, ' [background=,rho_background= (default ice)]'
    print, ' [xenergy=eV (otherwise 300 eV)]'
    print, ' [elenergy=keV (otherwise 100 keV)]'
    print, ' [s0= (default 4.12 nm^{-1}]'
    print, ' [inel_loss= (default 40 eV)]'
    print, ' [pixelsize=nm (otherwise 20 nm)]'
    print, ' [epsilon=(fraction of inel sneaking through; 0.01 default)]'
    print, ' [/phasecontrast for EM]'
    print, ' [/nophaseabsorption for perfect Zernike phase plate]'
    print, ' [titlesize=], [annotsize=]'
    return
ENDIF

IF (NOT keyword_set(titlesize)) THEN titlesize = 1.0

IF (NOT keyword_set(annotsize)) THEN annotsize = 1.0

IF (NOT keyword_set(xenergy)) THEN xenergy = 300.

IF (NOT keyword_set(elenergy)) THEN elenergy = 100.

IF (NOT keyword_set(pixelsize)) THEN pixelsize = 20.

IF (NOT keyword_set(epsilon)) THEN epsilon = 0.01

IF (NOT keyword_set(feature)) THEN feature = 'protein'

IF (NOT keyword_set(rho_feature)) THEN rho_feature = 1.

IF (NOT keyword_set(background)) THEN background = 'ice'

IF (NOT keyword_set(rho_background)) THEN rho_background = 1.

IF (NOT keyword_set(inel_loss)) THEN inel_loss = 40.

IF (NOT keyword_set(s0)) THEN s0 = 4.12

IF (NOT keyword_set(thickrange_nm)) THEN thickrange_nm = [20., 1500.]
IF (NOT keyword_set(doserange)) THEN doserange = [1., 1.e10]
IF (NOT keyword_set(xfrac)) THEN xfrac = 0.1

IF (NOT keyword_set(xinches)) THEN xinches = 6.
IF (NOT keyword_set(yinches)) THEN yinches = 4.

if (not keyword_set(xlog)) then xlog=0

IF keyword_set(epsfile) THEN BEGIN
    print, 'Writing PostScript file '+epsfile
    old_plot = !D.name
    old_font = !P.font
    old_charsize = !P.charsize
    set_plot, 'ps'
    !P.font = 0
    device, file = epsfile, /inch, /encap, xsize = xinches, ysize = yinches
ENDIF

feature_name = feature
background_name = background

n_thicknesses = 200
if (xlog eq 0) then begin
   thicknesses = thickrange_nm[0] + $
                 (thickrange_nm[1]-thickrange_nm[0])*$
                 findgen(n_thicknesses)/float(n_thicknesses-1)
endif else begin
   thicknesses = 10.^(alog10(thickrange_nm[0])+$
                      (alog10(thickrange_nm[1])-alog10(thickrange_nm[0]))*$
                      findgen(n_thicknesses)/float(n_thicknesses-1))
endelse

IF (epsilon EQ 1.) THEN epsilon = 0.
print, 'epsilon is ', epsilon

IF keyword_set(phasecontrast) THEN BEGIN
    phasestring = 'Electrons: phase contrast'
    IF (epsilon NE 0.) THEN phasestring = phasestring+','
    eldose, feature, rho_feature, pixelsize, $
      background, rho_background, thicknesses, elenergy, s0, inel_loss, $
      pixelsize, num_filt, num_nofilt, eldose_filt, eldose_nofilt, $
      epsilon = epsilon
    IF keyword_set(nophaseabsorption) THEN BEGIN
	xdose, feature, rho_feature, pixelsize, $
	  background, rho_background, thicknesses, pixelsize, xenergy, $
	  num_amp, num_phase, xdose_amp, xdose_phase, /nophaseabsorption
    ENDIF ELSE BEGIN
	xdose, feature, rho_feature, pixelsize, $
	  background, rho_background, thicknesses, pixelsize, xenergy, $
	  num_amp, num_phase, xdose_amp, xdose_phase
    ENDELSE
ENDIF ELSE BEGIN
    phasestring = 'Electrons: no phase contrast'
    IF (epsilon NE 0.) THEN phasestring = phasestring+','
    eldose, feature, rho_feature, pixelsize, $
      background, rho_background, thicknesses, elenergy, s0, inel_loss, $
      pixelsize, num_filt, num_nofilt, eldose_filt, eldose_nofilt, $
      epsilon = epsilon, /nophase
    IF keyword_set(nophaseabsorption) THEN BEGIN
	xdose, feature, rho_feature, pixelsize, $
	  background, rho_background, thicknesses, pixelsize, xenergy, $
	  num_amp, num_phase, xdose_amp, xdose_phase, /nophaseabsorption
    ENDIF ELSE BEGIN
	xdose, feature, rho_feature, pixelsize, $
	  background, rho_background, thicknesses, pixelsize, xenergy, $
	  num_amp, num_phase, xdose_amp, xdose_phase
    ENDELSE
ENDELSE

print,'xray num_amp min,max: ',min(num_amp),max(num_amp)
print,'xray num_phase min,max: ',min(num_phase),max(num_phase)
print,'xdose_amp min,max: ',min(xdose_amp),max(xdose_amp)
print,'xdose_phase min,max: ',min(xdose_phase),max(xdose_phase)
print,'doserange: ',doserange

elfiltmin = min(eldose_filt, max = elfiltmax)
elnofiltmin = min(eldose_nofilt, max = elnofiltmax)
xampmin = min(xdose_amp, max = xampmax)
xphasemin = min(xdose_phase, max = xphasemax)
ymin = min([elfiltmin, elnofiltmin, xampmin, xphasemin])
ymax = max([elfiltmax, elnofiltmax, xampmax, xphasemax])

xrange=[thicknesses[0],thicknesses[n_thicknesses-1]]
xstyle=xlog

plot, thicknesses, eldose_nofilt,$
      xlog=xlog,xrange=xrange,xstyle=xstyle,$
      /ylog,yrange=doserange,ystyle=1,$
      charsize = titlesize, linestyle = 1, $
      xtitle = 'Ice thickness in nm', $
      ytitle = 'Dose in Gray', $
      title = strtrim(string(pixelsize, format = '(i10)'), 2) + $
      ' nm '+feature_name+' in '+background_name+', '+ $
      strtrim(string(elenergy, format = '(i10)'), 2) + $
      ' keV electrons, '+ $
      strtrim(string(xenergy, format = '(i10)'), 2)+' eV x-rays'

IF (NOT keyword_set(nofilter)) THEN BEGIN
    oplot, thicknesses, eldose_filt
ENDIF
oplot, thicknesses, xdose_amp, linestyle = 2
oplot, thicknesses, xdose_phase, linestyle = 3

IF (NOT keyword_set(yfrac)) THEN BEGIN
    IF (ymin GT 1.e4) THEN BEGIN
	yfrac = 0.18
    ENDIF ELSE BEGIN
	yfrac = 0.78
    ENDELSE
ENDIF

xline = thickrange_nm(0) + $
  (thickrange_nm(1)-thickrange_nm(0))*[xfrac, xfrac+0.08]
yline = [1., 1.]*doserange(0)*10.^(yfrac*alog10(doserange(1)/doserange(0)))
xtit = thickrange_nm(0) + $
  (thickrange_nm(1)-thickrange_nm(0))*(xfrac+0.1)
ymult = 1./(0.07*alog10(doserange(1)/doserange(0)))
print, 'ymult = ', ymult
ymultt = 0.9

IF (NOT keyword_set(nofilter)) THEN BEGIN
    oplot, xline, yline*ymult^2.
ENDIF
oplot, xline, yline*ymult^3., linestyle = 1
oplot, xline, yline*ymult^4., linestyle = 2
oplot, xline, yline*ymult^5., linestyle = 3
IF (NOT keyword_set(nofilter)) THEN BEGIN
    xyouts, xtit, yline(0)*ymult^2.*ymultt, /data, 'Electrons, filter'
ENDIF
xyouts, xtit, yline(0)*ymult^3.*ymultt, /data, 'Electrons, no filter'
xyouts, xtit, yline(0)*ymult^4.*ymultt, /data, 'X-rays, amplitude'
xyouts, xtit, yline(0)*ymult^5.*ymultt, /data, 'X-rays, Zernike phase'
xyouts, xline(0), yline(0)*ymult^1.*ymultt, /data, phasestring
print, 'again, epsilon = ', epsilon
IF (epsilon NE 0.) THEN BEGIN
    xyouts, xline(0), yline(0)*ymultt, /data, $
      '  !Me!X='+strtrim(string(epsilon, format = '(f5.3)'), 2)
ENDIF

IF (n_elements(filename) GT 0) THEN BEGIN
    device, /close
    set_plot, old_plot
    !P.font = old_font
    !P.charsize = old_charsize
ENDIF

return
END

