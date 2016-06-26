pro elcatplot_si,energy=energy,epsfile=epsfile,normal=normal, $
                 compound_name=compound_name,density=density, $
                 s0=s0,inel_loss=inel_loss, $
                 max_thickness_nm=max_thickness_nm,$
                 titlesize=titlesize,annotsize=annotsize,debug=debug,help=help
  
  if keyword_set(help) then begin
     print,'elcatplot,energy,[epsfilename]'
     print,'  compound_name= (ice default), density= (1.0 default)'
     print,'  s0= (4.12 nm^{-1} default), inel_loss= (25 eV default)'
     print,'  /titlesize, /annotsize, /debug'
     return
  endif
  
  if keyword_set(epsfile) then begin
     epsfilename = 'elcatplot_si_300kv_raw.eps'
     print,'Writing PostScript file "'+epsfilename+'"'
     old_name=!d.name
     old_font=!p.font
     set_plot,'ps'
     device,file=epsfilename,/inch, $
            xsize=6.,ysize=5.,/encapsulated
     !p.font=0
  endif
  
  if (not keyword_set(energy)) then energy=300.
  if (not keyword_set(titlesize)) then titlesize=1.0
  if (not keyword_set(annotsize)) then annotsize=1.0
  if (not keyword_set(compound_name)) then compound_name='Si'
  if (not keyword_set(density)) then density=2.32
  
  x_title = compound_name+' thickness in nm'
  plot_title = strtrim(string(energy,format='(i5)'),2) + $
             ' kV'
  
  svec = size(max_thickness_nm)
  if (svec[svec[0]+1] eq 0) then begin
     if (energy le 200.) then begin
        max_thickness_nm = 1000.
     endif else begin
        max_thickness_nm = 2000.
     endelse
  endif
  thicknesses = max_thickness_nm*findgen(200)/199.
  
  if (not keyword_set(s0)) then s0=4.12
  if (not keyword_set(inel_loss)) then inel_loss=25.
  if keyword_set(debug) then begin
     elcategory,compound_name,density,thicknesses, $
                energy,s0,inel_loss, $
                i_noscat,i_one_el,i_mult_el,i_elout,i_inel,/debug
  endif else begin
     elcategory,compound_name,density,thicknesses, $
                energy,s0,inel_loss, $
                i_noscat,i_one_el,i_mult_el,i_elout,i_inel
  endelse
  
  normalization = i_noscat+i_one_el+i_elout+i_inel+i_mult_el
  
  if keyword_set(normal) then begin
     !y.title='Normalization check'
     plot,thicknesses,normalization, $
          yrange=[min(normalization),max(normalization)]
     print,'Normalization from ',min(normalization),' to ',max(normalization)
  endif else begin
     y_title = 'Fraction of intensity'
     plot_io,thicknesses,i_noscat,yrange=[0.001,1.],/nodata,$
             xtitle=x_title,ytitle=y_title,title=plot_title
     oplot,thicknesses,i_noscat,color=250
     oplot,thicknesses,i_one_el,linestyle=1,color=250
     oplot,thicknesses,i_mult_el,linestyle=2,color=250
     oplot,thicknesses,i_inel,linestyle=3,color=250
     oplot,thicknesses,i_elout,linestyle=4,color=250
     tlabel=['noscat','1el','multel','inel','elout']
     
     thisy=0.02
     xvec=[1100.,1400.]
     xtext=1500.
     
     yvec=[0.,0.]
     ymult=1.35
     ymultt=0.92
     
     for i=0,4 do begin
        oplot,xvec,thisy+yvec,linestyle=i
        xyouts,xtext,thisy*ymultt,/data,tlabel(i),charsize=annotsize
        thisy=thisy*ymult
     endfor
  endelse
  
  if keyword_set(epsfile) then begin
     device,/close
     set_plot,old_name
     !p.font=old_font
  endif
  
  return
end

