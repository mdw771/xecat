pro elcatplot,energy,filename,normal=normal, $
              compound_name=compound_name,density=density, $
              s0=s0,inel_loss=inel_loss, $
              eps=eps,encapsulated=encapsulated, $
              max_thickness_nm=max_thickness_nm,$
              titlesize=titlesize,annotsize=annotsize,debug=debug,help=help

if (keyword_set(help) or (n_elements(energy) eq 0)) then begin
  print,'elcatplot,energy,[filename]'
  print,'  compound_name= (ice default), density= (1.0 default)'
  print,'  s0= (4.12 nm^{-1} default), inel_loss= (25 eV default)'
  print,'  /eps or /encapsulated: Encapsulated PostScript'
  print,'  /titlesize, /annotsize, /debug'
  return
endif

if (n_elements(filename) gt 0) then begin
	print,'Writing PostScript file ',filename
	old_name=!d.name
	old_font=!p.font
	set_plot,'ps'
	!p.font=0
	if (keyword_set(eps) or keyword_set(encapsulated)) then begin
		device,file=filename,/inch,/portrait, $
			xsize=6.,ysize=5.,/encapsulated
	endif else begin
		device,file=filename,/inch,/landscape, $
			xsize=6.,ysize=5.
	endelse
endif

if (not keyword_set(titlesize)) then titlesize=1.0

if (not keyword_set(annotsize)) then annotsize=1.0

if (not keyword_set(compound_name)) then compound_name='ice'

if (not keyword_set(density)) then density=1.0

temp_compound_name=strupcase(strmid(compound_name,0,1)) + $
	strmid(compound_name,1,strlen(compound_name)-1)
!x.title=temp_compound_name+' thickness in nm'
!p.title=strtrim(string(energy,format='(i5)'),2) + $
	' keV, s!S!D0!R =4.1 nm!S!U-1!R'
!p.subtitle=''

svec = size(max_thickness_nm)
if (svec[svec[0]+1] eq 0) then begin
   if (energy le 200.) then begin
      max_thicness_nm = 1000.
   endif else begin
      max_thickness_nm = 3000.
   endelse
endif
thicknesses=max_thickness_nm*findgen(1000)/999.

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

normalization=i_noscat+i_one_el+i_elout+i_inel+i_mult_el

if keyword_set(normal) then begin
  !y.title='Normalization check'
  plot,thicknesses,normalization, $
	yrange=[min(normalization),max(normalization)]
  print,'Normalization from ',min(normalization),' to ',max(normalization)
endif else begin
  !y.title='Fraction of intensity'
  plot_io,thicknesses,i_noscat,yrange=[0.001,1.]
  oplot,thicknesses,i_one_el,linestyle=1
  oplot,thicknesses,i_mult_el,linestyle=2
  oplot,thicknesses,i_inel,linestyle=3
  oplot,thicknesses,i_elout,linestyle=4
  tlabel=['noscat','1el','multel','inel','elout']

  thisy=0.02
  xvec=[550.,750.]
  xtext=780.
  if (energy le 200.) then begin
  endif else begin
    xvec=xvec*3
    xtext=xtext*3
  endelse

  yvec=[0.,0.]
  ymult=1.35
  ymultt=0.92

  for i=0,4 do begin
    oplot,xvec,thisy+yvec,linestyle=i
    xyouts,xtext,thisy*ymultt,/data,tlabel(i),charsize=annotsize
    thisy=thisy*ymult
  endfor
endelse

if (n_elements(filename) gt 0) then begin
	device,/close
	set_plot,old_name
	!p.font=old_font
endif

return
end

