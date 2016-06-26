pro ellambda_test

elements,names,z,weights,densities

kev = 100.
inel_loss = 40.

el_array = fltarr(92)
inel_array = fltarr(92)
z = indgen(92)+1

for i=0,91 do begin
   el_array[i] = elsigel(z[i],kev)
   inel_array[i] = elsiginel(z[i],kev,inel_loss)
endfor

old_multi = !p.multi
!p.multi = [0,1,2,0,0]

plot,z,el_array,xtitle='Z',$
     ytitle='Elastic mean free path at '+$
     strtrim(string(kev,format='(i3)'),2)+' keV'
plot,z,inel_array,xtitle='Z',$
     ytitle='Inelastic mean free path at '+$
     strtrim(string(kev,format='(i3)'),2)+' keV'

end
