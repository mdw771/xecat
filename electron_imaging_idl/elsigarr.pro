pro elsigarr,zarr,enarr,s0,inel_loss,el,el_eta,inel,help=help

;+
; pro elsigarr,zarr,enarr,s0,inel_loss,el,el_eta,inel
;
;	Given zarr,enarr, makes arrays of el(z,en),
;		eleta(z,en,s0), inel(z,en,inel_loss)
;
;-

if (keyword_set(help) or (n_elements(zarr) eq 0)) then begin
	print,'elsigarr,zarr,enarr,s0,inel_loss,el,el_eta,inel'
	return
endif

numz=n_elements(zarr)
numen=n_elements(enarr)
el=fltarr(numz,numen)
el_eta=fltarr(numz,numen)
inel=fltarr(numz,numen)

for eni=0,(numen-1) do begin
	en=enarr(eni)
	el(0:(numz-1),eni)=elsigel(zarr,en)
	el_eta(0:(numz-1),eni)=el(0:(numz-1),eni)*eleta(en,s0)
	inel(0:(numz-1),eni)=elsiginel(zarr,en,inel_loss)
endfor

return
end
