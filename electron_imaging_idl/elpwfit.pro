pro elpwfit,nstart,nend

if (keyword_set(help) or (n_elements(nstart) eq 0)) then begin
	print,'elpwfit,nstart,nend'
	return
end

if (!version.os eq 'windows') then begin
	backslash='\'
	header=backslash+'idl'+backslash+'lib'+backslash+'local'+ $
		backslash+'x1'+backslash
endif else if (!version.os eq 'AIX') then begin
	header='~/idl/'
endif

read_mapper,header+'langpwall.map',zarr,pwall
read_mapper,header+'langpw2.map',zarr,pw2
read_mapper,header+'langpw5.map',zarr,pw5

lnz=alog(zarr)

rat2=(pw2/pwall)/0.8
rat5=(pw5/pwall)/0.5

ratall=0.5*rat2+0.5*rat5

plot,lnz,ratall,yrange=[0.,2.],psym=1

allz=findgen(99)+1.
for n=nstart,nend do begin
	fit=poly_fit(lnz,ratall,n)
	print,'Terms for polynomial order ',n
	print,'  0: ',fit(0,0)
	pfit=fit(0,0)
	for pn=1,n do begin
		pfit=pfit+fit(0,pn)*(lnz^float(pn))
		print,'  '+strtrim(string(pn,format='(i5)'),2)+': ',fit(0,pn)
	endfor
	oplot,lnz,pfit,linestyle=(n-nstart+1)
	print,'Again: ',fit
endfor

return
end
