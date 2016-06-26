PRO compound,compound_string,density,atwt,z_array,help=help

;+
; pro compound,compound_string,density,z_array
;
;	Type in a chemical formula and the density.  The routine
;	then figures out the atomic weight in g/mole of "molecules"
;	and how many of what Z atoms you have.	The
;	density isn't really needed here; it's just that if
;	what you type is defined in
;		(your directory) compound.dat
;               (a system directory; see below) compound.dat
;	then you can type
;		compound,'protein',density,atwt,z_array
;	and get both the density and z_array set for you.
;
;	Note that you can do "zshow,z_array" to check out the
;	composition.
;
; CJJ Feb 93
;
; BH 2005-07-12: add 'on_error, 2'
;-

  on_error, 2

  IF (keyword_set(help) OR (n_elements(compound_string) EQ 0)) THEN BEGIN
      print,'compound,compound_string,density,atwt,z_array'
      print,'  Looks in "compound.dat" first. The parameter'
      print,'  "atwt" is the atomic/molecular weight in g/mole.'
      z_array = 0
      return
  ENDIF
  
  tabchar=string(9B)
  backslash='\'
  z_array = 0

  osname = strupcase(strmid(!version.os,0,3))

  ;; changed BH 2004-06-28: use IDL function 'file_which' and get rid
  ;; of findpro etc. 
  ;; changed BH 2009-03-28: Rather than searching for compound.pro 
  ;; directly (which also exists in runtime_maps directory), search 
  ;; for compound.pro and then for compound.dat in the same directory
  compound_pro_path = file_which('compound.pro')
;  system_compound_file = file_which('compound.dat')
  system_compound_file = file_dirname(compound_pro_path, /mark_directory) + $
                         'compound.dat'


  ;; changed BH 2004-03-09 to make independent of absolute paths 
  ;; get path for this procedure -- compound.dat will be in same directory
;   findpro, 'compound.pro', /noprint, dirlist=dirlist
;   this_dir = dirlist[0]
;   IF (osname EQ 'WIN') THEN BEGIN
;       system_compound_file = this_dir+'\compound.dat'
;   ENDIF ELSE BEGIN
;       system_compound_file = this_dir+'/compound.dat'
;   ENDELSE
;   IF (osname EQ 'WIN') THEN BEGIN 
;       system_compound_file = !DIR+'\lib\idl_local\henke\compound.dat'
;   ENDIF ELSE BEGIN
;       system_compound_file=!DIR+'/lib/idl_local/henke/compound.dat'
;   ENDELSE
  
  on_ioerror,try_system_file
  openr,lun,'compound.dat',/get_lun
  goto,got_file
  
  try_system_file: on_ioerror,no_file
  openr,lun,system_compound_file,/get_lun
  goto,got_file
  
  no_file: print,'Could not open file "compound.dat" or "', $
    system_compound_file,'".'
  close,lun
  free_lun,lun
  z_array = 0
  return
  
  got_file: on_ioerror,nomatch
  input_string=''
  WHILE (NOT eof(lun)) DO BEGIN
      readf,lun,input_string
      input_string=strtrim(input_string,1)
      last_char_index=strlen(input_string)-1
      space_index=-1
      char_index=0
      test_char=strmid(input_string,char_index,1)
      WHILE ((char_index LT last_char_index) and (test_char NE ' ') AND $
             (test_char NE tabchar)) DO BEGIN
          char_index=char_index+1
          test_char=strmid(input_string,char_index,1)
      ENDWHILE
      
      ;; we should have found a space or a tabchar
      IF ((test_char NE ' ') AND (test_char NE tabchar)) THEN BEGIN
          message,'no white space to second word in compound.dat'
          z_array = 0
          return
      ENDIF ELSE BEGIN 
          space_index=char_index
      ENDELSE
      
      IF ((strmid(input_string,0,1) NE '!') AND $
          (space_index NE -1)) THEN BEGIN
          firstword=strmid(input_string,0,space_index)
          ;; check compound name case-insensitive, BH 2005-02-15
;          IF (firstword EQ compound_string) THEN BEGIN
          IF strcmp(firstword, compound_string, /fold) THEN BEGIN 
              last_char_index=strlen(input_string)-1
              secondword=strtrim(strmid(input_string,space_index, $
                                        (last_char_index-space_index+1)),1)
              
              last_char_index=strlen(secondword)-1
              space_index=-1
              char_index=0
              test_char=strmid(secondword,char_index,1)
              WHILE ((char_index LT last_char_index) AND $
                     (test_char NE ' ') AND (test_char NE tabchar)) DO BEGIN
                  char_index=char_index+1
                  test_char=strmid(secondword,char_index,1)
              ENDWHILE
              
              ;; we should have found a space or a tabchar
              IF ((test_char NE ' ') AND (test_char NE tabchar)) THEN BEGIN
                  message,'no white space to third word in compound.dat'
                  z_array = 0
                  return
              ENDIF ELSE BEGIN
                  space_index=char_index
              ENDELSE
              
              thirdword=strmid(secondword,space_index, $
                               (last_char_index-space_index+1))
              secondword=strtrim(strmid(secondword,0,space_index),2)
              thirdword=strtrim(thirdword,2)
              density=0.
              reads,thirdword,density
              compound_string=secondword
              zcompound,compound_string,z_array
              atwt=0.
              zatwt,atwt,z_array
              close,lun
              free_lun,lun
              return
          ENDIF
      ENDIF
  ENDWHILE 
  
  nomatch: close,lun
  free_lun,lun
  
  zcompound,compound_string,z_array
  atwt=0.
  zatwt,atwt,z_array
  
  return
  
confused:
  message,'got confused'
  close,lun
  free_lun,lun
  z_array = 0
  return
  
END

