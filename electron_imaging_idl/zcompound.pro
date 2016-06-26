;+
; pro zcompound,compound_string,z_array,paren_multiplier=paren_multiplier
;
; compound_string: like 'Si3N4'
;
; if called without paren_multiplier set to anything, it makes
;	zarr=fltarr(92)
;
; BH 2005-07-12: add 'on_error, 2'
;-

PRO zcompound,compound_string,z_array,paren_multiplier=paren_multiplier, $
              help=help

  on_error, 2
  
  if (keyword_set(help) or (n_elements(compound_string) eq 0) ) then begin
      print,'zcompound,compound_string,z_array'
      return
  endif
  
  if (strlen(compound_string) eq 0) then return
  
  last_char_index = strlen(compound_string) - 1
  
  if keyword_set(paren_multiplier) then begin
  endif else begin
      z_array=fltarr(92)
      paren_multiplier=1.
  endelse
  
  max_z_index=n_elements(z_array)+1
  
  ;; If we don't start off with a parenthesis, all we have to do
  ;; is strip off the first element and process it.  We then
  ;; call the routine over again to handle the next part of
  ;; the string...
  if (strmid(compound_string,0,1) ne '(') then begin
      ;; Look to see if the string has an element
      ;; like "C" or like "He".
      first_char=strmid(compound_string,0,1)
      second_char=strmid(compound_string,1,1)
      this_element_name = first_char
      if ((second_char ge 'a') and (second_char le 'z')) then begin
          this_element_name = this_element_name + second_char
          num_start_index = 2
      endif else begin
          this_element_name = this_element_name + ' '
          num_start_index = 1
      endelse
      ;; print,'this_element_name:',this_element_name,', num_start_index:', $
      ;;   num_start_index,', last_char_index:',last_char_index
      
      this_z=0
      case this_element_name of
          'H ': this_z=1
          'He': this_z=2
          'Li': this_z=3
          'Be': this_z=4
          'B ': this_z=5
          'C ': this_z=6
          'N ': this_z=7
          'O ': this_z=8
          'F ': this_z=9
          'Ne': this_z=10
          'Na': this_z=11
          'Mg': this_z=12
          'Al': this_z=13
          'Si': this_z=14
          'P ': this_z=15
          'S ': this_z=16
          'Cl': this_z=17
          'Ar': this_z=18
          'K ': this_z=19
          'Ca': this_z=20
          'Sc': this_z=21
          'Ti': this_z=22
          'V ': this_z=23
          'Cr': this_z=24
          'Mn': this_z=25
          'Fe': this_z=26
          'Co': this_z=27
          'Ni': this_z=28
          'Cu': this_z=29
          'Zn': this_z=30
          'Ga': this_z=31
          'Ge': this_z=32
          'As': this_z=33
          'Se': this_z=34
          'Br': this_z=35
          'Kr': this_z=36
          'Rb': this_z=37
          'Sr': this_z=38
          'Y ': this_z=39
          'Zr': this_z=40
          'Nb': this_z=41
          'Mo': this_z=42
          'Tc': this_z=43
          'Ru': this_z=44
          'Rh': this_z=45
          'Pd': this_z=46
          'Ag': this_z=47
          'Cd': this_z=48
          'In': this_z=49
          'Sn': this_z=50
          'Sb': this_z=51
          'Te': this_z=52
          'I ': this_z=53
          'Xe': this_z=54
          'Cs': this_z=55
          'Ba': this_z=56
          'La': this_z=57
          'Ce': this_z=58
          'Pr': this_z=59
          'Nd': this_z=60
          'Pm': this_z=61
          'Sm': this_z=62
          'Eu': this_z=63
          'Gd': this_z=64
          'Tb': this_z=65
          'Dy': this_z=66
          'Ho': this_z=67
          'Er': this_z=68
          'Tm': this_z=69
          'Yb': this_z=70
          'Lu': this_z=71
          'Hf': this_z=72
          'Ta': this_z=73
          'W ': this_z=74
          'Re': this_z=75
          'Os': this_z=76
          'Ir': this_z=77
          'Pt': this_z=78
          'Au': this_z=79
          'Hg': this_z=80
          'Tl': this_z=81
          'Pb': this_z=82
          'Bi': this_z=83
          'Po': this_z=84
          'At': this_z=85
          'Rn': this_z=86
          'Fr': this_z=87
          'Ra': this_z=88
          'Ac': this_z=89
          'Th': this_z=90
          'Pa': this_z=91
          'U ': this_z=92
          else: this_z=0
      endcase
      
      ;; print,'this_z=',this_z
      if (this_z eq 0) then begin
          print,'zcompound is confused: ',compound_string
          compound_string=''
          z_array = 0
          return
      endif
      
      ;; Find the next element or parenthesis, as
      ;; anything before it must be a number.
      postnum_index = num_start_index
      test_char=strmid(compound_string,postnum_index,1)
      while ( ((test_char eq '0') or (test_char eq '1') or $
               (test_char eq '2') or (test_char eq '2') or $
               (test_char eq '3') or (test_char eq '4') or $
               (test_char eq '5') or (test_char eq '6') or $
               (test_char eq '7') or (test_char eq '8') or $
               (test_char eq '9') or (test_char eq '.')) and $
              (postnum_index le last_char_index) ) do begin
          postnum_index=postnum_index+1
          if (postnum_index le last_char_index) then begin
              test_char=strmid(compound_string,postnum_index,1)
          endif else begin
              test_char=''
          endelse
      endwhile
      
      ;; print,'postnum_index:',postnum_index,', test_char:',test_char
      
      ;; is there more?
      if (num_start_index ne postnum_index) then begin
          number_string=strmid(compound_string,num_start_index, $
                               (postnum_index-num_start_index))
          num_multiplier=1.
          ;; print,'Trying to interpret ',number_string,' as a number.'
          if (strlen(number_string) ne 0) then begin
              reads,number_string,num_multiplier
          endif
      endif else begin
          num_multiplier=1.
      endelse
      
      ;; We've handled this element, so pop it into the
      ;; matrix and continue.
      if (this_z le max_z_index) then begin
          z_array(this_z-1) = z_array(this_z-1) + $
            num_multiplier * paren_multiplier
      endif else begin
          print,'zcompound: z_array smaller than ',max_z_index
          z_array = 0
          return
      endelse
      
      ;; And deal with what's left
      remaining_string=strmid(compound_string,postnum_index, $
                              (last_char_index-postnum_index+1))
      zcompound,remaining_string,z_array,paren_multiplier=paren_multiplier
  endif else begin
      ;; Now deal with the case where we start off with a left
      ;; parenthesis. Begin by finding the matching right
      ;; parenthesis.
      nest_level = 1            ;
      paren_end_index = 0       ;
      char_index = 1            ;
      test_char=strmid(compound_string,char_index,1)
      while ((nest_level gt 0) and (test_char ne '')) do begin
          if (test_char eq '(') then nest_level=nest_level+1
          if (test_char eq ')') then nest_level=nest_level-1
          if (char_index lt last_char_index) then begin
              char_index=char_index+1
              test_char=strmid(compound_string,char_index,1)
          endif else begin
              test_char=''
          endelse
      endwhile
      
      ;; Did we find a match?
      if (nest_level eq 0) then begin
          paren_end_index = char_index-1
      endif else begin
          print,'zcompound: could not find matching parenthesis in ', $
            compound_string
          z_array = 0
          return
      endelse
      
      ;; print,'paren_end_index=',paren_end_index,', nest_level=',nest_level
      
      ;; Do any multipliers come after the matching parenthesis?
      ;; If this is the end of the string, then evaluate what's
      ;; inside with no multiplier
      num_start_index = paren_end_index+1
      if (num_start_index gt last_char_index) then begin
          in_paren_multiplier=1.
          num_start_index=0
          postnum_index=0
      endif else begin
          ;; Find the next element or parenthesis, as
          ;; anything before it must be a number.
          postnum_index = num_start_index
          test_char=strmid(compound_string,postnum_index,1)
          while ( ((test_char eq '0') or (test_char eq '1') or $
                   (test_char eq '2') or (test_char eq '2') or $
                   (test_char eq '3') or (test_char eq '4') or $
                   (test_char eq '5') or (test_char eq '6') or $
                   (test_char eq '7') or (test_char eq '8') or $
                   (test_char eq '9') or (test_char eq '.')) and $
                  (postnum_index le last_char_index) ) do begin
              postnum_index=postnum_index+1
              if (postnum_index le last_char_index) then begin
                  test_char=strmid(compound_string, $
                                   postnum_index,1)
              endif else begin
                  test_char=''
              endelse
          endwhile
          
          ;; Got the number following parenthesis
          if (num_start_index ne postnum_index) then begin
              number_string=strmid(compound_string, $
                                   num_start_index, $
                                   (postnum_index-num_start_index))
              in_paren_multiplier=1.
              ;; print,'Trying to interpret ',number_string,' as a number.'
              if (strlen(number_string) ne 0) then begin
                  reads,number_string,in_paren_multiplier
              endif
          endif else begin
              in_paren_multiplier=1.
          endelse
      endelse
      
      ;; Factor in this multiplier, and crunch on what's
      ;; in the parenthesis.
      net_multiplier = in_paren_multiplier * paren_multiplier
      in_paren_string=strmid(compound_string,1,(paren_end_index-1))
      ;; print,'in_paren_string=',in_paren_string,', net_multiplier=',net_multiplier
      zcompound,in_paren_string,z_array, $
        paren_multiplier=net_multiplier
      
      ;; print,'postnum_index=',postnum_index
      ;; Crunch on whatever follows the multiplier. */
      if (postnum_index le last_char_index) then begin
          remaining_string=strmid(compound_string,postnum_index, $
                                  (last_char_index-postnum_index+1))
          ;; print,'remaining_string=',remaining_string
          zcompound,remaining_string,z_array, $
            paren_multiplier=paren_multiplier
      endif
  endelse
  
  return
END

