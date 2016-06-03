;+
; NAME:
;       CONSEC
;
; PURPOSE:
;       Finds sets of consecutive numbers in 
;       an array.
;
; CALLING SEQUENCE:
;       CONSEC, ARRAY, LO, HI [, NUM, /SAME]
;
; INPUTS:
;       ARRAY:  Umm, the array.
;
; KEYWORDS:
;       SAME:   Search for elements which are the same, rather than 
;               consecutive.
; OUTPUTS:
;       LO: The lower index limit of each set 
;           of consecutive numbers. Each set
;           goes from LO_i to HI_i
;           
;       HI: The upper index limit  
;
; OPTIONAL OUTPUTS:
;       NUM:          The number of sets.
;       DISTRIBUTION: Array containing the number of elements in each
;                     set
;
; OPTIONAL KEYWORDS:
;       SAME: Search for consecutive elements containing the same value.
;
; EXAMPLE:
;
; IDL> a = [findgen(3),findgen(3)+6]
; IDL> print,a
;   0.00000    1.00000    2.00000    6.00000    7.00000    8.00000
;
; IDL> consec,a,l,h,n
; IDL> print,l
;       0       3
; IDL> print,h
;       2       5
; IDL> print,n
;       2
; IDL> print,a[l[0]:h[0]]
;      0.00000      1.00000      2.00000
;
; MODIFICATION HISTORY:
;  Written in January 2003 by JohnJohn
;  5-13-2003 JohnJohn   Added SAME keyword
;  4-16-2005 JohnJohn   Found a bug where -1 was returned if there
;  were only two identical elemeents and /SAME. New logic introduced
;  to deal with this special case. Also added DISTRIBUTION keyword.
;
;-

pro consec, a, l, h, n, same=same, distribution=dist
on_error,2
nel = n_elements(a)
case nel of
    1: begin
        l = 0
        h = 0
        n = 1
    end
    2:  begin
        if keyword_set(same) then begin
            if a[1] - a[0] eq 0 then begin
                l = 0
                h = 1
                n = 1
            endif else begin
                l = -1
                h = -1
                n = 0
            endelse
        endif else begin
            if abs(a[1] - a[0]) eq 1 then begin
                l = 0
                h = 1
                n = 1
            endif else begin
                l = -1
                h = -1
                n = 0
            endelse
        endelse
    end
    else: begin
        if not keyword_set(same) then begin
            arr = [a[0],a,a[nel-1]] ;add fat
            cond1 = abs(arr - shift(arr,1)) eq 1
            cond2 = abs(arr - shift(arr,-1)) eq 1
            range = indgen(nel)+1 ;trim fat
        endif else begin
            arr = [a[0]+1,a,a[nel-1]-1] ;add fat
            cond1 = abs(arr - shift(arr,1)) eq 0
            cond2 = abs(arr - shift(arr,-1)) eq 0
            range = indgen(nel)+1 ;trim fat
        endelse
        l = where(cond2[range] and not cond1[range], nl)
        h = where(cond1[range] and not cond2[range], nh)
        if nh*nl eq 0 then begin 
            l = -1
            h = -1
            n = 0 
        endif else n = nh < nl
    end
endcase
if l[0] ne h[0] then dist = h-l+1 else dist = 0

end