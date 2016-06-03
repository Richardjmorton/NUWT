;+
; NAME:
;       STR
;
;
; PURPOSE:
;       Convert numerical values to a string and trim off
;       white space. Works better than IDL's STRING for quick and
;       dirty string formatting.
;
;
; CALLING SEQUENCE:
;       result = str(num)
;
;
; INPUTS:
;       num:   Scalar or vector numbers to be converted to string 
;               values.
;
;
; KEYWORD PARAMETERS:
;       LENGTH:  For integer values, this specifies the total string
;       length of the returned value including leading or trailing
;       characters. 
;
;       FORMAT: Conventional Fortran-style formatting string
;
;       TRAIL: By defaul, characters are added to the front of integer
;       values, e.g. str(1, char='0', length=3) returns 001
;
;       CHARACTER: The fill character to be used. Default is '0'
;
; OUTPUTS:
;       RESULT: The resulting string
;
; EXAMPLE:
;
;      IDL> print,str(1)                             
;      1
;      IDL> print,str(1, char='0', length=3)         
;      001
;      IDL> print,str(1, char='*', length=3, /trail)
;      1**
;      IDL> print,str(1.0, format='(f3.1)')
;      1.0
;
; MODIFICATION HISTORY:
; Created sometime in 2003 by JohnJohn
;-

function str,number,length=length,format=format,trail=trail,character=char
on_error,2
n = n_elements(number)
s = strtrim(string(number,format=format),2)
if 1-keyword_set(char) then char = '0'
if n_elements(length) gt 0 then begin
    ilen = strlen(s)
    for i = 0,n-1 do begin
        nz = length-ilen[i]
        if nz gt 0 then begin
            for j=0,nz-1 do begin
                if keyword_set(trail) then s[i] = s[i]+char else $
                  s[i] = char+s[i] 
            endfor
        endif
    endfor
    endif
return,s
end
