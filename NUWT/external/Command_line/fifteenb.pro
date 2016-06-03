;+
; NAME:
;      FIFTEENB
;
;
; PURPOSE:
;      Return the string version of the octal character
;      "15b which print's a new line to the screen. Used
;      by COUNTER.PRO This is a separate routine because the 
;      unclosed quotation mark messes up IDLWAVE in Emacs, and 
;      I can't have that!
;
; CATEGORY:
;
;
;
; CALLING SEQUENCE:
;
;      result = fifteenb()
;
; INPUTS:
;
;
;
; OPTIONAL INPUTS:
;
;
;
; KEYWORD PARAMETERS:
;
;
;
; OUTPUTS:
;
;      result is set to string("15b)
;
; OPTIONAL OUTPUTS:
;
;
;
; COMMON BLOCKS:
;
;
;
; SIDE EFFECTS:
;
;
;
; RESTRICTIONS:
;
;
;
; PROCEDURE:
;
;
;
; EXAMPLE:
;
;
;
; MODIFICATION HISTORY:
;
;-

function fifteenb
return,string("15b)
end