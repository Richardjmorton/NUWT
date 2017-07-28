;+
; ROUTINE:    pick
;
; PURPOSE:    Short version of cursor,x,y,/data,/down
;
; USEAGE:     pick,x,y
;
; INPUT:
;
; OPTIONAL INPUT: device - uses device coordinates
;                 doplot - plot cross at picked position
;                 noprint - suppress command line text
;
; OUTPUT:    x and y coords of chosen point
 
; Example:    
;             
; AUTHOR:     Peter T. Gallagher, May. '98
;             RJM added device, noprint and doplot keywords
;-

PRO pick, x, y,device=device,doplot=doplot,noprint=noprint 
    
 IF NOT keyword_set(noprint) THEN PRINT,'Click on active window' 
 
 IF NOT keyword_set(device) THEN CURSOR,x,y,/Data,/Down $
 ELSE CURSOR,x,y,/Data,/Down,/device
 
 IF NOT keyword_set(noprint) THEN PRINT, x,y

 IF keyword_set(doplot) THEN $
 IF NOT keyword_set(device) THEN plots,[x,y],psym=1 ELSE plots,[x,y],/device,psym=1
 
    
END
