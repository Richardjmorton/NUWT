


;Makes a diagonal line across data and reads data into an
;array of dimensions max(x,y) for all time
;
;INPUTS: data - NxMxT array with data in
;
;OPTIONAL INPUTS - x1,x2,y1,y2 - coordinates at either end of diagonal line    
;                  fn - frame number to plot
;                  xsize, ysize - size of graphics window (command for widgets)
;                  win_id - set to value of open graphics window you want to use. If not set, opens
;                           next free window.
;                  noopen - if set suppresses the opening of windows
;
;OUTPUTS: data_o - output array where data along diagonal line is stored
;         outvec - array of coordinates used for slit
;
;HISTORY: Created Richard Morton 2011
;         RM 2014 - Added moving plot line
;         RM NOV 2014 - added some graphics set up to use with widgets
;         RM APR 2015 - fixed rounding error with fix command
;         RM MAY 2015 - Introduced interpolation routine to extract data points - made diag_slit2 redundant!
;
;TO DO: Is current routine to set up windows making a meal of things?
;

PRO diag_slit, data,data_o,x1=x1,x2=x2,y1=y1,y2=y2,debug=debug,fn=fn,outvec=outvec,$
               xsize=xsize,ysize=ysize,win_id=win_id,noopen=noopen


IF NOT keyword_set(noopen) THEN BEGIN
   IF n_elements(fn) eq 0 THEN fn=0

   ;Set plotting elements
   ;writing to buffer pixmap reduces on-screen jitter!
   ;
   IF n_elements(win_id) EQ 0 THEN BEGIN
      determine_window_size,xsize=xsize,ysize=ysize,openwin=openwin,/new_win 
      window,openwin[0],xs=xsize,ys=ysize
      window_var=openwin[0]
   ENDIF ELSE BEGIN 
      wset,win_id
      window_var=win_id
   ENDELSE
   window,xs=xsize,ys=ysize,/pixmap,/free
   pixID=!d.window
   tvim,data[*,*,fn]
   wset,window_var
   device,copy=[0,0,!d.x_size,!d.y_size,0,0,pixid]
ENDIF

IF n_elements(x1) EQ 0. THEN BEGIN
   PRINT,'Select point with cursor'
   CURSOR,x1,y1,/Data,/Down
   print,  'x1='+strtrim(string(fix(x1+0.5)),1), ' y1='+strtrim(string(fix(y1+0.5)),1)
   PRINT,'Select point with cursor'

   IF !MOUSE.button EQ 1 THEN BEGIN
      CURSOR, do_nothing_x, do_nothing_y, /UP, /WAIT 
      !MOUSE.button=0
   ENDIF
   WHILE (!MOUSE.button NE 1) DO BEGIN
        CURSOR,x2,y2,/change
        window,xs=xsize,ys=ysize,/pixmap,/free
        pixID=!d.window
        tvim,data[*,*,fn]
        wset,window_var
        device,copy=[0,0,!d.x_size,!d.y_size,0,0,pixid]
        PLOTS,[x1,x2],[y1,y2]
   ENDWHILE
   print,  'x2='+strtrim(string(fix(x2+0.5)),1), ' y2='+strtrim(string(fix(y2+0.5)),1)

   x1=1.*fix(x1+0.5)
   x2=1.*fix(x2+0.5)
   y1=1.*fix(y1+0.5)
   y2=1.*fix(y2+0.5)

ENDIF ELSE BEGIN
   x1=1.*x1 & x2=1.*x2 & y1=1.*y1 & y2=1.*y2
ENDELSE

zsize=n_elements(data[0,0,*])
xsize=abs(x2-x1)
ysize=abs(y2-y1)


IF NOT keyword_set(noopen) THEN plots,[x1,x2],[y1,y2]

IF keyword_set(debug) THEN BEGIN
  IF x1 lt x2 THEN BEGIN
     x=x1+findgen(xsize)
     oplot,x,m*x+c,psym=1,nsum=3
  ENDIF ELSE BEGIN
     x=x1-findgen(xsize)
     oplot,x,m*x+c,psym=1,nsum=3
  ENDELSE
ENDIF


slit_lnth=round(sqrt(ysize^2+xsize^2))
mindat=min(data)
IF mindat GT 0 THEN mindat=-1. ELSE mindat=mindat-0.1*sqrt((moment(data))[1])

IF (xsize GT 0) AND (ysize GT 0) THEN BEGIN 

; Calculating gradient and intercept
m=(y2-y1)/((x2-x1))
c=-x1*(y2-y1)/(x2-x1)+y1

   IF xsize GE ysize THEN BEGIN
      data_o=fltarr(slit_lnth,zsize)

      FOR i=0,(zsize-1) DO BEGIN
          x=min([x1,x2])+findgen(xsize)
          y=m*x+c
          x=congrid(temporary(x),slit_lnth,cubic=-0.5)
	  y=congrid(temporary(y),slit_lnth,cubic=-0.5)
          IF y1 GT y2 THEN data_o[*,i]=reverse(interpolate(reform(data[*,*,i]),x,y,cubic=-0.5,missing=mindat)) $
          ELSE data_o[*,i]=interpolate(reform(data[*,*,i]),x,y,cubic=-0.5,missing=mindat)
		
      ENDFOR
    ENDIF

    IF ysize GE xsize THEN BEGIN
       data_o=fltarr(slit_lnth,zsize)

       FOR i=0,(zsize-1) DO BEGIN
         y=min([y1,y2])+findgen(ysize)
         x=(y-c)/m
	 x=congrid(temporary(x),slit_lnth,cubic=-0.5)
	 y=congrid(temporary(y),slit_lnth,cubic=-0.5)
         IF x1 GT x2 THEN data_o[*,i]=reverse(interpolate(reform(data[*,*,i]),x,y,cubic=-0.5,missing=mindat)) $
         ELSE data_o[*,i]=interpolate(reform(data[*,*,i]),x,y,cubic=-0.5,missing=mindat)
          
       ENDFOR
    ENDIF
  

ENDIF ELSE BEGIN
  
  IF xsize EQ 0 THEN BEGIN
     
     dum=[y1,y2]
     y1=min(dum) & y2=max(dum)
     data_o=fltarr(ysize,zsize)
     data_o[*,*]=reform(data[x1,y1:y2-1,*])
  ENDIF
  IF ysize EQ 0 THEN BEGIN
     dum=[x1,x2]
     x1=min(dum) & x2=max(dum)
     data_o=fltarr(xsize,zsize)
     data_o[*,*]=reform(data[x1:x2-1,y1,*])
  ENDIF 
ENDELSE

outvec=[x1,x2,y1,y2]

END


