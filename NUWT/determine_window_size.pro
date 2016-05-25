;To be used for advanced use of tvim.pro
;Determines scaled size of current graphic window
;or opens first free window and retrieves scaled window size
;
;OPTIONAL INPUTS: new_win - check to open new window otherwise currently open
;                           window is used  
;
;OUTPUTS: xsize, ysize - Obvious
;         openwin - array of windows that are un-open
;
;

PRO determine_window_size,xsize=xsize,ysize=ysize,openwin=openwin,new_win=new_win

  IF keyword_set(new_win) THEN BEGIN
     device,window_state=these_windows
     openwin=where(these_windows eq 0)
     window,openwin[0]

  ENDIF ELSE wset,!d.window

  plot, [0,1],[0,1],/nodata,xstyle=4,ystyle=4,charsize=pcharsize
  px=!x.window*!d.x_vsize
  py=!y.window*!d.y_vsize
  xsize=px(1)-px(0)
  ysize=py(1)-py(0)


END