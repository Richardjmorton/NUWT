;FUNCTION: 'Patchs up' the threads found from 'locate_things.pro', followed by
;          'follow_things.pro'. Follow things can skip forward time frames
;           leaving zero results. This uses interpolation to fix these gaps and 
;           weights these results with a larger error.
;

 

;PROCEDURE OUTLINE: The first step is to remove any examples where
;               the number of pixels in the thread is less than 2 and
;               any threads where more than half the pixels could not
;               be found with locate things (i.e. given the value 0). The routine then 
;               locates any points that have the value 0. These values are then 
;               replaced via linear interpolations given a large error (1
;               pixel). This is to ensure the weighted fitting
;               routines give limited significance to this value. 
;               
;                 
;
;INPUTS: fit_flag - 0-3 - defines which time series to work on (see follow_thread_fg.pro/moscill.pro)
;
;
;HISTORY: created 05/2014 R Morton
;
;TO DO OR THINK ABOUT: Is linear interpolation the best option?
;                      Should weighting be calculated using traditional error analysis?
;                     Does funny things if first elements in array is a 0 - need to look into this
;
;

pro patch_up_SA,threads=threads,fit_flag=fit_flag


;Designed for advanced fitting structures (_FG)
    IF fit_flag EQ 1 THEN BEGIN 
        dummy=threads.pos
        dummyerrs=threads.err_pos
    ENDIF
    IF fit_flag EQ 2 THEN BEGIN 
       dummy=threads.inten
       dummyerrs=threads.err_inten
       
    ENDIF
    IF fit_flag EQ 3 THEN BEGIN 
       dummy=threads.wid ;dummy arrays
       dummyerrs=threads.err_wid
    ENDIF
    


IF (size(dummy))[0] GT 1 THEN n_threads=(size(dummy))[2] ELSE n_threads=1

FOR j=0,(n_threads-1) DO BEGIN
   ;skips enteries with less than 2 positive values
   IF (n_elements(where(dummy[*,j] gt 0.))) GT 2. THEN BEGIN
            
          ;FIND CONSECUTIVE ZEROS IN THREAD
          b=where(dummy[*,j] EQ 0)
          IF b[0] GT -1 THEN consec,b,lo,hi,num
          plot,dummy[*,j]
          ;FILL IN USING LINEAR INTERPOLATION
          IF b[0] GT -1 THEN FOR i=0,num-1 DO BEGIN
                                  dum2=reform(dummy[*,j])
                                  linfill,dum2,b[lo[i]]-1,b[hi[i]]+1
				  dummy[*,j]=dum2
			     ENDFOR

          IF b[0] GT -1 THEN FOR i=0,num-1 DO dummyerrs[b(lo[i]):b(hi[i]),j]=1.
          
          
          ;FIND SINGLE ZEROS AND FILL VIA LINEAR INTERPOLATION
          b=where(dummy[*,j] EQ 0)
          IF b[0] GT -1 THEN FOR i=0,n_elements(b)-1 DO $
          dummy[b[i],j,0]=0.5*dummy[b[i]-1,j]+0.5*dummy[b[i]+1,j]
          IF b[0] GT -1 THEN FOR i=0,n_elements(b)-1 DO $
          dummyerrs[b[i],j]=1.
          oplot,dummy[*,j]
   ENDIF       
ENDFOR


IF fit_flag EQ 0 THEN BEGIN
;Designed for simple fitting structures
     threads.pos=dummy
     threads.err_pos=dummyerrs
ENDIF ELSE BEGIN

;Designed for advanced fitting structures (_FG)
    IF fit_flag EQ 1 THEN BEGIN 
       threads.pos=dummy
       threads.err_pos=dummyerrs
    ENDIF
    IF fit_flag EQ 2 THEN BEGIN 
       threads.inten=dummy
       threads.err_inten=dummyerrs
       
    ENDIF
    IF fit_flag EQ 3 THEN BEGIN 
       threads.wid=dummy ;dummy arrays
       threads.err_wid=dummyerrs
    ENDIF
    
ENDELSE

END
