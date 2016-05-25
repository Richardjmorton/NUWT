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

pro patch_up,fit_flag=fit_flag

;Loads common data generated from locate_things_fg
COMMON located_dat, located, nx, nt
COMMON threads_dat, threads

IF n_elements(fit_flag) EQ 0 THEN fit_flag=0
sz=size(threads)
n_threads=sz(1)

IF fit_flag EQ 0 THEN BEGIN
;Designed for simple fitting structures
     dummy=threads.pos
     dummyerrs=threads.err_pos
ENDIF ELSE BEGIN

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
    
ENDELSE



FOR j=0,(n_threads-1) DO BEGIN
   ;skips enteries with less than 2 positive values
   IF (n_elements(where(dummy[*,j] gt 0.))) GT 2. THEN BEGIN

      ;skips enteries where half the data points are set to zero, i.e. no
      ;value was obtained at fitting stage.
      IF (n_elements(where(dummy[*,j] EQ 0.))) LT 0.5*(n_elements(where(dummy[*,j] GE 0.))) THEN BEGIN
          
          ;FIND CONSECUTIVE ZEROS IN THREAD
          b=where(dummy[*,j] EQ 0)
          IF b[0] GT -1 THEN consec,b,lo,hi,num

          ;FILL IN USING LINEAR INTERPOLATION
          IF b[0] GT -1 THEN FOR i=0,num-1 DO linfill,dummy[*,j],b(lo[i])-1,b(hi[i])+1
          IF b[0] GT -1 THEN FOR i=0,num-1 DO dummyerrs[lo[i]:hi[i],j]=1.

          ;FIND SINGLE ZEROS AND FILL VIA LINEAR INTERPOLATION
          b=where(dummy[*,j] EQ 0)
          IF b[0] GT -1 THEN FOR i=0,n_elements(b)-1 DO $
          dummy[b[i],j,0]=0.5*dummy[b[i]-1,j]+0.5*dummy[b[i]+1,j]
          IF b[0] GT -1 THEN FOR i=0,n_elements(b)-1 DO $
          dummyerrs[b[i],j]=sqrt(dummy[b[i]-1,j]^2+dummy[b[i]+1,j]^2)
          
          
        
      ENDIF ELSE BEGIN
          ;Erases all entries where number of 0 elements gt 1/2 positive elements
          dummy[*,j]=-1.
      ENDELSE
 ENDIF ELSE BEGIN
   ;Erases all entries where number of positive elements lt 2
   dummy[*,j]=-1.
 ENDELSE
ENDFOR




END
