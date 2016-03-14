
;PURPOSE: Uses locate things to examine a time-distance diagram and
;         plot peaks. Mainly to be used for checking input parameters and 
;         plotting. 
;INPUTS: data - xt array
;
;OPTIONAL INPUTS:
;                 gauss - set keyword /gauss to enable sub-pixel location needs errors
;                 errors - array of uncertainties on intensity, e.g., photon noise, read noise, etc., 
;                          has to be same size as data array  
;                 meas_size - sets number of pixels used in Gaussian fit, default is 5, 
;                             odd numbers only
;                 cut_chisq - Controls cut off value for reduced chi squared, default 5. I.e. if reduced chi squared of the fit
;                             is greater than 5 then gaussian fit is ignored.
;                                              
;                 doplot - plot to file instead of screen
;                 savedir - directory and filename to save results to - input as directory/filename
; 		  check - plots gaussian fitting of cross-cuts to screen
;		  dx - pixel size
;		  dt - cadence
;
;
;
;HISTORY - created R Morton 2014
;          update  MAR 2016 R Morton added greater functionality
;
;



PRO quick_plot_ft,data,gauss=gauss,errors=errors,dx=dx,dt=dt,doplot=doplot,  $ 
                      savedir=savedir,check=check,meas_size=meas_size,cut_chisq=cut_chisq 

on_error,2
IF n_elements(data) EQ 0 THEN message,'You forgot to add some data!!'
IF keyword_set(savedir) THEN file=savedir ELSE file='~/'
IF keyword_set(gauss) THEN IF n_elements(errors) EQ 0 THEN message,'Uncertainties in intensities needed - errors keyword.'

;###############################################
;INITIAL HOUSE KEEPING
;###############################################

COMMON located_dat, located, nx, nt
COMMON threads_dat, threads

on_error,2
sz=size(data)
IF sz(0) EQ 3 THEN nslits=sz(3) ELSE nslits=1
nx=sz(1) & nt=sz(2)

IF NOT keyword_set(dx) THEN dx=1.
IF NOT keyword_set(dt) THEN dt=1.

IF n_elements(meas_size) EQ 0 THEN meas_size=7



IF NOT keyword_set(gauss) THEN print, 'Finding whole pixel values' ELSE print, 'Gaussian location enabled'

numg=1d
READ,numg,PROMPT='Enter gradient: '



;###############################################
;BEGIN WORKING ON EACH TIME-DISTANCE DIAGRAM
;###############################################

FOR k=0, nslits-1 DO BEGIN

    h=0 ; Counting device
    
    print,'&&&&&&&','Slit number:', k
      
    ;Fitting of intensity maxima
    IF not keyword_set(gauss) THEN BEGIN
      locate_things_min,data=data[*,*,k],grad=numg
    ENDIF ELSE BEGIN

         errorsi=errors[*,*,k]   
         locate_things_fg,data=data[*,*,k],grad=numg,meas_size=meas_size, $
               errors=errorsi,check=check,cut_chisq=cut_chisq 
    ENDELSE

    ;Plot found data points
    IF NOT keyword_set(doplot) THEN BEGIN
          window,0
          !p.multi=[0,2,1]

          tvim,data[*,*,k],xrange=[0,nx-1]*dx/1.e3,yrange=[0,nt-1]*dt,xtitle='Distance (Mm)',ytitle='Time (s)'
          tvim,located.peaks[*,*,1],xrange=[0,nx-1]*dx/1.e3,yrange=[0,nt-1]*dt,xtitle='Distance (Mm)',ytitle='Time (s)',title='Red crosses show non-Gaussian results'

          ;Overplots values where Gaussian fitting failed with red box
          
          in=where(located.errs[*,*,0] eq 0.5)
          IF in[0] ne -1 THEN BEGIN 	     
          	dum=located.errs[*,*,0]
          	dumx=(size(dum))[1]
                loadct,13,/silent
          	oplot, (in mod dumx)*dx/1.e3, (in/dumx)*dt, psym=1,thick=1,color=255
          	loadct,0,/silent
          ENDIF
          
          !p.multi=0
    ENDIF ELSE BEGIN
        set_plot,'ps'
        device,/encapsul,/color,filename=savedir+'slit'+strtrim(k,2)+'.eps'
        !p.multi=[0,1,2]
        med=median(data[*,*,k])
        sd=sqrt((moment(data[*,*,k]))[1])
        tvim,rotate(data[*,*,k],3)<(med+3.*sd)>(med-3.*sd),yrange=[0,nx-1]*dx/1.e3,xrange=[0,nt-1]*dt,ytitle='Distance (Mm)',xtitle='Time (s)',pcharsize=0.6
        image=located.peaks[*,*,0]
        image=image gt min(image)
        tvim,rotate(image,3),yrange=[0,nx-1]*dx/1.e3,xrange=[0,nt-1]*dt,ytitle='Distance (Mm)',xtitle='Time (s)',pcharsize=0.6
        device,/close
        !p.multi=0
    ENDELSE
ENDFOR

END