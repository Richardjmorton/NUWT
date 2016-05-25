;PURPOSE: Locates thread in x-t diagrams with sub pixel accuracy (if Gaussian fitting 
;option is selected) and fits sinusoidal functions to oscillations. Is essentially a 
;wrapper for called routines. Designed for measuring physical displacements of structure.
;
;
;
;INPUTS: data - xt array
;OUTPUTS: results - array of results
;         data_points - array of data points and error bars
;
;OPTIONAL INPUTS: /gauss - If set, performs Gaussian fit around maximum
;                          to find position with sub-pixel accuracy.
;                          If not set, whole pixel value of maximum taken
;                 algn_err - provide an estimate for error due to alignment, if given it 
;                            is add to uncertainty in Gaussian fit 
;                 meas_size - sets number of pixels used in Gaussian fit, default is 5, 
;                             odd numbers only                             
;                 errors - array of errors on intensity, i.e. photon noise, 
;                          has to be same size as data array
;                 damped - set if fit of damped sine function is required
;                 area  -  changes size of search box used in follow_thread, default is 5
;                          highly recommended to leave unchanged
;                 fname -  provide the directory to save to and the naming convention, i.e. directory/file
;                 simp_grad - Use of an analytic least-sqaures estimate for gradient - may be faster!  
;
;CALLS: locate_things.pro, locate_things_min.pro, follow_thread.pro, moscill.pro
;
;HISTORY - created R Morton OCT 2012
;          updated R Morton FEB 2013 - incorporated upgrades to locate_things.pro 
;                  R Morton MAR 2013 - incorporated optional fit of damped sin function
;                  R Morton APR 2013 - added area keyword
;
;                  R Morton APR 2014 - updated follow_thread.pro
;                  R Morton NOV 2014 - super update! Added structure format to remove all the arrays.
;                                      Also added COMMON variables so re-used values/structures are passed automatically
;                  R Morton MAR 2016 - Cleaned up!
;                  14 MAR 2016 - Release as version 1.0 of NUWT 
;                  R Morton MAY 2016 - minor changes to sub routines
;
;
;TO DO: - Make plot ready for plotting damped profiles + tdp
;       - area keyword still buggy and certain values cause the routine to crash
;      
;

pro track_rout, data,results=results,gauss=gauss,$ 
                algn_err=algn_err,meas_size=meas_size,errors=errors,$
                damped=damped,area=area,tdp=tdp,$
                fname=fname,check=check,cut_chisq=cut_chisq,$
                simp_grad=simp_grad

COMMON located_dat,located, nx, nt
COMMON threads_dat, threads
COMMON extended_dat, threads_fit_fg

sz=size(data)
IF sz(0) EQ 3 THEN nslit=sz(3) ELSE nslit=1

IF n_elements(algn_err) GT 0 THEN print,'Error on alignment is:', algn_err


numg=1d
READ,numg,PROMPT='Enter gradient: '

numt=1d
READ,numt,PROMPT='Enter minimum thread length: '



FOR i=0, nslit-1 DO BEGIN

    print,'#####################'
    print,'Slit number: '+strtrim(i+1,2)+' of '+strtrim(nslit,2)
    print,'#####################'
    IF not keyword_set(gauss) THEN BEGIN
      ;output is in common located_dat
      IF n_elements(errors) GT 0 THEN inerr=errors[*,*,i] ELSE inerr=0
      locate_things_min,data=data[*,*,i],grad=numg,errors=inerr,simp_grad=simp_grad
      print,'no Gaussian fit used'
    ENDIF ELSE BEGIN
         
      locate_things,data=data[*,*,i],grad=numg,meas_size=meas_size, $
                 errors=errors[*,*,i],check=check,cut_chisq=cut_chisq,simp_grad=simp_grad
    
      IF n_elements(algn_err) GT 0 THEN located.errs=located.errs+algn_err   ; adds alignment error to error from Gaussian fits

      print,'Gaussian fit used'
    ENDELSE

    window,0
    !p.multi=[0,2,1]
    tvim,data[*,*,i]
    tvim,located.peaks[*,*,1]
   
    
    ;output is in common thread_dat
    follow_thread,min_tlen=numt,area=area

    
   dsk='y'
   READ,dsk,PROMPT='Do you want to fit the transverse motions?: (y/n) '
   IF dsk EQ 'y' THEN BEGIN

      window,1
      !p.multi=[0,1,1]

      IF NOT keyword_set(damped) THEN BEGIN
         moscill,out=threads_fit
      ENDIF ELSE BEGIN
         IF NOT keyword_set(tdp) THEN BEGIN
 	        moscill,/tdp,out=threads_fit
         ENDIF ELSE BEGIN
               moscill,/damped,/tdp,out=threads_fit
         ENDELSE
      ENDELSE

      res=where(threads_fit[*].fit_result[1] ne 0)
      nfits=n_elements(res)
      print,'No. fits saved',nfits
   ENDIF

   
    yn='a'
    READ,yn,PROMPT='Save results? If the results are not saved now they will be overwritten! y/n '
    IF yn EQ 'y' THEN BEGIN
 
       IF n_elements(fname) EQ 0 THEN BEGIN
	  sv_dir='a'
       	   READ,sv_dir,PROMPT='Give directory + filename for saving data points and fits. Saved as directory/filename_N.idl: '
       ENDIF ELSE sv_dir=fname
            
      IF n_elements(threads_fit) EQ 0 THEN save,threads,filename=sv_dir+'_'+strtrim(i,2)+'.idl' $
      ELSE save,threads,threads_fit,filename=sv_dir+'_'+strtrim(i,2)+'.idl'

    ENDIF

ENDFOR

END
