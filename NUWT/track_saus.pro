;PURPOSE: Locates thread in x-t diagrams, then tries to fit a Gaussian function to cross-sectional flux profile. 
;It then fits sinusoidal functions to central positions of a thread, thread width & thread intensity. Designed for measuring physical 
;displacements of structure (kink) and width/intensity variations (sausage).
;
;
;
;INPUTS: data - xt array
;OUTPUTS: Can only be saved to file. The results of fitting each slit are saved to an .idlsave with a name of your choice.
;         Either specify naming system once at start (fname keyword) or each time (don't use fname keyword).
;
;OPTIONAL INPUTS:
;                 errors - array of errors on intensity, i.e. photon noise, 
;                          has to be same size as data array 
;                 algn_err - provide an estimate for error due to alignment, if given it 
;                            is added to uncertainty in Gaussian fit 
;                 meas_size - sets number of pixels used in Gaussian fit, default is 5, 
;                             odd numbers only                             
;                 damped - set if fit of damped sine function is required
;                 tdp    - set to fit a time-dependent damped sine function (do not recommend this option!!)
;                 area  -  changes size of search box used in follow_thread, default is 5
;                          highly recommended to leave unchanged
;                 fname - directory and filename to save results to - input as directory/filename
; 		   check - plots gaussian fitting of cross-cuts to screen
;                 patch - Performs a linear interpolation to fill in gaps in threads due to time jumps
;
;CALLS: locate_things_fg.pro, follow_thread_fg.pro, moscill.pro
;
;HISTORY - created R Morton APR 2014
;         R Morton NOV 2014 - super update! Added structure format to remove all the arrays.
;                            Also added COMMON variables so re-used values/structures are passed automatically
;         
;
;
;TO DO: Update plot section to use new structure
;
;

pro track_saus, data,algn_err=algn_err,meas_size=meas_size,$
                errors=errors,damped=damped,area=area,tdp=tdp,$
                fname=fname,check=check,cut_chisq=cut_chisq,patch=patch


COMMON located_dat, located, nx, nt
COMMON threads_dat, threads
COMMON extended_dat, threads_fit_fg

;#############################################
;Initial set up of variables
;#############################################
sz=size(data)
IF sz(0) EQ 3 THEN nslit=sz(3) ELSE nslit=1
IF n_elements(meas_size) EQ 0 THEN meas_size=9
IF n_elements(errors) EQ 0 THEN stop,'Uncertainties in intensities needed - errors keyword.'
IF n_elements(algn_err) GT 0 THEN print,'Error on alignment is:', algn_err


numg=1d
READ,numg,PROMPT='Enter gradient: '

numt=1d
READ,numt,PROMPT='Enter minimum thread length: '


;#############################################
;Begin search and collect routine
;#############################################
FOR i=0, nslit-1 DO BEGIN


    print,'#####################'
    print,'Slit number: ',strtrim(i+1,2), ' of ',strtrim(nslit,2)
    print,'#####################'
    
   ;#############################################
   ;Locate peaks and fit Gaussian
   ;#############################################
    IF n_elements(errors) gt 0. THEN errorsi=errors[*,*,i] 
    locate_things_FG,data=data[*,*,i],grad=numg,meas_size=meas_size, $
                 errors=errorsi,check=check,cut_chisq=cut_chisq
   
   ;++++++++++++++++++++++++++++++++++++++++++++
   ;Add alignment errors to error from Gaussian fits
   ;++++++++++++++++++++++++++++++++++++++++++++
    IF n_elements(algn_err) GT 0 THEN located.errs=located.errs+algn_err   ; adds alignment error to error from Gaussian fits
    
    window,0
    !p.multi=[0,2,1]
    tvim,data[*,*,i]
    tvim,located.peaks[*,*,1]
       
    window,1
    !p.multi=[0,1,1]
   ;++++++++++++++++++++++++++++++++++++++++++++
   ;Connect fitted peaks
   ;++++++++++++++++++++++++++++++++++++++++++++
    follow_thread_fg,min_tlen=numt,area=area
    print,'Threads found: '+strtrim((size(threads))(1)-1,2)

   ;Liner interpolation fill of the gaps
   IF keyword_set(patch) THEN patch_up,fit_flag=1
    
   dsk='y'
   READ,dsk,PROMPT='Do you want to fit the transverse motions?: (y/n) '
   IF dsk EQ 'y' THEN back=fitting_them(1,damped=damped,tdp=tdp)
   
   dsk='y'
   READ,dsk,PROMPT='Do you want to fit the intensity?: (y/n) '
   IF dsk EQ 'y' THEN back=fitting_them(2,damped=damped,tdp=tdp)
    
   dsk='y'
   READ,dsk,PROMPT='Do you want to fit the width?: (y/n) '
   IF dsk EQ 'y' THEN back=fitting_them(3,damped=damped,tdp=tdp)
      

   yn='a'
   READ,yn,PROMPT='Save results? If the results are not saved now they will be overwritten! y/n '
   IF yn EQ 'y' THEN BEGIN

      IF n_elements(fname) EQ 0 THEN BEGIN
         sv_dir='a'
       	 READ,sv_dir,PROMPT='Give directory + filename for saving data points and fits. Saved as directory/filename_N.idl: '
      ENDIF ELSE sv_dir=fname

      IF n_elements(threads_fit_fg) EQ 0 THEN save,threads,filename=sv_dir+'_'+strtrim(i,2)+'.idl' $
      ELSE save,threads,threads_fit_fg,filename=sv_dir+'_'+strtrim(i,2)+'.idl'
      

    ENDIF

ENDFOR

END
