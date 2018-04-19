;
;NUWT Version 2.0
;
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
;OPTIONAL INPUTS: /nearest_pixel - if set, will ALWAYS return the nearest whole-pixel 
;                                  location of the peaks without bothering with any
;                                  Gaussian fitting
;                 /invert - If set, will invert the intensity values of the input data 
;                           before finding peak locations.This results in the code 
;                           finding the local MINIMA in the original image.
;
;                 errors - array of errors on intensity, i.e. photon noise, 
;                          has to be same size as data array
;
;                 algn_err - provide an estimate for error due to alignment, if given it 
;                            is add to uncertainty in Gaussian fit 
;
;                 fname -  provide the directory to save to and the naming convention, 
;                          i.e. directory/file
;
;                 simp_grad - Use of an analytic least-squares estimate for gradient 
;                             - may be faster!
;
;                 user_func - string name of user defined function. Function must be 
;                             formatted in the same way as mysin.pro (see mpfit tutorial)
;
;                 start - initial guesses for parameter fits. Needs to be defined with 
;                         user_func. Used for initialisation - can edit during fitting.
;
;             used by locate_things
;                 
;                 bad_data - value used in the input data to indicate bad data points. 
;                            When inverting the image, all values <= bad_data will 
;                            be preserved and NOT have their values changed. This 
;                            prevents invalid data flagged with negative intensity 
;                            values from being selected as peaks. Only used when
;                            /invert is also set. Default is -999
; 
;               /despike - if set, will calculate the mean and stddev of the intensity 
;                          values in each timestep. All values with an intensity larger 
;                          than a selected number of stddev above the mean are then 
;                          replaced by the average of the two adjacent data points
;
;               spike_sigma - number of stddev above the mean to use for filtering out 
;                             large data spikes. Only used when /despike is also set. 
;                             Default is 3.0
;
;               num_search_bins - number of bins to look ON EITHER SIDE of a potential 
;                                peak for determining if it is a local maxima. 
;                                Elsewhere known as the 'order' of the search
;                                algorithm. Default is 5. Note, this will also set the
;                                minimum distance allowed between peaks.
;
;               percent_grad - if set to a value between 0 and 100, will automatically 
;                              rescale the intensity gradient to reject the specified 
;                              percent of found peaks.
;               meas_size - sets TOTAL number of pixels used in the Gaussian fit.
;                           Even inputs will have their value increased by one (so as 
;                           to be symmetric about the peak value). Default value is 7.
;
;               /simp_grad - set to use an analytic least-squares estimate for gradient.
;                            This may be faster for very large input arrays!
;
;              cut_chisq - cut off value for unacceptable chi squared. Default set at 3
;               sigma confidence level (based on number of points used in fit)
;              shift_cut - cut off value for maximum allowable peak shift. If the fit
;                          location shifts more than this, will instead use the peak
;                          location with an error of 0.5 pixels. default is 1.5
;              /full_gauss - set to output the 'full_gaussian' results (i.e. save the
;                           widths). Replicates the old output of 'locate_things_fg.pro'
;             /check - if set, will plot various variables about each good peak to a
;                      window. NOT recommended for long periods of observations.
;             /debug - if set, will save sub-pixel test parameters in the 'loc_debug_dat'
;                      COMMON block for later inspection
;
;
;             used by follow threads
;                 min_tlen - minimum thread length to be saved in the output structures
;                 max_dist_jump - maximum allowable distance (in pixels) between two 
;                                 peaks in the same thread. Default is 3 pixels
;                 max_time_skip - maximum allowable timesteps between two peaks in 
;                                 the same thread. Default is 4 time steps.
;                 scan_dir - sets the direction in which the program scans the search 
;                            box. Depending on the input array of peaks, this may 
;                            bias the results.
;                            May choose from the following four options:
;                            "outward" - [DEFAULT] bias towards smaller displacements
;                            "inward" - bias towards larger displacements
;                            "right" - bias towards features that slope to the LEFT
;                            "left" - bias towards features that slope to the RIGHT
;                 /check - if set, will plot each thread after saving it. Should only 
;                          be used for testing the code and NOT general operations.
;                 /debug - if set, will create an array of the same size as td-diagram 
;                          which assigns each pixel a number corresponding to its 
;                          thread (including rejected threads). Saves the array in the 
;                          'th_debug_dat' COMMON block              
;
;CALLS: locate_things.pro, follow_thread.pro, nuwt_moscill.pro
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
;                  R Morton 21/12/2017 Version 2.0 created
;                                      Removed pre-determined models
;                                      Add option for user defined function (user_func)
;                                      Renamed main routine
;                  27 MAR 2018 - Release as version 2.0 of NUWT
;
;
;TO DO: Enable parameter info to be passed to mpfit
;      
;

pro NUWT_TR, data,results=results,errors=errors,$ 
                algn_err=algn_err,meas_size=meas_size,$
                fname=fname,check=check,cut_chisq=cut_chisq,$
                simp_grad=simp_grad,start=start,user_func=user_func,_EXTRA=_extra,no_fill=no_fill

COMMON located_dat,located, nx, nt
COMMON threads_dat, threads
COMMON extended_dat, threads_fit_fg


;Check user function name - be sure it is a scalar string
;Borrowed from mpfit
IF n_elements(user_func) NE 0 THEN BEGIN 
    sz = size(user_func)
    IF sz[0] NE 0 THEN BEGIN
       FCN_NAME:
         errmsg = 'ERROR: user_func must be a scalar string'
        goto, TERMINATE
    ENDIF
    IF sz[sz[0]+1] NE 7 THEN goto, FCN_NAME

    ;Initial call to user defined function
    ;to see if it works
    print,'Checking user supplied function and start values'
    IF n_elements(start) EQ 0 THEN BEGIN
       print,'ERROR: Start variables undefined'
       goto,TERMINATE
    ENDIF  
    fnc=call_function(user_func,findgen(100),start)
ENDIF ELSE print,'Will use sinusoid model for wave fitting'



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
    ;IF not keyword_set(gauss) THEN BEGIN
      ;output is in common located_dat
      IF n_elements(errors) GT 0 THEN inerr=errors[*,*,i] ELSE inerr=0
    ;  locate_things_min,data=data[*,*,i],grad=numg,errors=inerr,simp_grad=simp_grad
    ;  print,'no Gaussian fit used'
    ;ENDIF ELSE BEGIN
         
    ;  locate_things,data=data[*,*,i],grad=numg,meas_size=meas_size, $
    ;             errors=errors[*,*,i],check=check,cut_chisq=cut_chisq,simp_grad=simp_grad
    ;print,invert
    locate_things,data[*,*,i],errors=inerr,res=res,grad=numg, _EXTRA=_extra
    
      IF n_elements(algn_err) GT 0 THEN located.errs=located.errs+algn_err   ; adds alignment error to error from Gaussian fits

    ;  print,'Gaussian fit used'
    ;ENDELSE



    window,0
    !p.multi=[0,2,1]
    tvim,data[*,*,i]
    tvim,located.peaks[*,*,1]
   
    
    ;output is in common thread_dat
    ;for accepted inputs see follow_thread
    follow_threads,min_tlen=numt,_EXTRA=_extra
    szt=size(threads)
    n_thread=szt(1)

    IF n_thread EQ 1 THEN BEGIN
       print,'No threads found'
       goto, TERMINATE
    ENDIF
    print,'Number of threads suitable to fit: '+strtrim(n_thread-1,2)
    
   dsk='y'
   READ,dsk,PROMPT='Do you want to fit the transverse motions?: (y/n) '
   IF dsk EQ 'y' THEN BEGIN

      window,1
      !p.multi=[0,1,1]

      nuwt_moscill,out=threads_fit,start=start,user_func=user_func,no_fill=no_fill
      

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

TERMINATE:

window,0
!p.multi=[0,1,1]
END
