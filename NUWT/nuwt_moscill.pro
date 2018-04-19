;FUNCTION: Fits sinusoidal function plus linear fit ('mysin.pro') to
;          the threads found from 'locate_things.pro', followed by
;          'follow_things.pro'. Interactive program to decide whether
;          to subtract higher order polynomials and whether to save fit.
 

;PROCEDURE OUTLINE: The first step is to remove any examples where
;               the number of pixels in the thread is less than 2 and
;               any threads where more than half the pixels could not
;               be found with locate things (i.e. given the value 0). The routine then locates
;               the start/ends of the threads (defined by -1). Any
;               points that have the value 0 are then set to the value
;               of the previous pixel and given a large error (5
;               pixels). This is to ensure the weighted fitting
;               routines give limited significance to this value. 
;               
;               The routine then asks whether any pixels should be removed
;               from the start/end of the thread. This is followed
;               by asking whether a pre-fit of a high order
;               polynomial (<2) should be performed. Allows for
;               repeat selections so a number of different
;               polynomial fits can be assessed.
;               Routine then fits sinusoidal function.    
;
;INPUTS: data - an array containing outputs from follow_thread,
;               contains thread position (data[*,*,0]) and errors (data[*,*,1])
;        fit_flags - Is a value from 0-3
;                      0 - used for working without full gaussian fit 
;                      1 - code for central position in full gaussian (_FG)
;                      2 - code for maximum intensity in full gaussian (_FG)
;                      3 - code for width in full gaussian (_FG)
;
;        total_flag - Value either 0 or 1
;                     0 - Fitting using only intensity information
;                     1 - Fitting using intensity and Doppler velocities 
;
;
;OPTIONAL: user_func - string name of user defined function. Function must be formatted in the same
;                      way as mysin.pro (see mpfit tutorial)
;          start - initial guesses for parameter fits. Needs to be defined with user_func.
;                   Really just used for initialisation - can edit during fitting process. 
;
;                      
;
;OUTPUTS: out - updated structure with fit details added in thread_fits.fit_result
;               - details of n-parameter model fit to waves, errors are 1-sigma errors
;               fit_result[0:n-1] - parameters
;               fit_result[n:2n-1] - formal 1 sigma errors on parameters
;               fit_result[2n]=chi^2 for fit
;               fit_result[2n+1]=start time from start of time-distance diagram
;               fit_result[2n+2]=end time
;               fit_result[2n+3]=pre-fitting used; (1)-non (2,3....)-
;               polynomial degree
;
;
;HISTORY: created 10/2012 R Morton
;         02/2013 added tdum+findgen(t1dum-tdum) for plot x-values - may cause a bug
;         03/2013 added /damping keyword which allows a damped sine function to be fitted
;         04/2014 updated 'out' array generation - out is now size of number of fitted 
;                 oscillations + added count + added datao
;         R Morton NOV 2014 - super update! Added structure format to remove all the arrays.
;                            Also added COMMON variables so re-used values/structures are passed automatically
;         R Morton 20 APR 2015 - Started a re-structure ready for widgetising!
;         R Morton 21/12/2017 Version 2.0 created
;                             Removed pre-determined models
;                             Add option for user defined function (user_func)
;                  27 MAR 2018 - Release as version 2.0 of NUWT
;
;TO DO: Enable parameter info to be passed to mpfit

FUNCTION find_range,indat
    mom=moment(indat)
    mn=min(indat)
    mx=max(indat)
    return,[mn-2.*sqrt(mom[1]),mx+2.*sqrt(mom[1])]
END

FUNCTION plotthreads,t0,t01,j,range,time,position,fit=fit

     IF n_elements(fit) gt 0 THEN input=position[t0:t01]-fit ELSE input=position[t0:t01]
     plot,time(t0:t01),input,yst=1,xst=1,yrange=[range[0],range[1]],$
          title='Thread No.'+strtrim(j,2)+' between frames t='+strtrim(t0,2)+' and t1='+$
          strtrim(t01,2)
     return,0
END


pro nuwt_moscill,fit_flag=fit_flag,total_flag=total_flag,out=out,start=start,user_func=user_func,no_fill=no_fill

COMMON located_dat,located
COMMON threads_dat,threads

sz=size(located.peaks)
nx=sz(1) & nt=sz(2)

IF n_elements(fit_flag) EQ 0 THEN fit_flag=0
sz=size(threads)
n_thread=sz(1)

IF n_elements(user_func) NE 0 THEN par_len=n_elements(start) ELSE par_len=5

;######################################################################################
;
;SETTING UP OUTPUT STRUCTURES FOR DIFFERENT SCENARIOS
;
;######################################################################################

IF fit_flag EQ 0 THEN BEGIN
;Designed for simple fitting structures
    position=threads.pos ;dummy arrays
    errors=threads.err_pos
    
    ;Create new structure to contain thread details and fit results
    newstruct={pos:fltarr(nt), err:fltarr(nt), fit_result:fltarr(2*par_len+4)}
    threads_fit=replicate(newstruct,n_thread)
    struct_assign,threads,threads_fit

ENDIF ELSE BEGIN

;Designed for advanced fitting structures (_FG)

    ;Define new COMMON variable are threads_fit will now be used multiple times
    COMMON extended_dat, threads_fit_fg

    IF fit_flag EQ 1 THEN BEGIN 
       position=threads.pos ;dummy arrays
       errors=threads.err_pos
    
       ;Create new structure to contain thread details and fit results
       newstruct={pos:fltarr(nt), err_pos:fltarr(nt), inten:fltarr(nt), err_inten:fltarr(nt), $
                  wid:fltarr(nt), err_wid:fltarr(nt), fit_result_pos:fltarr(2*par_len+4), $ 
                  fit_result_inten:fltarr(2*par_len+4), fit_result_wid:fltarr(2*par_len+4) }
       threads_fit_fg=replicate(newstruct,n_thread)
       struct_assign,threads,threads_fit_fg
    ENDIF
    IF fit_flag EQ 2 THEN BEGIN 
       position=threads.inten ;dummy arrays
       errors=threads.err_inten
       
    ENDIF
    IF fit_flag EQ 3 THEN BEGIN 
       position=threads.wid ;dummy arrays
       errors=threads.err_wid
    ENDIF
    
ENDELSE

temp_var=fltarr(2*par_len+4)
time=findgen(nt)

;######################################################################################
;
;MAIN FITTING ROUTINE STARTS
;
;######################################################################################

FOR j=1,(n_thread-1) DO BEGIN
print,'###########################################'
print,'Doing thread '+strtrim(j,2)+' of '+strtrim(n_thread-1,2)
print,'###########################################'

    
    IF NOT keyword_set(no_fill) THEN BEGIN
       ;if value missing set to same as last pixel and set
       zer=where(position[*,j] EQ 0.)
       IF zer[0] NE -1 THEN BEGIN
          FOR ii=0,n_elements(zer)-1 DO BEGIN
              position[zer[ii],j]=position[zer[ii]-1,j]
              errors[zer[ii],j]=1.
          ENDFOR
       ENDIF

    ENDIF

     ;Locate the start and end of thread
     in=where(position[*,j] ne -1.,count)
     t=in[0] & t1=in[-1]

     fitvalues=where(position[*,j] ne -1 and position[*,j] ne 0.)


     range=find_range(reform(position[fitvalues,j]))
     res=plotthreads(t,t1,j,range,time,reform(position[*,j]))
     ;plot,time(t:t1),reform(position[t:t1,j]),yst=1,xst=1,yrange=[range[0],range[1]]
     
     shallwe='n' 
     READ,shallwe,PROMPT='Fit this thread? - y/n '

     IF shallwe EQ 'y' THEN BEGIN
      ;FITTING ROUTINE
      endfit='y'
      ;while loop used to repeat whole fitting procedure
      WHILE (endfit EQ 'y') DO BEGIN
         dummyposition=reform(position[*,j])
         dummyerrors=reform(errors[*,j])
         tdum=t
         t1dum=t1               ;dummy variables for time
         theend=''
         theend='n'
         talk='new'
         fitvalues=where(dummyposition ne -1 and dummyposition ne 0.)
         
         ;while loop used to repeat polynomial fitting
         WHILE (theend EQ 'n') DO BEGIN
            
            range=find_range(dummyposition[fitvalues])
            ;plot,time(t:t1),dummyposition[tdum:t1dum],yst=1,xst=1,yrange=[range[0],range[1]]
	          res=plotthreads(tdum,t1dum,j,range,time,dummyposition)

            ;++++++++++++++++++++++++++++++++++++++++++++++
            ;'Does the length of the thread need changing?' loop
            IF talk EQ 'new' THEN BEGIN
               READ,talk,PROMPT='Change start end of thread? - y/n '

               IF talk EQ 'y' THEN BEGIN
                  print,'Current value of t='+strtrim(t,2)+' and t1='+strtrim(t1,2)
                  READ,tdum,PROMPT='Enter t value: '
                  READ,t1dum,PROMPT='Enter t1 value: '
               ENDIF
               
               fitvalues=fitvalues[where(fitvalues ge tdum and fitvalues le t1dum)]

               range=find_range(dummyposition[fitvalues])
              ; plot,time(tdum:t1dum),dummyposition(tdum:t1dum),yst=1,xst=1,yrange=[range[0],range[1]]
	             res=plotthreads(tdum,t1dum,j,range,time,dummyposition)
            ENDIF
            ;++++++++++++++++++++++++++++++++++++++++++++++


            ;Deciding whether a cubic or quadratic needs to fitted and subtracted
            ;before being sent to MPFITFUN
            talk3=1
            talk32='n'

            READ,talk32,PROMPT='Change polynomial degree? - y/n '

            IF talk32 EQ 'n' THEN theend='y'
       
            IF talk32 EQ 'y' THEN BEGIN
               READ,talk3,PROMPT='Fit trend - enter polynomial degree: '

               IF talk3 EQ 1 THEN BEGIN 
                  theend='y'
               ENDIF ELSE BEGIN

                 res=poly_fit(fitvalues,dummyposition[fitvalues],talk3,$
                 measure_errors=dummyerrors[fitvalues],yfit=fit,yband=yband,chisq=chisq)
                
                 range=find_range(dummyposition[fitvalues]-fit)
                 ;plot,time(tdum:t1dum),dummyposition[tdum:t1dum]-fit,yrange=[range[0],range[1]]
                 res=plotthreads(tdum,t1dum,j,range,time,dummyposition,fit=fit)
                 print,'Chi^2 of fit ',chisq

                 READ,theend , PROMPT='Is the fit good? y/n '

                 IF theend EQ 'y' THEN BEGIN
                    dummyposition[fitvalues]=dummyposition[fitvalues]-fit
                    dummyerrors[fitvalues]=sqrt(dummyerrors[fitvalues]^2+yband^2)
                 ENDIF

              ENDELSE
            ENDIF

          

         ENDWHILE
      
      ;++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
      ;++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      ;non-linear fit to thread
      ;perror is 1-sigma on parameters
      ;bestnorm is summed squared residuals (chi^2)
      
      ;Default fitting - a sinusoid
      IF n_elements(user_func) EQ 0 THEN BEGIN
          start=[dummyposition[tdum],1.,20.,0.5,1.]
          print,'Initial variables:',' Constant='+strtrim(start[0],2)+' Amplitude='+strtrim(start[1],2)+$
              ' Period='+strtrim(start[2],2)+' Phase='+strtrim(start[3],2)+' Linear='+strtrim(start[4],2)
          st=''
          READ,st,PROMPT='Change initial variables? y/n '

          IF st EQ 'y' THEN BEGIN
            new1=start[0]          ;READ,new1,PROMPT='Enter constant: '
            READ,new2,PROMPT='Enter amplitude: '
            READ,new3,PROMPT='Enter period: '
            new4=start[3]          ;READ,new4,PROMPT='Enter phase: '
            new5=start[4]          ;READ,new5,PROMPT='Enter linear: '
            start=[new1,new2,new3,new4,new5]
          ENDIF
        
          xplot=time[fitvalues] ;Set length for fitting/plot 
    
          res=mpfitfun('mysin',xplot,dummyposition[fitvalues],$
                     dummyerrors[fitvalues],start,perror=perror,bestnorm=bestnorm,/quiet)

      ENDIF ELSE BEGIN
          ;For use defined fitting function
          start[0]=dummyposition[tdum]
          print,'Initial variables:', start
          st=''
          READ,st,PROMPT='Change initial variables? y/n '

          IF st EQ 'y' THEN BEGIN
            new=fltarr(par_len)          
            READ,new,PROMPT='Enter '+strtrim(par_len,2)+' new start parameters (comma separated): '
            start=new
          ENDIF
        
          xplot=time[fitvalues];Set length for fitting/plot 
    
          res=mpfitfun(user_func,xplot,dummyposition[fitvalues],$
                     dummyerrors[fitvalues],start,perror=perror,bestnorm=bestnorm,/quiet)

      ENDELSE     

      print,'%'
      print,'%'
      print,'Fitted variables',res
      print,'Error on fits',perror
      print,'Chi^2', bestnorm

      oploterr,time[fitvalues],dummyposition[fitvalues],dummyerrors[fitvalues],psym=1
      
      IF n_elements(user_func) EQ 0 THEN BEGIN 
        pl_func=mysin(xplot,res)
        oplot,time[fitvalues],pl_func,linestyle=2
      ENDIF ELSE BEGIN
        pl_func=call_function(user_func,xplot,res)
        oplot,time[fitvalues],pl_func,linestyle=2
      ENDELSE
      

      READ,endfit,PROMPT='Repeat fitting procedure? - y/n ' 
     ENDWHILE


     sav=''
     READ,sav,PROMPT='Save details? -  y/n '

     IF sav EQ 'y' THEN BEGIN

     	temp_var[0:par_len-1]=res
      temp_var[par_len:2*par_len-1]=perror
     	temp_var[2*par_len]=bestnorm  ; chi^2
     	temp_var[2*par_len+1]=tdum      ;Start time
     	temp_var[2*par_len+2]=t1dum     ;end time
     	temp_var[2*par_len+3]=talk3     ; degree of polynomial fit used
     
     	
     	;######################################################################################
     	;
     	;SAVING OUTPUT STRUCTURES FOR DIFFERENT SCENARIOS
     	;
     	;######################################################################################

	IF fit_flag EQ 0 THEN BEGIN
        	threads_fit[j].fit_result=temp_var
        	threads_fit[j].pos[*]=position[*,j]
        	threads_fit[j].err[*]=errors[*,j]
    	ENDIF ELSE BEGIN
          IF fit_flag EQ 1 THEN BEGIN
             threads_fit_fg[j].fit_result_pos=temp_var
             threads_fit_fg[j].pos[*]=position[*,j]
             threads_fit_fg[j].err_pos[*]=errors[*,j]
          ENDIF
          IF fit_flag EQ 2 THEN BEGIN
             threads_fit_fg[j].fit_result_inten=temp_var
             threads_fit_fg[j].inten[*]=position[*,j]
             threads_fit_fg[j].err_inten[*]=errors[*,j]
          ENDIF
          IF fit_flag EQ 3 THEN BEGIN
             threads_fit_fg[j].fit_result_wid=temp_var
             threads_fit_fg[j].wid[*]=position[*,j]
             threads_fit_fg[j].err_wid[*]=errors[*,j]
          ENDIF
     	 ENDELSE
    ENDIF ;saving endif   

  ENDIF ;shallwe fit endif

  
 ;ENDIF
 ;ENDIF

ENDFOR

IF fit_flag EQ 0 THEN out=threads_fit ELSE out=threads_fit_fg


END
