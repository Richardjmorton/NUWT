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
;OPTIONAL: damped - set if damped sin function is wanted, otherwise constant sin function is fitted
;
;                      
;
;OUTPUTS: out - updated structure with fit details added in thread_fits.fit_result
;               - details of sinusoidal fit to waves, errors are 1-sigma errors
;               fit_result[0]=constant
;               fit_result[1]=amplitude
;               fit_result[2]=period
;               fit_result[3]=phase
;               fit_result[4]=linear coefficient
;               fit_result[5]=error on constant
;               fit_result[6]=error on amplitude
;               fit_result[7]=error on period
;               fit_result[8]=error on phase
;               fit_result[9]=error on linear coef
;               fit_result[10]=chi^2 for fit
;               fit_result[11]=start time from start of time-distance diagram
;               fit_result[12]=end time
;               fit_result[13]=pre-fitting used; (1)-non (2,3....)-
;               polynomial degree
;               fit_result[14]=damping term (if /damped used)
;               fit_result[15]=error on damping term
;	         fit_result[16]=time-dependent period term (if /tdp used)
;               fit_result[17]=error on time-dependent period term
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
;

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


pro moscill,damped=damped,tdp=tdp,fit_flag=fit_flag,total_flag=total_flag,out=out

COMMON located_dat,located,nx,nt
COMMON threads_dat,threads

IF n_elements(fit_flag) EQ 0 THEN fit_flag=0
sz=size(threads)
n_thread=sz(1)

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
    newstruct={pos:fltarr(nt), err:fltarr(nt), fit_result:fltarr(18)}
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
                  wid:fltarr(nt), err_wid:fltarr(nt), fit_result_pos:fltarr(18), fit_result_inten:fltarr(18), fit_result_wid:fltarr(18) }
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

temp_var=fltarr(18)
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

   
  ;THE FOLLOWING TWO IF STATMENTS HAVE BEEN MOVED TO FOLLOW_THREAD.PRO
   ;skips enteries with less than 2 positive values
 ;  IF (n_elements(where(position[*,j] gt 0.))) GT 2. THEN BEGIN

   ;skips enteries where half the data points are set to zero, i.e. no
   ;value was obtained at fitting stage.
 ;  IF (n_elements(where(position[*,j] EQ 0.))) LT 0.5*(n_elements(where(position[*,j] GE 0.))) THEN BEGIN

     ;if value missing set to same as last pixel and set
     zer=where(position[*,j] EQ 0.)
     IF zer[0] NE -1 THEN BEGIN
        FOR ii=0,n_elements(zer)-1 DO BEGIN
            position[zer[ii],j]=position[zer[ii]-1,j]
            errors[zer[ii],j]=1.
        ENDFOR
     ENDIF

     ;Locate the start and end of thread
     in=where(position[*,j] ne -1.,count)
     chac=n_elements(in)
     t=in[0] & t1=in[chac-1]


     range=find_range(reform(position[t:t1,j]))
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
         
         
         ;while loop used to repeat polynomial fitting
         WHILE (theend EQ 'n') DO BEGIN
            
            range=find_range(dummyposition[tdum:t1dum])
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
               
               range=find_range(dummyposition[tdum:t1dum])
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

                 res=poly_fit(findgen(t1dum-tdum+1),dummyposition[tdum:t1dum],talk3,$
                 measure_errors=dummyerrors[tdum:t1dum],yfit=fit,yband=yband,chisq=chisq)
                
                 range=find_range(dummyposition[tdum:t1dum]-fit)
                 ;plot,time(tdum:t1dum),dummyposition[tdum:t1dum]-fit,yrange=[range[0],range[1]]
                 res=plotthreads(tdum,t1dum,j,range,time,dummyposition,fit=fit)
                 print,'Chi^2 of fit ',chisq

                 READ,theend , PROMPT='Is the fit good? y/n '

                 IF theend EQ 'y' THEN BEGIN
                    dummyposition[tdum:t1dum]=dummyposition[tdum:t1dum]-fit
                    dummyerrors[tdum:t1dum]=sqrt(dummyerrors[tdum:t1dum]^2+yband^2)
                 ENDIF

              ENDELSE
            ENDIF

          

         ENDWHILE
      
      ;++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
      ;++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      ;non-linear fit to thread
      ;perror is 1-sigma on parameters
      ;bestnorm is summed squared residuals (chi^2)
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
      
      xplot=findgen(t1dum-tdum+1) ;Set length for fitting/plot 
  
      IF keyword_set(damped) THEN BEGIN
         IF keyword_set(tdp) THEN BEGIN
            start2=fltarr(7)
            start2[0:4]=start
            dmp=1d          
            READ,dmp,PROMPT='Enter damping term: '
            start2[5]=dmp
            tdp2=1d          
            READ,tdp2,PROMPT='Enter time-dependent term: '
            start2[6]=tdp2
       
             res=mpfitfun('mydampedsin_tdp',xplot,dummyposition[tdum:t1dum],$
                   dummyerrors[tdum:t1dum],start2,perror=perror,bestnorm=bestnorm,/quiet)
         
         ENDIF ELSE BEGIN
           start2=fltarr(6)
           start2[0:4]=start
           dmp=1d          
           READ,dmp,PROMPT='Enter damping term: '
           start2[5]=dmp

           res=mpfitfun('mydampedsin',xplot,dummyposition[tdum:t1dum],$
                   dummyerrors[tdum:t1dum],start2,perror=perror,bestnorm=bestnorm,/quiet)

          ENDELSE

      ENDIF ELSE BEGIN
         
         res=mpfitfun('mysin',xplot,dummyposition[tdum:t1dum],$
                   dummyerrors[tdum:t1dum],start,perror=perror,bestnorm=bestnorm,/quiet)
      ENDELSE     

      print,'%'
      print,'%'
      print,'Fitted variables',res
      print,'Error on fits',perror
      print,'Chi^2', bestnorm

      oploterr,time(tdum:t1dum),dummyposition[tdum:t1dum],dummyerrors[tdum:t1dum]
      
      
      IF keyword_set(damped) THEN BEGIN
         IF keyword_set(tdp) THEN BEGIN
            oplot,time(tdum:t1dum),res[0]+res[1]*sin(2.*!pi*xplot/(res[2]*(1+res[6]*xplot))-res[3])*exp(-res[5]*xplot)+xplot*res[4],linestyle=2

         ENDIF ELSE BEGIN
            oplot,time(tdum:t1dum),res[0]+res[1]*sin(2.*!pi*xplot/res[2]-res[3])*exp(-res[5]*xplot)+xplot*res[4],linestyle=2
         ENDELSE
      ENDIF ELSE BEGIN 
         oplot,time(tdum:t1dum),res[0]+res[1]*sin(2.*!pi*xplot/res[2]-res[3])+xplot*res[4],linestyle=2
      ENDELSE 

      READ,endfit,PROMPT='Repeat fitting procedure? - y/n ' 
     ENDWHILE


     sav=''
     READ,sav,PROMPT='Save details? -  y/n '

     IF sav EQ 'y' THEN BEGIN

     	temp_var[0]=res[0] ;constant
     	temp_var[1]=res[1] ;amplitude
     	temp_var[2]=res[2] ;period
     	temp_var[3]=res[3] ;phase
     	temp_var[4]=res[4] ; out[4,j]=res[4] ;linear coefficient
     	temp_var[5]=perror[0]  ;constant error
     	temp_var[6]=perror[1]  ;amplitude error
     	temp_var[7]=perror[2]  ;period error
     	temp_var[8]=perror[3]  ;phase error
     	temp_var[9]=perror[4]  ;linear error
     	temp_var[10]=bestnorm  ; chi^2
     	temp_var[11]=tdum      ;Start time
     	temp_var[12]=t1dum     ;end time
     	temp_var[13]=talk3     ; degree of polynomial fit used
     
     	IF keyword_set(damped) THEN BEGIN
        	temp_var[14]=res[5]    ;damping time
        	temp_var[15]=perror[5] ;damping time error
        
      		IF keyword_set(tdp) THEN BEGIN
             		temp_var[16]=res[6] ;value of time-dependence
             		temp_var[17]=perror[6] ;error on time-dependence
          	ENDIF
     	ENDIF


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
