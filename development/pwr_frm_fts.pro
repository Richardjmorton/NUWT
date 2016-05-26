
;PURPOSE - Takes a time distance diagram, picks out oscillating threads, then feeds them into an FFT to calculate the Power Spectral 
;          Density (PSD). Calculates the power for each individual thread and also the power averaged over all threads. Average
;          is calculated by creating frequency bins, taking the ln of power in frequency bins, plot the ln power PDF and fitting a
;          Gaussian. Centroid of Gaussian (mu) gives median value of distribution, i.e. exp(mu)=median.  
;
;
;INPUTS - data - an array of time-distance (t-d) diagrams of the form (x,t,z)
;
;OPTIONAL INPUTS - dx - spatial sampling of data. default is 1
;                - dt - temporal sampling of the data. default is 1
;                - signif_levels - sets significance levels for FFT results (default=0.95) 
;                - doplot - plots of threads and power spectra on screen for debug
;                
;                - gauss - uses a gaussian fit to locate centre of thread, provides sub-pixel measurements
;                - errors - estimated errors on intensity values, supplied to gaussian fitting routine.
;                           if not supplied a default value of 1% of intensity values is used
;                - if keyword set will remove linear trend from each thread. If not set mean value is removed
;                - smpad - SMooth And PAD - smooths the time-series using 3 pixel box-car average and
;                          pads with zeros to increase frequency resolution
;                - TNOG - Thread Number Guess, guess hows many threads will be found default is 100
;                         At present, once limit is hit routine will just move onto the next slit. 
;
;OUTPUTS - power - PSD results for each thread in the form (PSD, thread# , t-d#)
;          freq - corresponding frequencies for power array in the form (freq, thread#, t-d#)

;
;
;EXTERNAL CALLS - locate_things_min.pro, locate_things.pro, follow_thread.pro, cospd.pro
;
; 
;TO DO - Make histogram plot generic -i.e. no reference to directory
;      - Add method to generate arrays related to size of results rather than specifying at start
;        Will probably need to re-write routine structure so that save file is generated for each slit
;        and arrays are extended with each result.
;        Some redundancy in the routine - variable h counts threads, while the number is already known!
;      - TIDY UP!
;
;HISTORY: Created by R Morton 2014
;         R Morton MAR2015 - Updated to include structures 
;
;
;RESTRICTIONS: LOTS AT THE MOMENT. ROUTINE IS VERY MUCH DEPENDENT UPON USERS
;


; define temporal apodization
FUNCTION do_apod,tn,cpg
 
      apodt = fltarr(tn)+1
      apod=0.2
      apodrimt = tn*apod
      apodt[0] = (sin(!pi/2.*findgen(apodrimt)/apodrimt))^2
      apodt = apodt*shift(rotate(apodt,2),1) 
      CPG=total(apodt)/n_elements(apodt) ;Coherent Power Gain

RETURN,apodt
END

PRO pwr_frm_fts,data,dx=dx,dt=dt,doplot=doplot, gauss=gauss,errors=errors,ctrend=ctrend,$
                     fname=fname,signif_levels=signif_levels,smpad=smpad,tnog=tnog



file='data/SOT' ;'/users/richardmorton/data/ROSA/sep2010_res/velocity_maps'

;###############################################
;INITIAL HOUSE KEEPING
;###############################################

COMMON located_dat,located, nx, nt

sz=size(data)
IF sz(0) EQ 3 THEN nslits=sz(3) ELSE nslits=1

IF NOT keyword_set(smpad) THEN nf=sz(2)/2+1 ELSE nf=(sz(2)+1000)/2+1
IF n_elements(TNOG) EQ 0 THEN TNOG=100


power=fltarr(nf,TNOG,nslits)
displacement=fltarr(nf,TNOG,nslits)
freq=fltarr(nf,TNOG,nslits)
duration_series=fltarr(nf,TNOG,nslits)
phase2=fltarr(nf,TNOG,nslits)

displacement_max=[0.]
freq_max=[0.]
series_duration=[0.]

IF NOT keyword_set(dx) THEN dx=1.
IF NOT keyword_set(dt) THEN dt=1.
IF NOT keyword_set(signif_levels) THEN signif_levels=0.95

IF NOT keyword_set(gauss) THEN print, 'Finding whole pixel values' ELSE print, 'Gaussian location enabled'

numg=1d
READ,numg,PROMPT='Enter gradient: '

numt=1d
READ,numt,PROMPT='Enter thread length: '


;###############################################
;BEGIN WORKING ON EACH TIME-DISTANCE DIAGRAM
;###############################################

FOR k=0, nslits-1 DO BEGIN

    h=0 ; Counting device
    
    print,'&&&&&&&','Slit number:', k
      
    ;Fitting of intensity maxima
    IF not keyword_set(gauss) THEN BEGIN
      locate_things_min_test,data=data[*,*,k],grad=numg
    ENDIF ELSE BEGIN

         IF n_elements(errors) gt 0. THEN errorsi=errors[*,*,k]
 
         ;If errors are not supplied then error level set to default 10% of intensity
         IF n_elements(errors) eq 0 THEN BEGIN
            errorsi=fltarr(sz(1),sz(2))
            errorsi=abs(data[*,*,k])*0.1             
         ENDIF
  
         locate_things_fg_test,data=data[*,*,k],grad=numg,meas_size=meas_size, $
               errors=errorsi
    ENDELSE

    ;Plot found data points
    IF keyword_set(doplot) THEN BEGIN
          window,0
          !p.multi=[0,2,1]
          tvim,data[*,*,k]
          tvim,located.peaks[*,*,1]
    ENDIF
    
      
    follow_thread_test,min_tlen=numt,area=area
    COMMON threads_dat, threads

    sz=size(threads)
    n_threads=sz(1)
    print,'Number of threads',n_threads

    ;###############################################
    ;NOW CANDIDATE THREADS HAVE BEEN FOUND, SORT THEN FIT
    ;###############################################
    FOR j=0,(n_threads-1) DO BEGIN

            dummy=threads[j].pos ;dummy array
            terr=threads[j].err_pos

        ;BEGIN SORTING THREADS - THROW AWAY BAD ONES
       	 t=-1000. & t1=-1000. ;reset each time 

   	   ;skips enteries with less than 2 positive values
           IF (n_elements(where(dummy gt 0.))) GT 2. THEN BEGIN

               ;skips enteries where half the data points are set to zero, i.e. no
               ;value was obtained at fitting stage.
               IF (n_elements(where(dummy EQ 0.))) LT 0.5*(n_elements(where(dummy GE 0.))) THEN BEGIN

                   ;if value missing set to same as last pixel and set
                   zer=where(dummy EQ 0.)
                   IF zer[0] NE -1 THEN BEGIN
                     FOR ii=0,n_elements(zer)-1 DO BEGIN
                         dummy[zer[ii]]=dummy[zer[ii]-1]
                         terr[zer[ii]]=1.
                     ENDFOR
                   ENDIF

                   ;Locate the start and end of thread
                   in=where(dummy ne -1.)
                   chac=n_elements(in)
                   t=in[0] & t1=in[chac-1]


                   ;###############################################
                   ;FIT THE GOOD THREADS
                   ;###############################################
                   IF t gt 0 THEN BEGIN

                      dummy=dx*dummy[t:t1]
                                            
                      ; calculates linear trend or mean
                      IF keyword_set(ctrend) THEN res=poly_fit(findgen(t1-t+1),dummy,1,yfit=trend) $ 
                      ELSE trend=mean(dummy)  ;calculates mean of data
          
            

                      ;###############################################
                      ;Calculate power and calculates significance level of the signal
                      ; (see, e.g., Vaughan 2005)
                      ;###############################################
            
                      oscil=dummy-trend
       
           		s_len=n_elements(oscil)
           		sd=s_len*dt
           		apodt=do_apod(s_len,cpg) ; cpg is Coherent Power Gain - correction factor needed for power after apodisation
           		

           		IF NOT keyword_set(smpad) THEN BEGIN
                 
                  		oscil=oscil*apodt
                  		n_len=n_elements(oscil)  

                 		df=1./(s_len*dt)
                  		f=findgen(s_len/2)*df
                  		f=f(1:(s_len-1)/2-1)
                  		pow=(2.*(abs(fft(oscil,-1)))^2)[1:(s_len-1)/2-1]
                  		phase=atan((fft(oscil,-1))[1:(s_len-1)/2-1],/phase)

           		ENDIF ELSE BEGIN
                  
                  		oscil=smooth(oscil,3,/edge_truncate)
                  		oscil=[oscil*apodt,fltarr(1000)]
                  		n_len=n_elements(oscil)
                  
                  		df=1./(n_len*dt)
                  		f=findgen(n_len/2)*df
                  		f=f(1:(n_len-1)/2-1)
                  		pow=(2.*(abs(fft(oscil,-1)))^2)[1:(n_len-1)/2-1]
                  		phase=atan((fft(oscil,-1))[1:(n_len-1)/2-1],/phase)
                  
           		ENDELSE
           

            		;Doesn't let straight lines through
           		IF ((moment(pow))[1]) GT 1e-9 THEN BEGIN
            
             			;Significant tests for power spectra
             			;Default is Torrence & Compo 1998
             			;Alternative is Vaugan 2005
             			IF keyword_set(vaughn) THEN BEGIN
               				 ;Normalise power
             				  npw=s_len*pow/(moment(oscil))[1] 
                			  prob=1.-chisqr_pdf(2.*npw,2)
                			  nprob = MAKE_ARRAY(SIZE(prob, /DIM))
                			  log_pp = s_len * ALOG(DOUBLE(1.0 - prob))
                			  indx = WHERE(log_pp gt -30.0, count, COMPLEMENT=indx_c)
                			  IF (count gt 0) THEN nprob[indx] = exp(log_pp[indx])
                			  IF (count lt s_len) THEN nprob[indx_c] = 0.0
  
                			  sig_vals=where(nprob gt signif_levels)
  
             			ENDIF ELSE BEGIN

                			sig_lvl = SIGNIF_CONF(oscil,signif_levels)
                			sig_vals=where(pow gt sig_lvl)
             	
            			ENDELSE
           
             			power(sig_vals,h,k)=pow[sig_vals]*(n_len/s_len/CPG)^2 ;Right correction for power?
             			consec,sig_vals,lo,hi,num

             			nsv=n_elements(sig_vals)
             			totcons=total(hi-lo)+num
                         

             			IF (nsv GT 0) AND (sig_vals[0] ne -1) THEN BEGIN

                			IF nsv EQ totcons THEN BEGIN
                   				FOR ii=0,num-1 DO BEGIN
                      					 bran=sig_vals[lo[ii]] & tran=sig_vals[hi[ii]]
                      					 m=max(pow[bran:tran],loc)
                    
                      					 freq(ii,h,k)=f[bran+loc]
                       				 displacement(ii,h,k)=2.*sqrt(pow[bran+loc]/2.)*n_len/s_len/CPG
                       				 phase2[ii,h,k]=phase[bran+loc]
                       				 duration_series[ii,h,k]=sd
                       				 IF keyword_set(doplot) THEN print,displacement(ii,h,k),1./f[bran+loc]
                       				 IF displacement[ii,h,k] GT 2000 THEN print,displacement[ii,h,k],1./freq[ii,h,k]
                   				ENDFOR
                			ENDIF ELSE BEGIN
                    				IF totcons EQ 0 THEN BEGIN
                      					displacement[0:nsv-1,h,k]=2.*sqrt(pow[sig_vals]/2.)*n_len/s_len/CPG
                       				freq[0:nsv-1,h,k]=f[sig_vals]
                       				phase2[0:nsv-1,h,k]=phase[sig_vals]
                       				duration_series[0:nsv-1,h,k]=sd
                       
                   				ENDIF ELSE BEGIN
                       				FOR ii=0,num-1 DO BEGIN
                           					bran=sig_vals[lo[ii]] & tran=sig_vals[hi[ii]]
                           					m=max(pow[bran:tran],loc)
                    
                           				        freq(ii,h,k)=f[bran+loc]
                           					displacement(ii,h,k)=2.*sqrt(pow[bran+loc]/2.)*n_len/s_len/CPG
                           					phase2[ii,h,k]=phase[bran+loc]
                           					duration_series[ii,h,k]=sd
                           					;print,displacement(ii,h,k),1./f[bran+loc]

                           					sig_vals[lo[ii]:hi[ii]]=0.
                   
                       				ENDFOR
                       				nsig_vals=sig_vals[where(sig_vals gt 0.)]
                       				nnsv=n_elements(nsig_vals)                       
                       				displacement[num:num+nnsv-1,h,k]=2.*sqrt(pow[nsig_vals]/2.)*n_len/s_len/CPG
                       				freq[num:num+nnsv-1,h,k]=f[nsig_vals]
                       				phase2[num:num+nnsv-1,h,k]=phase[nsig_vals]
                       				duration_series[num:num+nnsv-1,h,k]=sd
                    				ENDELSE  

                			ENDELSE
             			ENDIF
             
           
             			;Calculate for maximum value   
             			m=max(2.*sqrt(pow/2.)*n_len/s_len/CPG,loc)
             			displacement_max=[displacement_max,m]          
             			freq_max=[freq_max,f(loc)]
             			series_duration=[series_duration,sd]

             			h=h+1

             			;Cheap way to stop routine crashing if number of threads found is
             			;larger than pre-determined array.
             			IF h EQ TNOG THEN BEGIN
               				 PRINT,'Capacity for threads met - may need to increase TNOG'
               				 goto,end_loop
             			ENDIF

             			;For debugging
             			IF keyword_set(doplot) THEN BEGIN
                  			 window,1
                   			!p.multi=[0,3,1]
                   			plot,findgen(s_len)*dt,(dummy-trend)*apodt,charsize=2.
                   			oplot,findgen(s_len)*dt,(dummy-trend),linestyle=2
                   
                   			fun=displacement[0,h-1,k]*cos(2.*!Pi*freq[0,h-1,k]*findgen(s_len)*dt+phase2[0,h-1,k])
                   			oplot,findgen(s_len)*dt,fun,linestyle=3,thick=2
		
                   			IF num GT 1 THEN BEGIN
                     				 fun2=displacement[1,h-1,k]*cos(2.*!Pi*freq[1,h-1,k]*findgen(s_len)*dt+phase2[1,h-1,k])
                      				 oplot,findgen(s_len)*dt,fun2,linestyle=4,thick=2
                   			ENDIF
                   			IF num GT 2 THEN BEGIN 
                    				fun3=displacement[2,h-1,k]*cos(2.*!Pi*freq[2,h-1,k]*findgen(s_len)*dt+phase2[2,h-1,k])
                   				 oplot,findgen(s_len)*dt,fun3,linestyle=5,thick=2
                  			ENDIF

                  			plot,f,pow,charsize=2.
                   			pow2=pow
                   			IF keyword_set(vaughn) THEN pow2[where(nprob lt 0.95)]=0. $
                   			ELSE pow2[where(pow lt sig_lvl)]=0.
                   			plot,f,pow2,charsize=2.
                   			pause
                                      !p.multi=0

             		        ENDIF

           		ENDIF
          
          	    ENDIF

               ENDIF
            ENDIF
        ENDFOR

           

;Routine jumps to here is array gets full!
end_loop:
print,'Number of features found and measured: ',h-1

ENDFOR   

savegenx,power,displacement,freq,duration_series,displacement_max,freq_max,series_duration,dt,dx,file=file+'/'+fname,/overwrite
!p.multi=0
set_plot,'x'
END

