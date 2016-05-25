;FUNCTION: Fits full Gaussian around an intensity peak in time-distance diagram 
;          
;
;PROCEDURE OUTLINE: Takes a (high-pass filtered)
;time-distance diagram and locates the peaks using a find_max crawling
;routine.The maximum then has to have a specific gradient (default
;>0.5) to be classed as a peak. Once a maximum is found a Gaussian fit
;is applied to the surrounding 9 pixels (four either side) to provide
;subpixel position for the maximum. Values for central intensity and structure
;width (Gaussian standard deviation) are also recorded.
;
;INPUTS: data - single time-distance diagram
;        grad - limits on gradient (default = 0.5)
;
;OPTIONAL INPUTS: meas_size - sets number of pixels used in Gaussian fit, 
;                             default value is 9, odd numbers only
;                 errors - errors on intensity for each pixel, i.e., estimates for Photon 
;                          noise, array should be same size as data array. Supplied to 
;                          Gaussian fitting routine. If not set default is 0 errors which will
;                          give meaningless errors on position, intensity, width measurements.
;
;                 shift_cut - Controls cut off value for Guassian fit, default 1.5. I.e. if gaussian peak position is
;                             more than 1.5 pixels from maximum intensity position the gaussian fit is ignored.
;
;                 cut_chisq - Controls cut off value for reduced chi squared, default 3 sigma from chi^2 
;                             distribution with appropriate degress of freedom. I.e. if reduced chi squared of the fit
;                             is greater than 5 then gaussian fit is ignored.
;
;                 check - 
;                 
;                 simp_grad - Use of an analytic least-sqaures estimate for gradient - may be faster!            
;
;OUTPUTS: located - structure containing both the positions and errors of found peaks.
;                   Saved as COMMON variable so not a routine output as such! 
;
;TO DO: Still need to work out the best way to describe intensity values for gaussian fit!
;
;HISTORY:
; Created- R. Morton Oct - 2014 - based on locate_things.pro
;          R Morton NOV 2014 - super update! Added structure format to remove all the arrays.
;                            Also added COMMON variables so re-used values/structures are passed automatically
;          R Morton MAR 2015 - fixed bug with definition of initial minimum value to set all located.peaks values to
;          R Morton APR 2015 - made the check keyword plot nice & informative
;          14 MAR 2016 - Release as version 1.0 of NUWT 
;          R Morton MAY 2016 - Added option for analytic LS calc for graident in hope of speed increase
;                              Also removed needless for loop at bottom!!
;                              Added calcuation of peak intensity
;                              Added better estimate for Chi^2 cutoff 
;
;CALLS: Requires the following routines to function: mpfitfun, fifteenb, 
;

pro locate_things_fg,data=data,grad=grad,meas_size=meas_size, $
                    errors=errors,check=check,cut_chisq=cut_chisq,shift_cut=shift_cut,simp_grad=simp_grad    

COMMON located_dat,located, nx, nt

;Default for measuring fibrils in ROSA data
IF NOT KEYWORD_SET(meas_size) THEN meas_size=7 
meas_size=fix(meas_size)
len=meas_size/2

;Set cut off value for unacceptable chi squared
;Set at value of 3 sigma confidence level
IF NOT KEYWORD_SET(cut_chisq) THEN cut_chisq=chisqr_cvf(1.-0.9973,meas_size-5) 

;Set cut off value for unacceptable peak shift
IF NOT KEYWORD_SET(shift_cut) THEN shift_cut=1.5 

;Sets default gradient if non specified
IF n_elements(grad) EQ 0 THEN grad=0.5

;Set up structures to store
sz=size(data)
nx=sz(1) & nt=sz(2)
located={peaks:fltarr(nx,nt,6), errs:fltarr(nx,nt,6)}
located.peaks(0:nx-1,0:nt-1,0:5)=(min(data)-10.)<(-10) ;Should ensure no confusion
;maxi={vals:fltarr(nx,6), errs:fltarr(nx,6)}

mini=min(located.peaks)



;Begin search in each time-slice
in=[-2,-1,0,1,2]
in2=-(meas_size/2)+indgen(meas_size)
IF len lt 5 THEN start_val=5 ELSE start_val=len
FOR j=0, (nt-1) DO BEGIN
       
       ;h=0
       image=reform(data[0:nx-1,j,0])
       err_dat=reform(errors[0:nx-1,j,0])

       ;Search each along x-direction, excludes edges
       FOR i=start_val, (nx-start_val-1) DO BEGIN
          
          ;Finds maximum
          mx=max(image[i-5:i+5],loc)
       
          ;Wait till maximum is at centre of search bar
          IF loc EQ 5 THEN BEGIN

             ;Finds gradients either side of the maximum
             IF NOT keyword_set(simp_grad) THEN BEGIN
                 res=poly_fit(in,image[i-4:i],1,yfit=fit,measure_errors=err_dat[i-4:i])
                 res2=poly_fit(in,image[i:i+4],1,yfit=fit2,measure_errors=err_dat[i:i+4])
                 m1=res[1] & m2=res2[1]
             ENDIF ELSE BEGIN
                 m1=total(in*(image[i-4:i]-total(image[i-4:i])/5.))/10. ;Analytic gradient
                 m2=total(in*(image[i:i+4]-total(image[i:i+4])/5.))/10.
             ENDELSE    

             ;If gradients greater than a certain value begin
             ;Gradient of quadratic evaluated at x=0
             IF (m1 GT grad) AND (m2 LT (-1.)*grad) THEN BEGIN

                ;find gaussian fit to surrounding points
                estimates=[image[i],0.,2.,min(image[(i-len):(i+len)]),0.1]
                coeff=mpfitfun('mygauss_plus_linear',in2,image[(i-len):(i+len)],err_dat[(i-len):(i+len)],estimates,$
                             perror=sigma,bestnorm=bestnorm,dof=dof,/quiet)
                chisq=bestnorm
                chisq_red=bestnorm/dof
                
               ;Plots useful stuff to a window
               IF keyword_set(check) THEN BEGIN
                  x=i-len+findgen(2*len+1)
                  toplo=image[(i-len):(i+len)]
                  erplo=err_dat[(i-len):(i+len)]
                  yran=[min(toplo)-0.01*abs(mean(toplo)),max(toplo)+0.01*abs(mean(toplo))]
                  xran=[min(x)-1,max(x)+1]
                  plot,x,toplo,thick=3,yst=1,title='Time frame '+strtrim(j,2),yrange=yran,xrange=xran,xst=1,$
                                           xtitle='Pixels',ytitle='Intensity units',position=[.1,.15,.8,.9],/norm,psym=1
                  oploterror,x,toplo,erplo,thick=3,psym=1
                  oplot,x,mygauss_plus_linear(in2,coeff),linestyle=2
                  oplot,x,coeff[3]+coeff[4]*in2
                  xyouts,[0.82,0.82],[0.85,0.85],'GAUSSIAN FIT PARAM',/norm,charsize=1.2
                  xyouts,[0.82,0.82],[0.8,0.8],'Center - '+strtrim(string(i+coeff[1],format='(3f0.3)'),2),/norm,charsize=1.2
                  xyouts,[0.82,0.82],[0.76,0.76],'Half width - '+strtrim(string(coeff[2],format='(3f0.3)'),2),/norm,charsize=1.2
            		  xyouts,[0.82,0.82],[0.72,0.72],'Peak - '+strtrim(string(coeff[0],format='(3f0.3)'),2),/norm,charsize=1.2
            		  xyouts,[0.82,0.82],[0.68,0.68],'Chi!e2!n!dv!n - '+strtrim(string(chisq_red,format='(3f0.3)'),2),/norm,charsize=1.2
            		  xyouts,[0.82,0.82],[0.64,0.64],'Center error - '+string(sigma[1],format='(3f0.2)'),/norm,charsize=1.2

            		  xyouts,[0.82,0.82],[0.55,0.55],'ACCEPTABLE LIMITS',/norm,charsize=1.2
                  xyouts,[0.82,0.82],[0.5,0.5],'Center- '+string(i,format='(3f0.1)')+'pm'+string(shift_cut,format='(3f0.1)'),/norm,charsize=1.2
                  xyouts,[0.82,0.82],[0.46,0.46],'Half width- '+'<'+string(meas_size,format='(3f0.1)'),/norm,charsize=1.2
            		  xyouts,[0.82,0.82],[0.42,0.42],'Peak- '+'>0',/norm,charsize=1.2
            		  xyouts,[0.82,0.82],[0.38,0.38],'Chi!e2!n!dv!n- '+'<'+string(cut_chisq/dof,format='(3f0.1)'),/norm,charsize=1.2
            		  xyouts,[0.82,0.82],[0.34,0.34],'Center error - '+string(1.5,format='(3f0.1)'),/norm,charsize=1.2

                  IF (abs(coeff[1]) LT shift_cut) AND (sigma[1] LT 1.5) AND (coeff[2] lt meas_size) $
                  AND (chisq lt cut_chisq) AND (sigma[1] GT 0.) AND (coeff[0] gt mini) THEN $
                  xyouts,[0.82,0.82],[0.20,0.20],'Meets criteria',/norm,charsize=1.2  
		
                  clearline=fifteenb()
                  form="($,'pause',a,a)"
                  print, form=form, '         ',clearline
	                pause,/quiet
               ENDIF                  
               
                ;For Gaussian fit results to be used the coefficients have to
                ;be less than one pixel from maximum and with an error less than
                ;one pixel. Otherwise position of maximum is used with 0.5 pixel error.
                
                IF (abs(coeff[1]) LT shift_cut) AND (sigma[1] LT 1.5) AND (coeff[2] lt meas_size) $
                AND (chisq lt cut_chisq) AND (sigma[1] GT 0.) AND (coeff[0] gt mini) THEN BEGIN
                    ;Calculate value at peak
                    peak_val=mygauss_plus_linear(coeff[1],coeff)
                    erpeak=mygauss_plus_linear(coeff[1]+sigma[1],coeff)
                    sig_peak=peak_val-erpeak ;estimate of error

                    located.peaks(round(i+coeff[1]),j,0:5)=[i+coeff[1],peak_val,abs(coeff[2]),coeff[3],$
                                                            coeff[4],i]
                    located.errs(round(i+coeff[1]),j,0:5)=[sigma[1],sig_peak,sigma[2],sigma[3],$
                                                            sigma[4],1.]

                   ; maxi.vals(h,0)=i+coeff[1]    & maxi.errs(h,0)=sigma[1] 
                   ; maxi.vals(h,1)=coeff[0]+coeff[3]      & maxi.errs(h,1)=sqrt(sigma[0]^2+sigma[3]^2)
                   ; maxi.vals(h,2)=abs(coeff[2]) & maxi.errs(h,2)=sigma[2]
                   ; maxi.vals(h,3)=coeff[3]      & maxi.errs(h,3)=sigma[3]
                   ; maxi.vals(h,4)=coeff[4]      & maxi.errs(h,4)=sigma[4]
                   ; maxi.vals(h,5)=i             & maxi.errs(h,5)=1. ;????? Check Error
               
                ENDIF ELSE BEGIN
                    located.peaks(i,j,0:2)=[i,mx,0.]
                    located.errs(i,j,0:2)=[0.5,err_dat[i],0.]
                   ; maxi.vals(h,0)=i             & maxi.errs(h,0)=0.5
                   ; maxi.vals(h,1)=mx            & maxi.errs(h,1)=errors[i,j,0]
                   ; maxi.vals(h,2)=0.            & maxi.errs(h,2)=0.
                                    
                ENDELSE
                
                ;h=h+1

             ENDIF

          ENDIF

       ENDFOR


       ;Saves found values
      ; FOR k=0,h-1 DO BEGIN
      ;     located.peaks(maxi.vals(k,0),j,0:5)=maxi.vals(k,0:5)
      ;     located.errs(maxi.vals(k,0),j,0:5)=maxi.errs(k,0:5)
      ; ENDFOR
       
ENDFOR


END
