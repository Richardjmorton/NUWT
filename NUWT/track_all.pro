; PURPOSE: Locates thread in x-t diagrams, then tries to fit a Gaussian function to cross-sectional flux profile to find
;        signatures of oscillatory behaviour.
;
;
; 
; PROCEDURE: Locates intensity maxima then fits Gaussian to surrounding data points (see locate_things_fg.pro for further 
; details). The individual fits are then joined into threads based on the central position from the Gaussian fit (see 
; follow_thread_fg.pro for further details). Any gaps (for max allowable gap see follow_thread_fg.pro) in threads central 
; positions are joined via linear interpolation and given large uncertainties (see patch_up.pro).
; 
; The central positions are then used to extract the value of Doppler shift at the centre of the thread from the Doppler 
; x-t map.
; 
; The gaps in the width measurements are then filled in with the average value of the measured width for the structure.
; This is temporary and may change. 
;
; The central positions and widths are then used to determine a boundary for the structure. The average Doppler velocity 
; over the cross-section is then determined. Also the Doppler velocity on the left and right hand side of the central peak.
; These will be used in some combination to explore oscillatory features in Doppler.
;  
; The obtained variables (central positions, width & intensity) are then fitted with sinusoidal functions. Designed for 
; measuring physical displacements of structure (kink) and width/intensity variations (sausage).
;
; If doppler x-t diagram included then will also provide average Doppler variations
; across the cross-section and also the Doppler variations on each side of structure
;
;
;
; INPUTS:  data - xt array
; OUTPUTS: results - array of results
;          data_points - array of data points and error bars
;
; OPTIONAL INPUTS: dop_data - xt doppler array
;                  algn_err - provide an estimate for error due to alignment, if given it 
;                            is added to uncertainty in Gaussian fit 
;                  meas_size - sets number of pixels used in Gaussian fit, default is 9, 
;                             odd numbers only                             
;                  errors - array of errors on intensity, i.e. photon noise, 
;                          has to be same size as data array
;                  damped - set if fit of damped sine function is required
;                  area  -  changes size of search box used in follow_thread, default is 5
;                          highly recommended to leave unchanged
;                 
;
;
;CALLS: locate_things_FG.pro, follow_thread_FG.pro, patch_up.pro moscill.pro, dop_central_interp.pro
;
;HISTORY - created R Morton APR 2014
;         R Morton NOV 2014 - super update! Added structure format to remove all the arrays.
;                            Also added COMMON variables so re-used values/structures are passed automatically
;         
;
;
;TO DO: * calculate SD for ave_dop, add calculation of Dop_left & Dop_right + SD
;       * Make decisions about linearly interpolating missing width data points
;       * Decisions needed about averaging for doppler data
;       * How best to use Doppler data?
;       * Doppler data is not yet saved
;       * Plot section needs updating to use structures


PRO display_orig_dat,sl_dt=sl_dt,fnd_dat=fnd_dat,dopps=dopps,xs=xs

window,0
    !p.multi=[0,3,1]
    tvim,sl_dt,title='Slit data',pcharsize=2.0
   ; xyouts,-0.2*xs(1),xs(2)/2,'Slit Data',orientation=90,charsize=2.0

    tvim,fnd_dat,title='Central positions of fitted Gaussians',pcharsize=2.0
  ; xyouts,-0.2*xs(1),xs(2)/4,'Central positions of fitted Gaussians',orientation=90,charsize=1.2

   loadct,70,/silent
   IF n_elements(dopps) GT 0 THEN tvim,dopps,title='Doppler info at central positions',pcharsize=2.0 $
   ELSE plot,sl_dt[*,0],/nodata,xrange=[0,xs(1)],yrange=[0,xs(2)],xst=1,yst=1
   ;xyouts,-0.2*xs(1),xs(2)/4,'Doppler info at central positions',orientation=90,charsize=1.3
   loadct,0,/silent
END

FUNCTION sgn_of,in
   return,fix(in gt 0.)-fix(in lt 0.)
END



pro track_all, data,dop_data=dop_data,algn_err=algn_err,meas_size=meas_size,errors=errors,$
                damped=damped,plot=plot,area=area,tdp=tdp,check=check,fname=fname,kl=kl,dll=dll,dlc=dlc,dlr=dlr

COMMON located_dat,located, nx, nt
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
READ,numt,PROMPT='Enter thread length: '



;#############################################
;Begin search and collect routine
;#############################################
FOR i=0, nslit-1 DO BEGIN
    

    print,'#####################' & print,'Slit number:', i+1,'/'+strtrim(nslit,2) & print,'#####################'
    
   
   IF n_elements(errors) gt 0. THEN errorsi=errors[*,*,i] 
   
   ;++++++++++++++++++++++++++++++++++++++++++++
   ;Locate peaks and fit Gaussian
   ;++++++++++++++++++++++++++++++++++++++++++++
   locate_things_FG,data=data[*,*,i],grad=numg,meas_size=meas_size, $
                 errors=errorsi,check=check

     
   ;++++++++++++++++++++++++++++++++++++++++++++
   ;Add alignment errors to error from Gaussian fits
   ;++++++++++++++++++++++++++++++++++++++++++++
   IF n_elements(algn_err) GT 0 THEN located.errs=located.errs+algn_err

  
   ;++++++++++++++++++++++++++++++++++++++++++++
   ;Connect fitted peaks
   ;++++++++++++++++++++++++++++++++++++++++++++
   follow_thread_fg,min_tlen=numt,area=area
  
   ;++++++++++++++++++++++++++++++++++++++++++++
   ;Fills in 0 elements in the central positions
   ;++++++++++++++++++++++++++++++++++++++++++++ 
   patch_up,fit_flag=0
   
    
  
   ;NEW BIT FOR KINK DOPPLER & ALFVEN WAVE TRACKING!!
   ;++++++++++++++++++++++++++++++++++++++++++++ 
   ;Get values from the Doppler data
   ;++++++++++++++++++++++++++++++++++++++++++++ 
   IF n_elements(dop_data) GT 0 THEN BEGIN
   
     nthreads=n_elements(threads)
     doppler={ave_dop:fltarr(nt), ave_dop_err:fltarr(nt), dop_cen:fltarr(nt), dop_left:fltarr(nt), dop_right:fltarr(nt)}
     phold={ave_dop:fltarr(nt), ave_dop_err:fltarr(nt), dop_cen:fltarr(nt), dop_left:fltarr(nt), dop_right:fltarr(nt)}
     
     FOR ia=1,nthreads-1 DO BEGIN
         
         pos=threads[ia].pos 
         wid=threads[ia].wid 
         
         ;Find time coordinates where there are gaussian fit results
         points=where(pos GT 0.)
         
         ;++++++++++++++++++++++++++++++++++++++++++++ 
         ;Picks out Doppler data at central position - uses interpolation 
         ;to get velocity values at 0.1 pixel accuracy 
         ;++++++++++++++++++++++++++++++++++++++++++++ 
         FOR ij=0,n_elements(points)-1 DO BEGIN

         ;!!!!!!!!!!!###################!!!!!!!!!!!!!!!!!!!!!!!!
         ;STILL NEEDS TESTING HERE - seems okay though
         ;!!!!!!!!!!!###################!!!!!!!!!!!!!!!!!!!!!!!!
         	dop_central_interp,dop_data[*,points[ij]],pos[points[ij]],out_val
         	phold.dop_cen[ij]=out_val
       
         ENDFOR
         ;++++++++++++++++++++++++++++++++++++++++++++ 
         ;Sets 0 width data points to mean value of width
         ;Should this be linear interpolated as with centres?
         ;DOES THIS WORK PROPERLY?
         ;++++++++++++++++++++++++++++++++++++++++++++ 
         
         wid[where(wid EQ 0)]=mean(wid[points(where(wid[points] NE 0))])
         
        
         ;++++++++++++++++++++++++++++++++++++++++++++ 
         ;Define the number of points to the right and left of the
         ;central position to average over - at present 1 sigma
         ;++++++++++++++++++++++++++++++++++++++++++++  
         csecr=round(pos[points]+wid[points])
         csecl=round(pos[points]-wid[points])
         
         
         ;++++++++++++++++++++++++++++++++++++++++++++ 
         ;Calculate average doppler shift over cross-section
         ;Is it best to average in each half or read off a particular value?
         ;++++++++++++++++++++++++++++++++++++++++++++ 
         FOR ip=0,n_elements(points)-1 DO $ 
         phold.ave_dop[ip]=total(dop_data[csecl[ip]:csecr[ip],points[ip]],1)/(csecr[ip]-csecl[ip]+1)

       
         ;++++++++++++++++++++++++++++++++++++++++++++ 
         ;Calculate Maximum doppler shift in section of loop cross-section
         ;Subtract average Doppler shift? 
         ;++++++++++++++++++++++++++++++++++++++++++++ 
         
         FOR ip=0,n_elements(points)-1 DO BEGIN 
 
          ;Calculates the coordinates of normalised cross-section
          ;using PDF for normal distribution.
          f=0.4
          xf=wid[points[ip]]*sqrt(-2.*alog(f));+pos[points[ip]]
          xf2=wid[points[ip]]*sqrt(-2.*alog(f-0.3));+pos[points[ip]]
          mu=(pos[points[ip]])

          mx=max(abs(dop_data[round(mu+floor(xf)):round(mu+ceil(xf2)),points[ip]]),loc)
          phold.dop_right[ip]=sgn_of(dop_data[round(mu+floor(xf))+loc,points[ip]])*mx-dop_cen[ip,ia]
          
          
          indp1=round(mu-ceil(xf2)) & indp2=round(mu-floor(xf))
          IF indp1 LT 0 THEN indp1=0
            
          mx=max(abs(dop_data[indp1:indp2,points[ip]]),loc)
          phold.dop_left[ip]=sgn_of(dop_data[round(mu-ceil(xf2))+loc,points[ip]])*mx-dop_cen[ip,ia]
          
          ;plot,dop_data[*,points[ip]]
          ;oplot,data[*,points[ip]]
          ;plots,[mu+floor(xf),mu+floor(xf)],[0,1000]
          ;plots,[mu+ceil(xf2),mu+ceil(xf2)],[0,1000],linestyle=2
          ;plots,[round(mu-ceil(xf2)),round(mu-ceil(xf2))],[0,1000]
          ;plots,[round(mu-floor(xf)),round(mu-floor(xf))],[0,1000],linestyle=2

          
         ;pause
         ENDFOR    
         doppler=[temporary(doppler),phold]

        plot,doppler[ia].dop_cen,yrange=[-20,20]
        oplot,doppler[ia].ave_dop,psym=1
        oplot,doppler[ia].dop_left,linestyle=1
        oplot,doppler[ia].dop_right,linestyle=3
        pause

         
     ENDFOR
  
       
   
  
   dop_check=dop_data
   dop_check[*]=0.
   her=where(located.pos gt 0.)
   dop_check[her]=dop_data[her]

  ENDIF

   ;kl=pout_ed
   ;dll=dop_left
   ;dlr=dop_right
   ;dlc=dop_cen   
   
    ;display_orig_dat,sl_dt=data[*,*,i],fnd_dat=out[*,*,i],dopps=dop_check,xs=sz

    window,1
    !p.multi=[0,1,1]



   dsk='y'
   READ,dsk,PROMPT='Do you want to fit the transverse motions?: (y/n) '
   IF dsk EQ 'y' THEN back=fitting_them(1,damped=damped,tdp=tdp)
   
   dsk='y'
   READ,dsk,PROMPT='Do you want to fit the intensity?: (y/n) '
   IF dsk EQ 'y' THEN back=fitting_them(2,damped=damped,tdp=tdp)
    
   dsk='y'
   READ,dsk,PROMPT='Do you want to fit the width?: (y/n) '
   IF dsk EQ 'y' THEN back=fitting_them(3,damped=damped,tdp=tdp)
   
   


    ;PLOTTING ROUTINE
    IF keyword_set(plot) THEN BEGIN

         px=1d
         READ,px,PROMPT='Enter pixel size in km: '
         tx=1d
         READ,tx,PROMPT='Enter time in (s): '
         file='a'
         READ,file,PROMPT='Provide save directory + save filename for plotting, i.e., directory/filename: '  
    
     FOR ii=1,nfits-1 DO BEGIN
         sz=size(data)   
         
         t=tx*findgen(sz(2)-1)
         t1=res[11,ii]
         t2=res[12,ii]
         time_len=t2-t1
         px2=px/1000.
         outno=count[ii]
         wave=px2*(res[0,ii]+res[4,ii]*findgen(30))+ $
              px2*res[1,ii]*sin(2.*!pi*findgen(time_len)/res[2,ii]-res[3,ii])
         p2p=max(wave)-min(wave)       
         IF P2P GT 2.*px2 THEN mx=max(wave)+P2P ELSE mx=max(wave)+2.*px2
         IF P2P GT 2.*px2 THEN mn=min(wave)-P2P ELSE mn=min(wave)-2.*px2

         set_plot,'ps'
         !p.font=0
         device,/encapsul,filename=file+'slit_'+strtrim(i,2)+'_fit_'+strtrim(ii,2)+'.eps',helvetica=1
         
         plot,t(t1:t2),wave,xst=1,yst=1,yrange=[mn,mx],xtitle='Time (s)',ytitle='Distance (Mm)'
         oplot,t(t1:t2),px2*out3[t1:t2,outno,0],psym=1 
         oploterr,t(t1:t2),px2*out3[t1:t2,outno,0],px2*out3[t1:t2,outno,1],psym=1
            
         device,/close
         set_plot,'x'

      ENDFOR 
    ENDIF

    yn='a'
    READ,yn,PROMPT='Save results? If the results are not saved now they will be overwritten! y/n '
    IF yn EQ 'y' THEN BEGIN
      IF n_elements(threads_fit_fg) EQ 0 THEN BEGIN
         IF n_elements(fname) EQ 0 THEN BEGIN
		  sv_dir='a'
       		   READ,sv_dir,PROMPT='Give directory + filename for saving data points and fits. Saved as directory/filename_N.idl: '
    	 ENDIF ELSE sv_dir=fname
         save,threads,filename=sv_dir+'_'+strtrim(i,2)+'.idl'

      ENDIF ELSE BEGIN
      	 IF n_elements(fname) EQ 0 THEN BEGIN
		  sv_dir='a'
       		   READ,sv_dir,PROMPT='Give directory + filename for saving data points and fits. Saved as directory/filename_N.idl: '
    	 ENDIF ELSE sv_dir=fname
         save,threads,threads_fit_fg,filename=sv_dir+'_'+strtrim(i,2)+'.idl'
      ENDELSE

    ENDIF

ENDFOR

END
