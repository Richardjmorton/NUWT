
;PURPOSE - Create a sample time-distance diagram of a coronal loop that looks like data from imagers. 
;          No physics is included. Structures recreated using measured properties from observed coronal loops. 
;          Intended for testing accuracy of Gaussian fitting technique.
;
;INPUTS - maxint - maximum intensity for in DN for loop  
;         width - standard deviation from Gaussian fit to loop cross-section 
;         amp - amplitude of wave in pixels
;         period - period of wave in seconds
;         rsltn - resolution in arcseconds
;         cad - cadence in seconds
;         stand_dev - sigma value for noise (noise is randomly picked from normally distributed values)
;         mea - mean value of noise
;         seed - sets a seed value normally distributed random numbers 
;         aligner - set sigma of alignment error (alignment error is randomly picked from normally distributed values)         
;
;OUTPUTS - karr - provides 2d array of kink wave
;          middle - returns central position of gaussian (20+middle)
;
;OPTIONAL INPUTS - clean - set keyword to set noise to a low level
;                  plot - set keyword to plot the generated array to screen
;                  background - set keyword to allow prompt and setting of background function
;
;TO DO - add some way to change initial array and output array that depend on resolution and cadence
;      - add optional extra threads, i.e., multi-thread
;
; Created R Morton - 2013



;###############################################
;Create large array to create a well resolved Gaussian profile whose central position varies with time
;Can add additional 'sausage' mode which introduces intensity and pi out of phase width variation
;
;############################################### 
FUNCTION do_the_kink,maxint,xk,width,y,rndarr,saus_amp,wid_amp,period_saus,saus
  
  kink=fltarr(2000,600)
  
  IF n_elements(saus) LT 1  THEN BEGIN
         FOR i=0, 599 do kink[*,i]=maxint*exp(-((findgen(2000)/50.-(xk[i]+20.))/width)^2/2.)
  ENDIF ELSE BEGIN

         SA=maxint+saus_amp*sin(2.*!Pi*(y)/period_saus)
         WD=width+wid_amp*sin(2.*!Pi*(y)/period_saus+!pi)

         FOR i=0, 599 do kink[*,i]=SA[i]*exp(-((findgen(2000)/50.-(xk[i]+20.))/WD[i])^2/2.)
  ENDELSE 

  RETURN,kink
END

;###############################################
;Create large array to create a well resolved uniform Doppler shifts whose central position varies with time
;############################################### 

FUNCTION do_the_dop,xk,width,y,rndarr,xk2,pc_area
  
 
  main=fltarr(2000,600)
  FOR i=0,599 DO BEGIN
      main[*,i]=exp(-((findgen(2000)/50.-(xk[i]+20.))/width)^2/2.)
      ind=where(main[*,i] gt pc_area)
      main[*,i]=0
      main[ind,i]=1.
     
      main[ind,i]=xk2[i]*main[ind,i]     
      ;exp(-((findgen(2000)/50.-(xk[i]+20.))/width)^2/2.)
  ENDFOR 

  dop_vel=main
 
  RETURN,dop_vel
END

FUNCTION do_the_alf,xk,width,y,rndarr,xTA,pc_area
  
  main=fltarr(2000,600)
    
  FOR i=0,599 DO BEGIN
      main[*,i]=exp(-((findgen(2000)/50.-(xk[i]+20.))/width)^2/2.)
      ind=where(main[*,i] gt pc_area)
      main[*,i]=0
      main[ind,i]=1.
      D_across=cos(!pi*findgen(n_elements(ind))/(n_elements(ind)-1))
      main[ind,i]=d_across*xTA[i]      
   
 ENDFOR  
  
  dopta=main    
      

  RETURN,dopTA
END

pro create_kink, maxint=maxint,width=width,amp=amp,period=period,rsltn=rsltn,cad=cad,stand_dev=stand_dev,$
                 mea=mea,seed=seed,algnerr=algnerr,background=background,kink=kink,$  
                 middle=middle,plot=plot,clean=clean,dop_vel=dop_vel,dopTA=dopTA,do_saus=do_saus,$
                 pc_area=pc_area,comb_dop=comb_dop



;Properties of data/structure
IF n_elements(maxint) EQ 0 THEN maxint=700.0 ;Value for a loop in AIA 171A (800+300 of bck)
IF n_elements(width) EQ 0 THEN width=1.2    ; Value for a loop in AIA 171A

IF n_elements(rsltn) EQ 0 THEN rsltn=0.6 ; AIA 171 plate scale
IF n_elements(cad) EQ 0 THEN cad=12. ; Typical AIA 171 cadence
IF n_elements(stand_dev) EQ 0 THEN stand_dev=18.37  ; Results taken from fits to noise 
IF n_elements(mea) EQ 0 THEN mea=-1.59              ; around a loop - see details given in 'data/kink_sim/'
IF n_elements(seed) EQ 0 THEN seed=1L
IF n_elements(algnerr) EQ 0 THEN algnerr=0.05                            
IF keyword_set(clean) THEN stand_dev=1.  ; Sets a 'clean' value, i.e., very low noise
IF KEYWORD_SET(do_saus) THEN saus=1.   



;print,'Amplitude - '+strtrim(amp,2),' Period - '+strtrim(period,2),' Resolution-',rsltn,' Cadence-',cad  







y=findgen(600)
rndarr=0;algnerr*randomn(2.*seed,600)


;###############################################
;Properties of Kink mode in plane

if not keyword_set(amp) then amp=1.0
if not keyword_set(period) then period=300.

saus_amp=300.
wid_amp=0.5
period_saus=150.

!p.multi=[0,2,1]

xk=amp*sin(2.*!Pi*(y+rndarr)/period)
kink=do_the_kink(maxint,xk,width,y,rndarr,saus_amp,wid_amp,period_saus,saus)


;Middle position of Gaussian
middle=fltarr(50)
for i=0,49 do middle[i]=xk(12.*i)

;###############################################
;Properties of Kink mode in LOS
ampk2=4.
periodk2=100.
xk2=ampk2*sin(2.*!Pi*(y+rndarr)/periodk2)

IF n_elements(pc_area) EQ 0 THEN pc_area=0.2 ;Related to how much of the Gaussian shape to use

dop_vel=do_the_dop(xk,width,y,rndarr,xk2,pc_area)


!p.multi=0

;###############################################
;Properties of Torsional Alfven mode
dopTA=fltarr(2000,600)
ampk2=20.
periodTA=200.
xTA=ampk2*sin(2.*!Pi*(y+rndarr)/periodTA)

dopTA=do_the_alf(xk,width,y,rndarr,xTA,pc_area)



;size for new array
sz=intarr(2)
sz(0)=40
sz(1)=50

rndarr2=stand_dev*randomn(seed,sz(0),sz(1))+mea

comb_dop=dop_vel+dopTA

kink=congrid(temporary(kink),sz(0),sz(1))
dop_vel=congrid(temporary(dop_vel),sz(0),sz(1))
dopTA=congrid(temporary(dopTA),sz(0),sz(1))
comb_dop=congrid(temporary(comb_dop),sz(0),sz(1))

;POISSON NOISE - calculated in DN. For Poisson distribution mean=variance=F/Gain 
;pnoise=poidev((kink2/gain))  
;error=sqrt(pnoise+addnoise)
kink=kink+smooth(rndarr2,3,/edge_wrap)




if keyword_set(background) then begin
    
    print,'#################################'
    print,'Configure the background details'
    print,'%'
    print,'%'
 
    endlp='y'
    ;while loop used to repeat whole fitting procedure
    while (endlp eq 'y') do begin
      
      
      plod=1d
      READ,plod,PROMPT='Enter polynomial order (1-3): '
     
      bckg=fltarr(plod+1)
     
      for i=0,plod do begin

       num=1d
       print,'x^',i 
       READ,num,PROMPT='Enter coefficient: '
       bckg[i]=num

      endfor   
     
      x=findgen(sz(0))
      if plod eq 1 then begin
         plot,kink[*,20]+bckg[1]*x+bckg[0]
         for i=0,sz(1)-1 do kink[*,i]=kink[*,i]+bckg[1]*x+bckg[0]
      endif     
      if plod eq 2 then begin
         plot,kink[*,20]+bckg[2]*x^2+bckg[1]*x+bckg[0]
         for i=0,sz(1)-1 do kink[*,i]=kink[*,i]+bckg[2]*x^2+bckg[1]*x+bckg[0]
      endif
      if plod eq 3 then begin
         plot,kink[*,20]+bckg[3]*x^3+bckg[2]*x^2+bckg[1]*x+bckg[0]
         for i=0,sz(1)-1 do kink[*,i]=kink[*,i]+bckg[3]*x^3+bckg[2]*x^2+bckg[1]*x+bckg[0]
      endif

      READ,endlp,PROMPT='Change the background? - y/n '
           
     endwhile
    
      
endif 

if keyword_set(plot) then begin
    
    set_plot,'x'
    window,0
    tvim,rotate(kink,3),/noaxis
    loadct,13,/silent
    oplot,19-middle,color=250
    loadct,0,/silent
    axis,0,yaxis=0,yrange=rsltn*0.735*[0,40],yst=1,ytitle='Distance (Mm)'
    axis,0,xaxis=0,xrange=cad*[0,50],xst=1,xtitle='Time (s)'
    
endif

end
