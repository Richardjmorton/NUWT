;FUNCTION: Locates peaks in time-distance plot with pixel accuracy -
;          reduced version of 'locate_min.pro'

;PROCEDURE OUTLINE: Takes an unsharp masked (and high-pass filtered)
;time-distance diagram and locates the peaks using a find_max crawling
;routine.The maximum then has to have a specific gradient (default
;>0.5) to be classed as a peak. 

;INPUTS: data - time-distance diagram
;        grad - limits on gradient (default = 0.5)
;
;OPTIONAL INPUTS: 
;                 simp_grad - Use of an analytic least-sqaures estimate for gradient - may be faster!   
;
;OUTPUTS: out - array showing nearest whole pixel positions of located
;               maximum
;               values (out[*,*,0]) and value of maximum (out[*,*,1]).
;         errs - array showing errors on fits to position
;
;
;
;Author: R. Morton Oct - 2012
;        R Morton NOV 2014 - super update! Added structure format to remove all the arrays.
;                            Also added COMMON variables so re-used values/structures are passed automatically
;        R Morton MAR 2015 - fixed bug with definition of initial minimum value to set all located.peaks values to
;        14 MAR 2016 - Release as version 1.0 of NUWT 
;        R Morton MAY 2016 - Added option for analytic LS calc for graident in hope of speed increase
;                            Also removed needless for loop at bottom!!

pro locate_things_min,data=data,grad=grad,simp_grad=simp_grad,errors=errors

COMMON located_dat,located, nx, nt

sz=size(data)
nx=sz(1) & nt=sz(2)
located={peaks:fltarr(nx,nt,2), errs:fltarr(nx,nt)}
located.peaks(*,*,*)=(min(data)-10.)<(-10) ;Should ensure no confusion
;maxi={vals:fltarr(100,2), errs:fltarr(100,2)}


;Sets default gradient if non specified
IF n_elements(grad) EQ 0 THEN grad=0.5

;IF NO ERRORS GIVEN USE WEIGHTS EQAUL TO 1
IF n_elements(errors) LE 1 THEN BEGIN
   errors=fltarr(nx,nt)
   errors[0:nx-1,0:nt-1]=1.
ENDIF
;Begin search in each time-slice
in=[-2,-1,0,1,2]

FOR j=0, (nt-1) DO BEGIN
    ;h=0
    image=reform(data[0:nx-1,j,0])
    err_dat=reform(errors[0:nx-1,j,0])

    ;Search each along x-direction, excludes edges
    FOR i=6, (nx-6) DO BEGIN

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
	   IF (m1 GT grad) AND (m2 LT (-1.)*grad) THEN BEGIN
	      ;maxi.vals(h,0:1)=[i,mx]
	      ;maxi.errs(h,0)=0.5
	      ;h=h+1
          located.peaks(i,j,0:1)=[i,mx]
	      located.errs(i,j)=0.5

	   ENDIF

	ENDIF

     ENDFOR

     ;Saves found values
    ; FOR k=0,h-1 DO BEGIN
	;located.peaks(maxi.vals(k,0),j,0)=maxi.vals(k,0)
	;located.peaks(maxi.vals(k,0),j,1)=maxi.vals(k,1)
	;located.errs(maxi.vals(k,0),j)=maxi.errs(k,0) 
     ;ENDFOR
ENDFOR



END
