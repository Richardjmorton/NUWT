;FUNCTION: Locates peaks in time-distance plot with pixel accuracy -
;          reduced version of 'locate_min.pro'

;PROCEDURE OUTLINE: Takes an unsharp masked (and high-pass filtered)
;time-distance diagram and locates the peaks using a find_max crawling
;routine.The maximum then has to have a specific gradient (default
;>0.5) to be classed as a peak. 

;INPUTS: data - time-distance diagram
;        grad - limits on gradient (default = 0.5)
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


pro locate_things_min,data=data,grad=grad

COMMON located_dat,located, nx, nt

sz=size(data)
nx=sz(1) & nt=sz(2)
located={peaks:fltarr(nx,nt,2), errs:fltarr(nx,nt)}
located.peaks(*,*,*)=(min(data)-10.)<(-10) ;Should ensure no confusion
maxi={vals:fltarr(100,2), errs:fltarr(100,2)}


;Sets default gradient if non specified
IF n_elements(grad) EQ 0 THEN grad=0.5

;Begin search in each time-slice
FOR j=0, (nt-1) DO BEGIN
    h=0
    image=data[0:nx-1,j,0]

    ;Search each along x-direction, excludes edges
    FOR i=6, (nx-6) DO BEGIN

        ;Finds maximum
        mx=max(image[i-5:i+5,0,0],loc)

        ;Wait till maximum is at centre of search bar
	IF loc EQ 5 THEN BEGIN
	   
  	   ;Finds gradients either side of the maximum
	   in=[-2,-1,0,1,2]
	   res=poly_fit(in,image[i-4:i,0,0],2,yfit=fit)
	   res2=poly_fit(in,image[i:i+4,0,0],2,yfit=fit2)

	   ;If gradients greater than a certain value begin
	   IF (res[1] GT grad) AND (res2[1] LT (-1.)*grad) THEN BEGIN

	      maxi.vals(h,0)=i
	      maxi.vals(h,1)=mx
	      maxi.errs(h,0)=1.
	      h=h+1

	   ENDIF

	ENDIF

     ENDFOR

     ;Saves found values
     FOR k=0,h-1 DO BEGIN
	located.peaks(maxi.vals(k,0),j,0)=maxi.vals(k,0)
	located.peaks(maxi.vals(k,0),j,1)=maxi.vals(k,1)
	located.errs(maxi.vals(k,0),j)=maxi.errs(k,0) 
     ENDFOR
ENDFOR



END
