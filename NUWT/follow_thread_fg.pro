;
;PURPOSE: Joins peaks in intensity found with locate_things_fg.pro to create 'threads'
;
;INPUTS: - Called via COMMON block defined in locate_things_fg
;
;OUTPUTS: - threads - structure containing all the information on the found threads
;           threads.pos(*,0) & threads.errs(*,0) - arrays containing central positions and errors
;           threads.pos(*,1) & threads.errs(*,1) - arrays containing central intensities and errors 
;           threads.pos(*,2) & threads.errs(*,2) - arrays containing Gaussian widths and errors
;                Threads beginning and end marked by minimum intensity value of the array minus 10
;                Where routine misses pixel in thread a 0. value is recorded
;
;OPTIONAL INPUTS: - min_tlen - change minimum thread length, default is 30
;	          - area - change size of search area used for locating neighbouring pixels
;                          default is 5 by 5 box
;
;HISTORY: 2012 - Created by RJM
;         2014 APR - updated code so 'out' array is same size as No. features found
;         R Morton NOV 2014 - super update! Added structure format to remove all the arrays.
;                            Also added COMMON variables so re-used values/structures are passed automatically

;
;REMINDER: located.peaks[*,*,0] is central positions, located.peaks[*,*,1] is maximum intensities, located.peaks[*,*,2] is widths
;
;

pro follow_thread_fg,min_tlen=min_tlen,area=area,debug=debug

;Loads common data generated from locate_things_fg
COMMON located_dat, located, nx, nt
COMMON threads_dat, threads

mini=min(located.peaks[0:nx-1,0:nt-1,0])


;Set minimum thread length if nothing specified
IF n_elements(min_tlen) EQ 0 THEN min_tlen=30

;Sets size of search area
IF NOT keyword_set(area) THEN area=5

; set up dummy arrays
image=located.peaks[0:nx-1,0:nt-1,0];Select binary x-t array of central positions for thread following
im_err=located.errs[0:nx-1,0:nt-1,0]

;Format of structures and arrays for saving
threads={pos:fltarr(nt), err_pos:fltarr(nt), inten:fltarr(nt), err_inten:fltarr(nt), wid:fltarr(nt), err_wid:fltarr(nt)}
posi=fltarr(nt,3) & pos_err=fltarr(nt,3)
tnum=0


kk=0

;Search over time
FOR j=0,nt-(area+1) DO BEGIN 
        ; Search over x
	FOR i=area,nx-(area+1) DO BEGIN 

	    tlen=0.
            ; reset temporary loop arrays 
            posi[0:nt-1,0:2]=0. & pos_err[0:nt-1,0:2]=0.
          
            ;Values before thread are set equal to -1
            posi(0:j,0:2)=-1. & pos_err(0:j,0:2)=-1.
            h=j+1

	    IF image[i,j] GT mini THEN BEGIN
                
		val=fix(area)/2 ; old - area-2

		;create box that's (2*val+1) by area-1
		a=image[(i-val):(i+val),h:h+(area-1)]
		k=i

		;set up exit from loop - if no pixel with value greater than minimum
		;                        in box then end
		WHILE max(a) GT mini DO BEGIN

		      ;determines how many points have values greater than minimum
		      b=where(a gt mini)

		      IF b(0) LT 0 THEN BEGIN
			 a=mini ;if no points then set a to minimum value so loop is 
                                ;broken
		      ENDIF ELSE BEGIN

			;find coordinates of first non-minimum point in the box
		      	xm=b(0) mod (val*2+1)
		      	ym=b(0)/(val*2+1)

		      	;saves coordinates then erases them from dummy array 
                      	;so values not used twice
                     
		      	k=k-val+xm

                       posi(h+ym,0)=reform(image[k,h+ym,0]) & pos_err(h+ym,0)=reform(im_err[k,h+ym,0])
                       posi(h+ym,1)=reform(located.peaks[k,h+ym,1]) & pos_err(h+ym,1)=reform(located.errs[k,h+ym,1])
                       posi(h+ym,2)=reform(located.peaks[k,h+ym,2]) & pos_err(h+ym,2)=reform(located.errs[k,h+ym,2])
                      
		      	image[k,h+ym]=mini
		      	im_err[k,h+ym]=mini

			IF ym EQ 0 THEN ym=1
			h=h+ym


			;stops loop from crashing as it reaches end and
			;sets box (a) for next loop
                      
			IF h GT nt-(area-1) THEN a=mini ELSE a=image[(k-val):(k+val),h:h+(area-2)]

			tlen=tlen+1

		     ENDELSE

		ENDWHILE

               ;Any values after the thread are set to -1
		posi(h:nt-1,0:2)=-1 & pos_err(h:nt-1,0:2)=-1
   		

		IF tlen GT min_tlen THEN BEGIN
                 phold={pos:posi(*,0),err_pos:pos_err(*,0), inten:posi(*,1), err_inten:pos_err(*,1), $
                        wid:posi(*,2), err_wid:pos_err(*,2)}
                 threads=[temporary(threads),phold]
		   IF keyword_set(debug) THEN BEGIN
                   kk=kk+1
                   !p.multi=[0,1,3]
                   window,1
                   plot,threads[kk].pos[*],xtitle='Time',ytitle='Central position',charsize=2
                   oploterror,threads[kk].pos[*],threads[kk].err_pos[*],psym=1
                   plot,threads[kk].inten[*],xtitle='Time',ytitle='Maximum Intensity',charsize=2
                   oploterror,threads[kk].inten[*],threads[kk].err_inten[*],psym=1
                   plot,threads[kk].wid[*],xtitle='Time',ytitle='Width',charsize=2
                   oploterror,threads[kk].wid[*],threads[kk].err_wid[*],psym=1
                   print,threads[kk].err_inten[*]
                   pause
                  ENDIF
		ENDIF

	  ENDIF
     ENDFOR
ENDFOR

!p.multi=0

END
