FUNCTION fitting_them,fit_flag,damped=damped,tdp=tdp

    COMMON located_dat, located, nx, nt
    COMMON threads_dat, threads
    COMMON extended_dat, threads_fit_fg

    IF fit_flag EQ 0 THEN BEGIN

    ENDIF ELSE BEGIN

   	 s=''
    	CASE fit_flag OF
    	   1:s='NOW FITTING POSITION VARIATIONS'
   	   2:s='NOW FITTING INTENSITY VARIATIONS'
   	   3:s='NOW FITTING WIDTH VARIATIONS'
   	   4: ;Something to do with doppler data
       	   ELSE:res=0
    	ENDCASE
     	print,'##################################'
     	print,s
     	print,'##################################'
    
        IF NOT keyword_set(damped) THEN BEGIN
       		 moscill,fit_flag=fit_flag
        ENDIF ELSE BEGIN
           IF NOT keyword_set(tdp) THEN BEGIN
 	          moscill,fit_flag=fit_flag,/damped
           ENDIF ELSE BEGIN
                moscill,fit_flag=fit_flag,/damped,/tdp
           ENDELSE
        ENDELSE

     	CASE fit_flag OF
       		1:res=where(threads_fit_fg[*].fit_result_pos[1] ne 0)
       		2:res=where(threads_fit_fg[*].fit_result_inten[1] ne 0)
      	 	3:res=where(threads_fit_fg[*].fit_result_wid[1] ne 0)
      	 	4: ;Something to do with doppler data
      	 	ELSE:res=0
     	ENDCASE
     
     nfits=n_elements(res)
     print,'No. fits saved',nfits-1
    ENDELSE

return,nfits
END