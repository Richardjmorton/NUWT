;+
; NAME:
;          COUNTER
; PURPOSE:
;          Print a progress status to the screen on a single
;          line. This is WAY cooler than it sounds.
;
; CALLING SEQUENCE:
;          COUNTER, NUMBER, OUTOF  [,INFOSTRING, /PERCENT, WAIT_TIME=variable]
;
; INPUTS:
;          NUMBER:  The current number. Usually the loop index variable
;          OUTOF:   The total number of iterations. The last loop index
;
; OPTIONAL INPUTS:
;
;          INFOSTRING: A string telling the user what is being
;                      counted e.g. 'Flat '
;
; KEYWORD PARAMETERS:
;         
;          PERCENT: Set to output update in percent completed 
;                   percent = rount(number/outof) * 100
;
;          TIMELEFT:  Set to append estimated time remaining.
;          STARTTIME= Used in conjunction w/ /TIMELEFT. Named variable 
;                     that stores the start time of the loop, used for 
;                     calculation of time remaining
;
;          WAIT_TIME:  Used for test and demo purposes only. See
;                      example below.
;
; OUTPUTS:
;          Status is printed to the screen and updated on a single line.
;
; SIDE EFFECTS:
;         This program takes much longer than a simple 
;         print statement. So use COUNTER judiciously. 
;         If your loop consists of only a couple 
;	  of relatively quick commands, updating the 
;	  status with this program could take up a 
;	  significant portion of the loop time!

; PROCEDURE:
;          Put counter statement inside your loop, preferably at the end.
;
; PROCEDURES CALLED: fifteenb.pro
;            
; EXAMPLE:
;          Try this to see how it works:
;
;          IDL> for i = 0,4 do counter,i,4,'test ',wait=.5
;
;
; MODIFICATION HISTORY:
;      Written by JohnJohn, Berkeley 06 January 2003
;  07-Apr-2008 JohnJohn: Finally fixed /TIMELEFT and STARTTIME keywords
;-

pro counter,num,outof,infostring $
            ,wait_time = waittime $
            ,percent=percent      $
            ,clear=clear $
            ,timeleft=timeleft $
            ,starttime=starttime
on_error,2
clearline = fifteenb()          ;get "15b character to create a fresh line
if n_elements(infostring) eq 0 then infostring = 'Number ' 
if keyword_set(clear) then begin
    if keyword_set(timeleft) then begin
        timeinit = {t0: systime(/sec), tot: 0.}
        defsysv,'!time',timeinit
        return
    endif else begin
        len = strlen(infostring)
        print,clearline,format='('+strtrim(len,2)+'x, a)'
        return
    endelse
endif 

case 1 of 
    keyword_set(timeleft): begin
        if n_elements(starttime) eq 0 then begin
            starttime = systime(/sec) 
        endif else begin
            tottime = (systime(/sec) - starttime)
            tave = tottime / float(num)
            tleft = sixty((outof-num) * tave/3600.)
            tleft = strjoin(str(fix(tleft),len=2), ':')
            len = strtrim(strlen(strtrim(tleft,2)),2)
            lenst = strtrim(strlen(infostring),2)
            leni = strtrim(strlen(strtrim(num,2)),2)
            leno = strtrim(strlen(strtrim(outof,2)),2)
            form = "($,a"+lenst+",i"+leni+",' of ',i"
            form += leno+",' Estimated time remaining: ',a"+len+",a,a)"
            print, form=form, infostring, num, outof $
                   , tleft, '         ', clearline
        endelse
    end
    keyword_set(percent) : begin
        per = strtrim(round(float(num)*100./outof),2)
        lenp = strtrim(strlen(strtrim(per,2)),2)
        form="($,a"+lenp+",' % Completed',a,a)"
        print, form=form, per, '         ', clearline
    end
    else : begin
        lenst = strtrim(strlen(infostring),2)
        leni = strtrim(strlen(strtrim(num,2)),2)
        leno = strtrim(strlen(strtrim(outof,2)),2)
        form="($,a"+lenst+",i"+leni+",' of ',i"+leno+",a,a)"
        print, form=form, infostring, num, outof, '         ',clearline
    end
endcase
if n_elements(waittime) gt 0 then wait,waittime
end