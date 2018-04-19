;+
;NAME: FOLLOW_THREADS
;
;PURPOSE:
;   Find and follow "threads" of local intensity peaks identified by "locate_things.pro"
;
;INPUTS:
;   Main input is the "located_dat" COMMON block with the peak locations
;
;OPTIONAL INPUTS:
;   min_tlen - minimum thread length to be saved in the output structures
;   max_dist_jump - maximum allowable distance (in pixels) between two peaks in
;                   the same thread. Default is 3 pixels
;   max_time_skip - maximum allowable timesteps beweewn two peaks in the same
;                   thread. Default is 4 time steps.
;   scan_dir - sets the direction in which the program scans the search box.
;              Depending on the input array of peaks, this may bias the results.
;              May choose from the following four options:
;                "outward" - [DEFAULT] bias towards smaller displacements
;                "inward" - bias towards larger displacements
;                "right" - bias towards features that slope to the LEFT
;                "left" - bias towards features that slope to the RIGHT
;   /check - if set, will plot each thread after saving it. Should only be used
;            for testing the code and NOT general operations.
;   /debug - if set, will create an array of the same size as td-diagram which assigns
;            each pixel a number corresponding to its thread (including rejected threads).
;            Saves the array in the 'th_debug_dat' COMMON block
;
;OUTPUTS:
;   threads - array of structures containing the threads found and saved to the
;             "threads_dat" COMMON block. Each thread structure (array element)
;             has the following format:
;               .pos - [nt] long array with the thread position at each timestep.
;                      Values of -1 indicate timesteps outside of the thread and
;                      a 0 indicates timesteps skipped due to no nearby peaks
;               .err_pos - [nt] long array with the thread position errors
;               .bin_flags - [nt] array with flags indicating the type of peak
;                            found in each timestep of the thread. Possible values are:
;                               -1 : time-step not part of thread
;                                0 : data gap inside thread
;                                1 : lower quality data found (not currently used)
;                                    (may be used for filling with rejected peaks)
;                                2 : higher quality data found
;               .start_bin - timestep bin with the first thread position
;               .end_bin - timestep bin with the last thread position
;               .length - total length (in timesteps) of the thread
;             If "locate_things.pro" was run with the /FULL_GAUSS option set,
;             the output structures will also include the following tags:
;               .inten - [nt] long array with the peak intensity values
;               .err_inten - [nt] long array with the estimated intensity errors
;               .wid - [nt] long array with the gaussian width found by the
;                      subpixel fitting method of "locate_things.pro"
;               .err_wid - [nt] long array with the gaussian fit errors
;
;HISTORY: Name---------Date---------Description
;         R Morton  OCT, 2012 - Initial coding
;         R Morton  ???, 2014 - updated code so 'out' array is same size as No. features found
;         R Morton  NOV, 2014 - super update! Added structure format to remove all
;                               the arrays. Also added COMMON variables so re-used
;                               values/structures are passed automatically
;         R Morton  APR, 2015 - moved some IF statements from moscill.pro to here
;         R Morton  14 MAR, 2016 - Released as version 1.0 of NUWT
;         M Weberg  JULY, 2016 - Modified the program to avoid saving blank results
;                                as the first "thread"
;         M Weberg  SEPT, 2016 - moderate update:
;                                 - fixed a few small bugs
;                                 - if multiple peaks are found in the same timestep
;                                   of the search box, will now select the peak with
;                                   the smallest shift in position
;                                 - added tracking of start_bin, end_bin, and length
;         M Weberg   JAN, 2017 - added the 'debug' keyword flag
;         M Weberg   DEC, 2017 - added the 'scan_dir' options
;         M Weberg   JAN, 2018 - merged in capabilities of "follow_threads_fg.pro"
;                                Will automatically detect when "locate_things.pro"
;                                was run using the /FULL_GAUSS option and use the
;                                correct output array
;         M Weberg   MAR, 2018 - Improved the behavior of the tracking method when
;                                it nears the sides or end of the td-diagram
;                    27 MAR 2018 - Release as version 2.0 of NUWT
;
;TO DO:
;   - none at the moment
;-

PRO FOLLOW_THREADS, min_tlen=min_tlen, max_dist_jump=max_dist_jump, $
                         max_time_skip=max_time_skip, scan_dir=scan_dir, $
                         check=check, debug=debug

COMPILE_OPT IDL2
;@defaults######################################################################
;SETTING DEFAULT VALUES
;###############################################################################
;Loads common data generated from locate_things
COMMON located_dat, located
COMMON threads_dat, threads

mini = min(located.peaks[*,*,0]) ; min data val possible (normally the 'invalid data' fill value)

;Set minimum thread length to output
IF N_ELEMENTS(min_tlen) EQ 0 THEN min_tlen = 20 ;old default of 30 was too high

;Setting variables which control the search box size
IF NOT KEYWORD_SET(max_dist_jump) THEN max_dist_jump = 3
IF NOT KEYWORD_SET(max_time_skip) THEN max_time_skip = 4
IF max_time_skip LT 1 THEN max_time_skip = 1
IF typename(max_dist_jump) NE 'INT' THEN max_dist_jump = fix(max_dist_jump)
IF typename(max_time_skip) NE 'INT' THEN max_time_skip = fix(max_time_skip)
search_box_width = 2*max_dist_jump + 1
search_box_height = max_time_skip

;Setting default scan direction
IF NOT KEYWORD_SET(scan_dir) THEN scan_dir = 'outward'
scan_dir = STRLOWCASE(scan_dir)
IF TOTAL(STRMATCH(['outward', 'inward', 'right', 'left'], scan_dir)) LE 0 THEN BEGIN
    PRINT, 'WARNING: unknown scan direction! Defaulting to "outward".'
    scan_dir = 'outward'
ENDIF

;Set up dummy arrays
image = located.peaks[*,*,0]
im_err = located.errs[*,*,0]

sz = size(located.peaks)
nx = sz[1]
nt = sz[2]

num_pk_layers = sz[3] ;num vals output by 'locate_things.pro' at each peak
IF num_pk_layers EQ 6 THEN full_gauss = 1 ELSE full_gauss = 0

;@int_output####################################################################
;INITIALIZE THE "OUTPUT" ARRAY OF STRUCTURES AND FILLING WITH DEFAULT VALUES
;###############################################################################
;Initialize the output array of structures (prevents crashes if no threads are found)
IF full_gauss EQ 0 THEN BEGIN
    ;[Default] Only track peak locations
    threads = {pos:fltarr(nt), err_pos:fltarr(nt), $
               bin_flags:intarr(nt), $
               start_bin:-1, end_bin:-1, length:-1}

    ;Temporary structure to hold the results before appending to output array (or not)
    temp_th = {pos:fltarr(nt), err_pos:fltarr(nt), $
               bin_flags:intarr(nt), $
               start_bin:-1, end_bin:-1, length:1}
ENDIF ELSE IF full_gauss EQ 1 THEN BEGIN
    ;For use with the "/full_gauss" mode of "locate_things.pro"
    threads = {pos:fltarr(nt), err_pos:fltarr(nt), $
               inten:fltarr(nt), err_inten:fltarr(nt), $
               wid:fltarr(nt), err_wid:fltarr(nt), $
               bin_flags:intarr(nt), $
               start_bin:-1, end_bin:-1, length:-1}

    ;Temporary structure to hold the results before appending to output array (or not)
    temp_th = {pos:fltarr(nt), err_pos:fltarr(nt), $
               inten:fltarr(nt), err_inten:fltarr(nt), $
               wid:fltarr(nt), err_wid:fltarr(nt), $
               bin_flags:intarr(nt), $
               start_bin:-1, end_bin:-1, length:-1}
ENDIF

;Quick key for values in '.bin_flags':
;   -1 : timestep not part of thread
;    0 : data gap inside thread
;    1 : lower quality data found (currently not used but, in the future,
;        may be used for filling with lower quality or rejected peaks)
;    2 : higher quality data found

IF KEYWORD_SET(debug) THEN BEGIN
    COMMON th_debug_dat, th_debug
    th_debug = {raw_th_num:fltarr(nx, nt)}
ENDIF

;@set_scan_dir##################################################################
;SET THE SEARCH DIRECTION
;###############################################################################
;Reorder the indices for scanning the search box values.
IF scan_dir EQ 'outward' THEN BEGIN
    ;Scan FROM the center TOWARDS the outside edges
    ;Biased towards selecting whichever peak gives the SMALLEST displacement
    reorder_cols = INTARR(search_box_width)
    sub_ind = INDGEN(max_dist_jump)*2 + 1 ;every other index index = 1
    abs_col_shift = INDGEN(max_dist_jump) + 1 ;relative col shift from the center
    reorder_cols[0] = max_dist_jump ;set the middle col (no dist change) to the START
    reorder_cols[sub_ind] = max_dist_jump - abs_col_shift ;cols to the left of center
    reorder_cols[sub_ind+1] = max_dist_jump + abs_col_shift ;cols to the right
ENDIF ELSE IF scan_dir EQ 'inward' THEN BEGIN
    ;Scan FROM the outside edges TOWARDS the center
    ;Biased towards selecting whichever peak gives the LARGEST displacement
    reorder_cols = INTARR(search_box_width)
    sub_ind = INDGEN(max_dist_jump)*2 ;every other index starting at index = 0
    abs_col_shift = INDGEN(max_dist_jump) + 1 ;relative col shift from the center
    reorder_cols[-1] = max_dist_jump ;set the middle col (no dist change) to the END
    reorder_cols[sub_ind] = REVERSE(max_dist_jump - abs_col_shift) ;cols to the left of center
    reorder_cols[sub_ind+1] = REVERSE(max_dist_jump + abs_col_shift) ;cols to the right
ENDIF ELSE IF scan_dir EQ 'right' THEN BEGIN
    ;Scan FROM the left TO the right
    ;Biased towards selecting features that slope TOWARDS the LEFT
    reorder_cols = INDGEN(search_box_width)
ENDIF ELSE IF scan_dir EQ 'left' THEN BEGIN
    ;Scan FROM the right TO the left
    ;Biased towards selecting features that slope TOWARDS the RIGHT
    reorder_cols = REVERSE(INDGEN(search_box_width))
ENDIF

;Create look-up arrays for finding pixel shift values (cleaner and easier to follow)
dist_pxl_shift = INDGEN(search_box_width) - max_dist_jump
dist_pxl_shift = REBIN(dist_pxl_shift, search_box_width, search_box_height)
dist_pxl_shift = dist_pxl_shift[reorder_cols, *]

time_pxl_shift = INDGEN(search_box_height) + 1
time_pxl_shift = REBIN(REFORM(time_pxl_shift, 1, search_box_height), $
                       search_box_width, search_box_height)

;@track#########################################################################
;TRACK THREADS IN TIME
;###############################################################################
last_timestep = max([max_time_skip, min_tlen])
; - no need to keep looking if the longest thread you could find is too short!
num_saved_threads = 0
;a useful counter, primarily used to check if the first thread has been saved yet
raw_th_count = 0 ;counter for total number of possible threads (including rejected threads)
FOR j=0L, nt-(last_timestep+1) DO BEGIN ;Search over time
    FOR i=max_dist_jump, nx-(max_dist_jump+1) DO BEGIN ; Search over x

        IF image[i,j] GT mini THEN BEGIN
            ;Found the start of a possible thread

            ;reset values in temp structure
            temp_th.pos[0:-1] = 0.0
            temp_th.pos[0:j] = -1.0
            temp_th.err_pos[0:-1] = 0.0
            temp_th.err_pos[0:j] = -1.0
            temp_th.bin_flags[0:-1] = 0.0
            temp_th.bin_flags[0:j] = -1
            temp_th.start_bin = j
            temp_th.end_bin = j
            temp_th.length = 1

            ;Copy over the first data point
            ;(no need to clear values from the dummy array since, by definition,
            ; start points do not fall within the search boxes of any earlier
            ; data points)
            temp_th.pos[j] = image[i,j]
            temp_th.err_pos[j] = im_err[i,j]
            temp_th.bin_flags[j] = 2

            IF full_gauss EQ 1 THEN BEGIN
                temp_th.inten[0:-1] = 0.0
                temp_th.inten[0:j] = -1.0
                temp_th.err_inten[0:-1] = 0.0
                temp_th.err_inten[0:j] = -1.0
                temp_th.wid[0:-1] = 0.0
                temp_th.wid[0:j] = -1.0
                temp_th.err_wid[0:-1] = 0.0
                temp_th.err_wid[0:j] = -1.0

                temp_th.inten[j] = located.peaks[i,j, 1]
                temp_th.err_inten[j] = located.errs[i,j, 2]
                temp_th.wid[j] = located.peaks[i,j, 2]
                temp_th.err_wid[j] = located.errs[i,j, 2]
            ENDIF

            ;Select the initial search box:
            ;   width is equal to 2*max_dist_jump + 1 (in space)
            ;   height is equal to max_time_skip (in time)
            a = image[(i-max_dist_jump):(i+max_dist_jump),j+1:j+max_time_skip]
            a = a[reorder_cols, *]

            h = j ;current peak timestep
            k = i ;current peak distance

			WHILE MAX(a) GT mini DO BEGIN
                ;keep following the thread as long as there is a pixel with a
                ;value greater than the minimum in the current search box

                b = WHERE(a GT mini) ;returns -1 when it fails to find anything

                IF b[0] LT 0 THEN BEGIN
                    ;no peaks in box; set "a" to min value so the loop is broken
                    a = mini
                ENDIF ELSE BEGIN
                    ;find the coordinates (relative to the current peak) of the
                    ;nearest non-minimum point in the box
                    xm = dist_pxl_shift[b[0]]
                    ym = time_pxl_shift[b[0]]

                    k = k+xm ;update current distance
                    h = h+ym ;update current timestep

                    ;saves values at the new coords then erases them from the
                    ;dummy array so values are not used twice
                    temp_th.pos[h] = image[k,h]
                    temp_th.err_pos[h] = im_err[k,h]
                    temp_th.bin_flags[h] = 2
                    image[k,h] = mini
                    im_err[k,h] = mini

                    IF full_gauss EQ 1 THEN BEGIN
                        temp_th.inten[h] = located.peaks[k,h, 1]
                        temp_th.err_inten[h] = located.errs[k,h, 2]
                        temp_th.wid[h] = located.peaks[k,h, 2]
                        temp_th.err_wid[h] = located.errs[k,h, 2]
                    ENDIF

                    ;assign each pixel in a given thread the same ID number
                    IF KEYWORD_SET(debug) THEN th_debug.raw_th_num[k,h+ym] = raw_th_count

                    ;Find the indices of the next search box
                    left_ind = k - max_dist_jump
                    right_ind = k + max_dist_jump
                    sooner_ind = h + 1
                    later_ind = h + max_time_skip
                    IF later_ind GE (nt-1) THEN later_ind = nt-1

                    ;Update the search box to prepare for the next loop
                    IF (k LE 0) OR (k GE nx-1) OR (h GE nt-1) THEN BEGIN
                        ;quit the loop if the last selected peak is on the very
                        ;edge of the image (since more than half the next search
                        ;box will lie outside of the image).
                        a = mini
                    ENDIF ELSE BEGIN
                        IF left_ind LT 0 THEN BEGIN
                            ;adjust for missing cols on the left
                            first_real_col = ABS(left_ind)
                            a = FLTARR(search_box_width, later_ind-sooner_ind+1)
                            a[*,*] = mini
                            a[first_real_col:-1, *] = image[0:right_ind,sooner_ind:later_ind]
                            a = a[reorder_cols, *]
                        ENDIF ELSE IF right_ind GE nx THEN BEGIN
                            ;adjust for missing cols on the right
                            last_real_col = 2*max_dist_jump - (right_ind-nx+1)
                            a = FLTARR(search_box_width, later_ind-sooner_ind+1)
                            a[*,*] = mini
                            a[0:last_real_col, *] = image[left_ind:(nx-1),sooner_ind:later_ind]
                            a = a[reorder_cols, *]
                        ENDIF ELSE BEGIN
                            ;[Most common case]
                            ;Entire width of the search box is inside the image
                            a = image[left_ind:right_ind,sooner_ind:later_ind]
                            a = a[reorder_cols, *]
                        ENDELSE
                    ENDELSE
                    ; IF k LT max_dist_jump OR k GE nx-max_dist_jump THEN BEGIN
                    ;     ;stops loop from crashing as it reaches the side edges
                    ;     a = mini
                    ; ENDIF ELSE BEGIN
                    ;     IF h GT nt-(max_time_skip+2) THEN BEGIN
                    ;         a = mini
                    ;     ENDIF ELSE BEGIN
                    ;         ;update the search box
                    ;         a = image[(k-max_dist_jump):(k+max_dist_jump),h+1:h+max_time_skip]
                    ;         a = a[reorder_cols, *]
                    ;     ENDELSE
                    ; ENDELSE

                    temp_th.end_bin = h
                    temp_th.length = 1 + h - j

                ENDELSE

            ENDWHILE

            raw_th_count = raw_th_count + 1

            IF h LT (nt-1) THEN BEGIN
                ;fill remaining timesteps of output array that are NOT in the given thread
                temp_th.pos[h+1:-1] = -1
                temp_th.err_pos[h+1:-1] = -1
                temp_th.bin_flags[h+1:-1] = -1

                IF full_gauss EQ 1 THEN BEGIN
                    temp_th.inten[h+1:-1] = -1
                    temp_th.err_inten[h+1:-1] = -1
                    temp_th.wid[h+1:-1] = -1
                    temp_th.err_wid[h+1:-1] = -1
                ENDIF
            ENDIF

            ;@save_th###########################################################
            ;saving tracked thread to the output array (assuming its long enough)
            IF temp_th.length GE min_tlen THEN BEGIN
                num_nonzero_pts = N_ELEMENTS(where(temp_th.pos GT 0.0))
                num_zero_vals = N_ELEMENTS(where(temp_th.pos EQ 0.0, /NULL))

                ;skip threads with less than two real, positive values
                IF num_nonzero_pts GT 2.0 THEN BEGIN

                    ;skip threads where half the data points are set to zero,
                    ;i.e. no value was obtained at fitting stage (data gaps).
                    IF num_zero_vals LT 0.5*temp_th.length THEN BEGIN

                        IF num_saved_threads EQ 0 THEN BEGIN
                            threads = [temp_th]
                            num_saved_threads = 1
                        ENDIF ELSE BEGIN
                            threads = [temporary(threads), temp_th]
                            num_saved_threads = num_saved_threads + 1
                        ENDELSE
                    ENDIF
                ENDIF

                IF KEYWORD_SET(check) AND full_gauss EQ 1 THEN BEGIN
                    kk = num_saved_waves - 1 ;index of the current thread
                    !p.multi=[0,1,3]
                    window,1
                    plot,threads[kk].pos[*],xtitle='Time',ytitle='Central position',charsize=2
                    oploterror,threads[kk].pos[*],threads[kk].err_pos[*],psym=1
                    IF full_gauss EQ 1 THEN BEGIN
                        plot,threads[kk].inten[*],xtitle='Time',ytitle='Maximum Intensity',charsize=2
                        oploterror,threads[kk].inten[*],threads[kk].err_inten[*],psym=1
                        plot,threads[kk].wid[*],xtitle='Time',ytitle='Width',charsize=2
                        oploterror,threads[kk].wid[*],threads[kk].err_wid[*],psym=1
                        print,threads[kk].err_inten[*]
                    ENDIF
                    pause
                ENDIF
            ENDIF

        ENDIF

    ENDFOR
ENDFOR
END
