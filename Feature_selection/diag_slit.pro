;+
;NAME: DIAG_SLIT2
;
;PURPOSE:
;   Extracts a time-distance slice along a diagonal line in a series of images.
;   Can operate in an interactive mode. Basic code framework was adapted from
;   diag_slit.pro
;
;INPUTS:
;   in_data - 2D or 3D array containing image data
;
;OPTIONAL INPUTS:
;   index - index structure containing the .fits file header information
;           (currently only accepts SDO/AIA index structures). If provided,
;           many parameters (solar radius, center, image res., ect...) will be
;           extracted from the index structure.
;   x1,x2,y1,y2 - X & Y coordinates at either end of the diagonal line
;   res - spatial resolution of data (in [arcsec]). Default is 0
;   cad - temporal cadence of the data (in [s]). Default is 0
;   km_per_arcsec - ratio of [km]/[arcsec]. Defaults to 725.27 which is the
;                   mean scale of the solar surface.
;   /debug - if set, will overplot the data slit that that is to be extracted based
;            on only the start x value, total width, and calculated slope. Useful
;            for debugging issues.
;   fn - frame number to plot
;   /log_plot - if set, will log-scale the example image intensity values
;   xsize, ysize - size of graphics window (command for widgets)
;   win_id - set to value of open graphics window you want to use. If not set,
;            opens next free window.
;   /noopen - if set, suppresses the opening of plot windows
;
;OUTPUTS:
;   out_data - output array where data along the diagonal line is stored
;   meta_out - structure containing various metadata about the extracted slit
;   outvec - array of [x1,x2,y1,y2] parameters used for the circle
;
;HISTORY: Name---------Date---------Description
;         R Morton  ???, 2011  Created by Richard Morton in 2011
;         R Morton  ???, 2014  Added moving plot line
;         R Morton  NOV, 2014  added some graphics set up to use with widgets
;         R Morton  APR, 2015  fixed rounding error with fix command
;         R Morton  MAY, 2015  Introduced interpolation routine to extract data
;                              points - made diag_slit2 redundant!
;         M Weberg  MAR, 2017  Forked code and created a new "diag_slit2" program.
;                              Main differnce is the addition of metadata handling
;                              and output.
;
;TO DO:
;   - Test speed of procedure and optomise if needed
;   - Confirm that the data slit is ALWAYS extracted running FROM (x1,y1) TO (x2,y2)
;-

;helper function (taken from calc_mean_data_cadence.pro but renamed to avoid issues)
FUNCTION CALC_MEAN_INDEX_CAD, fits_index, trim=trim, time_diff=time_diff
;WARNING! currently only works if all of the data is from the same day
COMPILE_OPT IDL2

num_timesteps = N_ELEMENTS(fits_index.date_obs)

IF KEYWORD_SET(trim) THEN BEGIN
    trim_locs = STRSPLIT(fits_index[0].date_obs, '.') ;find where the "." is
    trim_fits_date = STRMID(fits_index.date_obs, 0, LONG(trim_locs[-1]-1)) ;cuts off the fractional seconds
ENDIF ELSE BEGIN
    trim_fits_date = fits_index.date_obs ;keeps the full timestamp
ENDELSE

split_date = STRSPLIT(trim_fits_date, 'T:', /EXTRACT) ; YYYY-MM-DD, Hr, Min, Sec

datetime_arr = split_date.ToArray() ;array of form [timestep_col, time_element_row]
time_in_sec = float(datetime_arr[*,1])*3600.0+float(datetime_arr[*,2])*60.0+float(datetime_arr[*,3])

sec_diff = time_in_sec - SHIFT(time_in_sec, 1)

mean_cadence = MEAN(sec_diff[1:-1]) ;ignore 0th element since SHIFT wraps around
time_diff = sec_diff

RETURN, mean_cadence
END

;@main##########################################################################
;###############################################################################
;###############################################################################
PRO DIAG_SLIT2, in_data, out_data, index=index, x1=x1, x2=x2, y1=y1, y2=y2, $
                res=res, cad=cad, km_per_arcsec=km_per_arcsec, $
                debug=debug, outvec=outvec, fn=fn, log_plot=log_plot, $
                xsize=xsize, ysize=ysize, win_id=win_id, noopen=noopen

COMPILE_OPT IDL2
;@int_meta######################################################################
;INITIALIZE METADATA STRUCTURE AND COMMON BLOCK
;###############################################################################
COMMON slit_meta_dat, slit_meta

slit_meta = {date_extracted:systime(), dataset:'unknown', wavelength:-1, $
             start_date:'unknown', end_date:'unknown', notes:'', $
             source_size:[0,0,0], res:0.0d, cad:0.0d, $
             km_per_arcsec:0.0d, ref_center:[0.0,0.0], r_sun:0.0, $
             slit_type:'diag', start_coord:[0.0,0.0], end_coord:[0.0,0.0], width:1}

;@set_units#####################################################################
;UNIT CONVERSIONS AND TAGGING
;###############################################################################
IF NOT KEYWORD_SET(res) THEN res = 0d
IF NOT KEYWORD_SET(cad) THEN cad = 0d
IF NOT KEYWORD_SET(km_per_arcsec) THEN km_per_arcsec = 725.27

;@check_inputs##################################################################
;CHECKING DIMENSIONS OF INPUT DATA
;###############################################################################
sz = size(in_data)
IF sz[0] EQ 3 THEN BEGIN
    nx_data = sz[1]
    ny_data = sz[2]
    nt_data = sz[3]
ENDIF ELSE BEGIN
    IF sz[0] EQ 2 THEN BEGIN
        PRINT, 'WARNING: input data is only 2D! Output array will be 1D.'
        nx_data = sz[1]
        ny_data = sz[2]
        nt_data = 1
        fn = 0
    ENDIF ELSE BEGIN
        MESSAGE, 'ERROR: please input a valid 2D or 3D array'
    ENDELSE
ENDELSE

;@plot_ex########################################################################
;PLOT EXAMPLE IMAGE FOR REFERENCE AND/OR INTERACTIVE SELECTION
;###############################################################################
IF NOT keyword_set(noopen) THEN BEGIN
    IF n_elements(fn) EQ 0 THEN fn = 0
    IF fn GT (nt_data-1) THEN fn = nt_data-1

    ex_img = in_data[*,*,fn]
    IF KEYWORD_SET(log_plot) THEN BEGIN
        min_ex = MIN(ex_img)
        ex_img = ex_img + (min_ex + 1)
        ex_img = ALOG10(ex_img)
    ENDIF

    ;Set plotting elements
    ;writing to buffer pixmap reduces on-screen jitter!
    IF n_elements(win_id) EQ 0 THEN BEGIN
        determine_window_size, xsize=xsize, ysize=ysize, openwin=openwin, /new_win
        window, openwin[0], xs=xsize, ys=ysize
        window_var = openwin[0]
    ENDIF ELSE BEGIN
        wset, win_id
        window_var = win_id
    ENDELSE
    window, xs=xsize, ys=ysize, /pixmap, /free
    pixID = !d.window
    tvim, ex_img
    wset, window_var
    device, copy=[0,0,!d.x_size,!d.y_size,0,0,pixid]
ENDIF

;@get_meta######################################################################
;IF GIVEN AN INDEX STRUCTURE WITH IMAGE METADATA, LOOKUP VALUES AS NEEDED
;###############################################################################
IF KEYWORD_SET(index) THEN BEGIN
    index_tag_names = TAG_NAMES(index[0])
    test_for_telescope_tag = WHERE(STRUPCASE(index_tag_names) EQ 'TELESCOP')
    IF test_for_telescope_tag GE 0 THEN BEGIN
        IF index[0].TELESCOP.StartsWith('SDO', /FOLD_CASE) THEN BEGIN
            ;SDO / AIA data dectected!
            slit_meta.dataset = index[0].TELESCOP+'_'+index[0].INSTRUME+'_lvl_'+string(index[0].LVL_NUM, format='(f4.2)')
            slit_meta.wavelength = index[0].WAVELNTH
            slit_meta.start_date = index[0].DATE_OBS
            slit_meta.end_date = index[-1].DATE_OBS
            res = index[0].CDELT1
            cad = CALC_MEAN_INDEX_CAD(index)
            km_per_arcsec = (index[0].RSUN_REF/index[0].RSUN_OBS)/1000.0
            ;note: the AIA FITS coordinates system labels the lower left corner
            ;of a full-disk image as [1,1] not [0,0] as in IDL (NO SHIFT INCLUDED!)
            slit_meta.ref_center = [(index[0].CRPIX1), (index[0].CRPIX2)]
            slit_meta.r_sun = index[0].R_SUN
        ENDIF
    ENDIF
ENDIF

;@query#########################################################################
;QUERY USER FOR SLIT START AND END COORDINATES AS NEEDED
;###############################################################################
IF N_ELEMENTS(x1) EQ 0 OR N_ELEMENTS(y1) EQ 0 THEN BEGIN
    IF NOT KEYWORD_SET(noopen) THEN BEGIN
        PRINT, 'Pleast select the slit start point with the cursor'
        CURSOR, x1, y1, /Data, /Down
        PRINT, 'x1='+STRTRIM(FIX(x1+0.5), 2), ' y1='+STRTRIM(FIX(y1+0.5), 2)
    ENDIF ELSE BEGIN
        READ, x1, PROMPT='Please enter x1: '
        READ, y1, PROMPT='Please enter y1: '
    ENDELSE

    ;make the inputs floating point variables (but with whole pixel values)
    x1 = 1.0*FIX(x1+0.5)
    y1 = 1.0*FIX(y1+0.5)
ENDIF ELSE BEGIN
    x1 = 1.0*x1
    y1 = 1.0*y1
ENDELSE

IF N_ELEMENTS(x2) EQ 0 OR N_ELEMENTS(y2) EQ 0 THEN BEGIN
    IF NOT KEYWORD_SET(noopen) THEN BEGIN
        PRINT,'Please select the slit end point with the cursor'
        IF !MOUSE.button EQ 1 THEN BEGIN
            CURSOR, do_nothing_x, do_nothing_y, /UP, /WAIT
            !MOUSE.button = 0
        ENDIF
        WHILE (!MOUSE.button NE 1) DO BEGIN
            CURSOR, x2, y2, /change
            window, xs=xsize, ys=ysize, /pixmap, /free
            pixID = !d.window
            tvim, ex_img
            WSET, window_var
            DEVICE, copy=[0,0,!d.x_size,!d.y_size,0,0,pixid]
            PLOTS, [x1,x2], [y1,y2]
        ENDWHILE
        PRINT, 'x2='+STRTRIM(FIX(x2+0.5), 2), ' y2='+STRTRIM(FIX(y2+0.5), 2)
    ENDIF ELSE BEGIN
        READ, x2, PROMPT='Please enter x2: '
        READ, y2, PROMPT='Please enter y2: '
    ENDELSE

    ;make the inputs floating point variables (but with whole pixel values)
    x2 = 1.0*FIX(x2+0.5)
    y2 = 1.0*FIX(y2+0.5)

ENDIF ELSE BEGIN
    x2 = 1.0*x2
    y2 = 1.0*y2
ENDELSE

;@validate######################################################################
;VALIDATING SLIT VALUES AND RESETTING IF NEEDED (prevents unexpected behavior)
;###############################################################################
IF x1 LT 0 THEN x1 = 0.0
IF x1 GE (nx_data) THEN x1 = 1.0*nx_data - 1
IF y1 LT 0 THEN y1 = 0.0
IF y1 GE (ny_data) THEN y1 = 1.0*ny_data - 1
IF x2 LT 0 THEN x2 = 0.0
IF x2 GE (nx_data) THEN x2 = 1.0*nx_data - 1
IF y2 LT 0 THEN y2 = 0.0
IF y2 GE (ny_data) THEN y2 = 1.0*ny_data - 1

IF NOT KEYWORD_SET(noopen) THEN PLOTS, [x1,x2], [y1,y2] ;plot validated slit

;@get_slit######################################################################
;EXTRACT THE DATA SLIT
;###############################################################################
x_diff = ABS(x2-x1)
y_diff = ABS(y2-y1)

slit_lnth = ROUND(SQRT(y_diff^2 + x_diff^2))
out_data = FLTARR(slit_lnth,nt_data)

mindat = MIN(in_data)
IF mindat GT 0 THEN mindat = -1.0 ELSE mindat = mindat-  0.1*STDDEV(in_data)

IF (x_diff GT 0) AND (y_diff GT 0) THEN BEGIN
    ;True diagonal slits

    ; Calculating gradient and intercept
    m = (y2-y1)/((x2-x1))
    c = -x1*(y2-y1)/(x2-x1)+y1

    IF KEYWORD_SET(debug) THEN BEGIN
      IF x1 LT x2 THEN BEGIN
         x = x1 + FINDGEN(x_diff) ;note the sign differnces (+)
         OPLOT, x, m*x+c, psym=1, nsum=3
      ENDIF ELSE BEGIN
         x = x1 - FINDGEN(x_diff) ;note the sign differnces (-)
         OPLOT, x, m*x+c, psym=1, nsum=3
      ENDELSE
    ENDIF

    IF x_diff GT y_diff THEN BEGIN
        FOR i=0, (nt_data-1) DO BEGIN
            x = MIN([x1,x2]) + FINDGEN(x_diff)
            y = m*x+c
            x = CONGRID(TEMPORARY(x), slit_lnth,cubic=-0.5)
            y = CONGRID(TEMPORARY(y), slit_lnth,cubic=-0.5)
            IF y1 GT y2 THEN BEGIN
                out_data[*,i] = REVERSE(INTERPOLATE(REFORM(in_data[*,*,i]), x, y, cubic=-0.5, missing=mindat))
            ENDIF ELSE BEGIN
                out_data[*,i] = INTERPOLATE(REFORM(in_data[*,*,i]), x, y, cubic=-0.5, missing=mindat)
            ENDELSE
        ENDFOR
    ENDIF

    IF y_diff GE x_diff THEN BEGIN
        FOR i=0, (nt_data-1) DO BEGIN
            y = MIN([y1,y2]) + FINDGEN(y_diff)
            x = (y-c)/m
            x = CONGRID(TEMPORARY(x), slit_lnth, cubic=-0.5)
            y = CONGRID(TEMPORARY(y), slit_lnth, cubic=-0.5)
            IF x1 GT x2 THEN BEGIN
                out_data[*,i] = REVERSE(INTERPOLATE(REFORM(in_data[*,*,i]), x, y, cubic=-0.5, missing=mindat))
            ENDIF ELSE BEGIN
                out_data[*,i] = INTERPOLATE(REFORM(in_data[*,*,i]), x, y, cubic=-0.5, missing=mindat)
            ENDELSE
        ENDFOR
    ENDIF

ENDIF ELSE BEGIN

    IF x_diff EQ 0 THEN BEGIN
        ;Vertical slits
        ;Note: currently fixed to extract the data from the bottom to the top
        dum = [y1,y2]
        y1 = min(dum)
        y2 = max(dum)
        out_data = fltarr(y_diff, nt_data)
        out_data[*,*] = reform(in_data[x1, y1:y2-1, *])
    ENDIF

    IF y_diff EQ 0 THEN BEGIN
        ;Horizontal slits
        ;Note: currently fixed to extract data from left to right
        dum = [x1,x2]
        x1 = min(dum)
        x2 = max(dum)
        out_data = fltarr(x_diff, nt_data)
        out_data[*,*] = reform(in_data[x1:x2-1, y1, *])
    ENDIF
ENDELSE

;@save_meta#####################################################################
;SAVING INFORMATION TO THE SLIT_META STRUCTURE (AS WELL AS OPTIONAL OUTPUTS)
;###############################################################################
slit_meta.source_size = [nx_data, ny_data, nt_data]
slit_meta.start_coord = [x1,y1]
slit_meta.end_coord = [x2,y2]
slit_meta.res = res
slit_meta.cad = cad
slit_meta.km_per_arcsec = km_per_arcsec

;Output optional vector
outvec = [x1,x2,y1,y2]

END
