;+
;NAME: LOCATE_THINGS
;
;PURPOSE:
;   Locates peaks in time-distance plot with sub-pixel accuracy
;
;PROCEDURE OUTLINE:
;   Takes an input time-distance diagram (can be optionally unsharp masked and/or
;   high-pass filtered) and locates the peaks using a find_max crawling routine.
;   The maximum then has to have a specific gradient (default > 0.5) to be classed
;   as a peak. Once a maximum is found a Gaussian fit (with linear terms) is
;   applied to the surrouding pixels to provide subpixel position for the maximum.
;   Alternatively, the code can be set to simply return the nearest whole-pixel
;   peak location.
;
;EXTERNAL CALLS: mpfitfun, mygauss_plus_linear, fifteenb (only with /check)
;
;INPUTS:
;   data - time-distance diagram
;   errors - errors on intensity for each pixel, i.e., estimates for Photon
;            noise, array should be same size as data array. Supplied to
;            Gaussian fitting routine. If not set, default is 1
;   grad - limits on gradient (default = 0.5) - Key parameter to modify!
;
;OPTIONAL INPUTS:
;   res - spatial resolution of data in [arcsec]. Defaults to a value of 1 [pixel]
;   cad - temporal cadence of the data in [s]. Defaults to a value of 1 [timestep]
;   km_per_arcsec - ratio of [km]/[arcsec]. Defaults to 725.27 which is the
;                   mean scale of the solar surface.
;   set_dx_units - units of distance that result from multiplying pixel
;                  locations by the output value "dx". Defaults to units of
;                  [arcsec] if a "res" value is given or [pixels] if not.
;                  May choose from "pixels", "arcsec", "m", "km", or "Mm".
;   set_dt_units - units of time that results from multiplying timestep number
;                  output value "dt". Defaults to units of [min] if a "cad"
;                  value is given or [timesteps] if not. May choose from
;                  "timesteps", "s", "min", or "hr".
;   /invert - If set, will invert the intensity values of the input data before
;             finding peak locations.This results in the code finding the local
;             MINIMA in the original image.
;   bad_data - value used in the input data to indicate bad data points. When
;              when inverting the image, all values <= bad_data will be preserved
;              and NOT have thier values changed. This prevents invalid data
;              flagged with negative intensity values from being selected as peaks.
;              Only used when /invert is also set. Default is -999
;   /despike - if set, will calculate the mean and stddev of the intensity values
;              in each timestep. All values with an intensity larger than a
;              selected number of stddev above the mean are then replaced by the
;              average of the two adjacent data points
;   spike_sigma - number of stddev above the mean to use for filtering out large
;                 data spikes. Only used when /despike is also set. Default is 3.0
;   num_search_bins - number of bins to look ON EITHER SIDE of a potential peak
;                     for determining if it is a local maxima. Elsewhere known
;                     as the 'order' of the search algorithm. Default is 5. Note,
;                     this will also set the minimum distance allowed between peaks.
;   percent_grad - if set to a value between 0 and 100, will automatically rescale
;                  the intensity gradient to reject the specified percent of found
;                  peaks.
;   /nearest_pixel - if set, will ALWAYS return the nearest whole-pixel location
;                    of the peaks without bothering with any gaussian fitting
;   meas_size - sets TOTAL number of pixels used in the Gaussian fit.
;               Even inputs will have their value increased by one (so as to be
;               symmetric about the peak value). Default value is 7.
;   /simp_grad - set to use an analytic least-sqaures estimate for gradient.
;                This may be faster for very large input arrays!
;   /weighted_mean - [EXPERIMENTAL!] if set, will use a weighted mean to find the
;                    sub-pixel location of each peak. This may be considerably
;                    faster than gaussian fitting but initial tests do not show
;                    much difference from simply picking the nearest pixel.
;   cut_chisq - cut off value for unacceptable chi squared. Default set at 3
;               sigma confidence level (based on number of points used in fit)
;   shift_cut - cut off value for maximum allowable peak shift. If the fit
;               location shifts more than this, will instead use the peak
;               location with an error of 0.5 pixels. default is 1.5
;   /full_gauss - set to output the 'full_gaussian' results (i.e. save the
;                 widths). Replicates the old output of 'locate_things_fg.pro'
;   /check - if set, will plot various variables about each good peak to a
;            window. NOT recommended for long periods of observations.
;   /debug - if set, will save sub-pixel test parameters in the 'loc_debug_dat'
;            COMMON block for later inspection
;
;OUTPUTS:
;   grad_out - actual intnesity graident used to filter peaks (may be auto-scaled)
;   located - struture containing result arrays that is then saved to the
;             "located_dat" COMMON block. Format is as follows:
;       .peaks - [nx,nt,2] or [nx,nt,6] array with sub-pixel locations and
;                maximum values at each located peak (in [pixels])
;              - [*,*,0] - sub-pixel location of peak
;              - [*,*,1] - maximum intensity value at peak
;              - [*,*,2] - gaussian width (/full_gauss only)
;              - [*,*,3] - constant coeff (/full_gauss only)
;              - [*,*,4] - slope of linear term (/full_gauss only)
;              - [*,*,5] - i index of maximum data point (/full_gauss only)
;       .errs - [nx,nt] or [nx,nt,6] array with errors on the fit values.
;               Format is similar to .peaks (but for the errors instead)
;       .allpeaks - [nx,nt] array containing integer flags for ALL local maxima.
;                   The values are as folows:
;                       1 - local peak rejected by choosen gradient
;                       2 - gaussian fit failed. Defaulted to whole pixels
;                       3 - gaussian fit obtained for subpixel resolution
;       .grad_left - [nx,nt] array with the gradient values on the
;                            lefthand side of ALL potential peaks
;       .grad_right - same as as the above but for right-hand side gradients
;       .td_img - td-diagram AFTER any inversion or despiking operations are applied
;       .inverted - binary flag indicating if the image was inverted
;       .despiked - binary flag indicating id the image was despiked
;       .dx & .units_dx - [pixels] to [.units_dx] conversion factor
;       .dt & .units_dt - [timesteps] to [.units_dt] conversion factor
;       .res, .cad, & .km_per_arcsec - source data resolution and cadence info
;                                      See descriptions above for more details
;
;HISTORY: Name---------Date---------Description
;         R Morton  OCT, 2012 - Initial coding
;         R Morton  NOV, 2014 - super update! Added structure format to remove all
;                               the arrays. Also added COMMON variables so re-used
;                               values/structures are passed automatically
;         R Morton  MAR, 2015 - fixed bug with definition of initial minimum value
;                               to set all located.peaks values to
;         R Morton  MAR, 2016 - updated so similar to locate_things_fg
;         R Morton  14 MAR, 2016 - Release as version 1.0 of NUWT
;         R Morton  MAY, 2016 - Added option for analytic LS calc for graident in
;                               hope of speed increase.
;                               Also removed needless for loop at bottom!!
;                               Added calcuation of peak intensity
;                               Added better estimate for Chi^2 cutoff
;         M Weberg  MAY, 2016 - Implemented faster inital peak finding. Still requires
;                               looping over the potential peaks for computing gradients
;         M Weberg  JUL, 2016 - Added arrays to output stucture for saving the left- and
;                               right-hand gradient values.
;         M Weberg  JUL, 2016 - IMPORTANT! merged the program 'locate_things.pro' with
;                               'locate_things_fg.pro' and added the /full_gauss keyword
;         M Weberg  JAN, 2017 - merged in 'debug' option and experimental 'weighted_mean'
;                               method from alternative versions of the code
;         M Weberg  NOV, 2017 - added the '/invert' keyword for finding local minima.
;                               Also added the 'bad_data', '/despike', and 'spike_sigma'
;                               keywords and included more metadata in the output structure
;         M Weberg  DEC, 2017 - finally merged in the program 'locate_things_min.pro'
;                               and added the '/nearest_pixel' option. This is now
;                               the only version of locate_things that needs to
;                               be maintained.
;         M Weberg  FEB, 2018 - Added the percent_grad mode. This included splitting
;                               the gradient calcuation and peak filtering into
;                               two seperate step. this enables new modes to be
;                               easily added in the future.
;                  27 MAR 2018 - Release as version 2.0 of NUWT
;
;TO DO:
;   - add better handling of flat-topped local maxima. Currently, if two pixels
;     have the same max value, the program will select the left-most as the peak
;   - Still need to work out the best way to describe intensity values for gaussian fit!
;-

PRO LOCATE_THINGS, input_data, errors=errors, $
                        res=res, cad=cad, km_per_arcsec=km_per_arcsec, $
                        set_dx_units=set_dx_units, set_dt_units=set_dt_units, $

                        invert=invert, bad_data=bad_data, $
                        despike=despike, spike_sigma=spike_sigma, $

                        num_search_bins=num_search_bins, grad=grad, $
                        percent_grad=percent_grad, out_grad=out_grad, $
                        nearest_pixel=nearest_pixel, meas_size=meas_size,$
                        simp_grad=simp_grad, weighted_mean=weighted_mean, $

                        cut_chisq=cut_chisq, shift_cut=shift_cut, $
                        full_gauss=full_gauss, check=check, debug=debug

COMPILE_OPT IDL2
;@defaults######################################################################
;SETTING DEFAULT VALUES
;###############################################################################
COMMON located_dat, located

IF KEYWORD_SET(invert) THEN invert = 1 ELSE invert = 0


;Limit value used in the input array to indicate bad data (only used if /invert is set)
IF N_ELEMENTS(bad_data) EQ 0 THEN bad_data = -999

IF KEYWORD_SET(despike) THEN despike = 1 ELSE despike = 0
;Sets the number of standard deviations above the time-slice mean used to dectected
;   large data spikes (only used if /despike it set)
IF N_ELEMENTS(spike_sigma) EQ 0 THEN spike_sigma = 3.0

;Sets the number of bins to search on either side of a potential peak
IF NOT KEYWORD_SET(num_search_bins) THEN num_search_bins = 5

;Sets default gradient for filtering peaks
IF N_ELEMENTS(grad) EQ 0 THEN grad = 0.5

;Sets total number of pixels to use for the gaussian fit
;Note 1: a default value of 7 was also used for measuring fibrils in ROSA data
;Note 2: 5 data points were previosly fit by Thurgood et al. HOWEVER this unfortunatly
;        will not work well with our chosen "guassian plus linear" fit function which
;        attempts to find 5 different fit parameters (therby reducing the degrees
;        of freedom to zero if we use just 5 bins).
IF NOT KEYWORD_SET(meas_size) THEN meas_size = 7
IF meas_size MOD 2 GT 0 THEN meas_size = meas_size + 1 ;ensures num of total pxls is odd
meas_size = fix(meas_size)
half_meas = floor(meas_size/2.0)

;Set cut off value for unacceptable chi squared
;[DEFAULT] Set at value of 3 sigma confidence level
;[2017-Feb] Note: tests with simulated data yield better results when we effectively
;   "remove" the chisq cutoff by setting it to a very large number (such as 1e12).
IF NOT KEYWORD_SET(cut_chisq) THEN cut_chisq = chisqr_cvf(1.0 - 0.9973, meas_size-5)

;Set cut off value for unacceptable peak shift
IF NOT KEYWORD_SET(shift_cut) THEN shift_cut = 1.5

;Find size of the input data
sz = size(input_data)
nx = sz[1]
nt = sz[2]

;@units#########################################################################
;UNIT CONVERSIONS AND LABELS
;###############################################################################
;Set defaults and validate input units (if given)
IF NOT KEYWORD_SET(km_per_arcsec) THEN km_per_arcsec = 725.27
IF NOT KEYWORD_SET(set_dx_units) THEN set_dx_units = 'arcsec'
IF NOT KEYWORD_SET(set_dt_units) THEN set_dt_units = 'min'
IF TOTAL(STRMATCH(['pixels', 'arcsec', 'm', 'km', 'Mm'], set_dx_units)) EQ 0 THEN BEGIN
    print, 'WARNING: invalid value for "set_dx_units"! Using default units of "pixels" or "arcsec" (as appropriate)'
    set_dx_units = 'arcsec'
ENDIF
IF TOTAL(STRMATCH(['timesteps', 's', 'min', 'hr'], set_dt_units)) EQ 0 THEN BEGIN
    print, 'WARNING: invalid values for "set_dt_units"! Using default units of "timesteps" or "min" (as appropriate)'
    set_dt_units = 'min'
ENDIF

IF KEYWORD_SET(res) THEN BEGIN
    ;if res is known, PROHIBIT setting the dx units to 'pixels'
    IF set_dx_units EQ 'pixels' THEN BEGIN
        print, 'WARNING: set_dx_units =/= "pixels" when "res" is known. Defaulting to "arcsec"'
        set_dx_units = 'arcsec'
    ENDIF
ENDIF ELSE BEGIN
    ;[DEFAULT] res is unknown / not given
    res = 1.0
    set_dx_units = 'pixels'
ENDELSE

IF KEYWORD_SET(cad) THEN BEGIN
    ;if cad is known, PROHIBIT setting the dt units to 'timesteps'
    IF set_dt_units EQ 'timesteps' THEN BEGIN
        print, 'WARNING: set_dt_units =/= "timesteps" when "cad" is known. Defaulting to "min"'
        set_dt_units = 'min'
    ENDIF
ENDIF ELSE BEGIN
    ;[DEFAULT] cad is unknown / not given
    cad = 1.0
    set_dt_units = 'timesteps'
ENDELSE

;distance units conversions
IF set_dx_units EQ 'pixels' THEN dx = 1.0
IF set_dx_units EQ 'arcsec' THEN dx = res
IF set_dx_units EQ 'm' THEN dx = res*km_per_arcsec*1000.0
IF set_dx_units EQ 'km' THEN dx = res*km_per_arcsec
IF set_dx_units EQ 'Mm' THEN dx = res*km_per_arcsec/1000.0

;time unit conversions
IF set_dt_units EQ 'timesteps' THEN dt = 1.0
IF set_dt_units EQ 's' THEN dt = cad
IF set_dt_units EQ 'min' THEN dt = cad/60.0
IF set_dt_units EQ 'hr' THEN dt = cad/3600.0

;@invert_img####################################################################
;INVERT THE IMAGE IF LOOKING FOR LOCAL MINIMA RATHER THAN LOCAL MAXIMA
;###############################################################################
IF KEYWORD_SET(invert) THEN BEGIN
    ;If set, will find the local MINIMA in the input image
    max_intensity_val = MAX(input_data, /NAN)
    min_intensity_val = MIN(input_data, /NAN)
    IF min_intensity_val LE bad_data THEN BEGIN
        ;preserves bad data values equal to or less than a given limit [default is -999]
        loc_neg_vals = WHERE(input_data LE bad_data)
        td_diagram = -input_data[*,*] + max_intensity_val
        td_diagram[loc_neg_vals] = input_data[loc_neg_vals]
    ENDIF ELSE BEGIN
        ;simple inversion of the data
        td_diagram = -input_data[*,*] + max_intensity_val
    ENDELSE
ENDIF ELSE BEGIN
    ;[DEFAULT] Finds the local MAXIMA in the input image
    td_diagram = input_data[*,*]
ENDELSE

;@despike#######################################################################
;REMOVES DATA SPIKES ABOVE A CERTAIN THRESHOLD
;###############################################################################
IF KEYWORD_SET(despike) THEN BEGIN
    ;If set, will remove large data values
    ts_mean = MEAN(td_diagram, dimension=1, /NAN)
    ts_stddev = STDDEV(td_diagram, dimension=1, /NAN)
    spike_limit = ts_mean + spike_sigma*ts_stddev ;default limit is 3.0 sigma
    spike_limit = REBIN(REFORM(spike_limit, 1, nt), nx, nt) ;recast as a 2D array
    loc_spikes = WHERE(td_diagram GT spike_limit, /NULL)
    IF N_ELEMENTS(loc_spikes) GT 0 THEN BEGIN
        ;replace data spikes with the average of the two neighboring pixels
        right_vals = SHIFT(td_diagram, -1, 0)
        left_vals = SHIFT(td_diagram, 1, 0)
        td_diagram[loc_spikes] = (left_vals[loc_spikes] + right_vals[loc_spikes])/2.0
    ENDIF
    right_vals = 0 ;free up memory
    left_vals = 0
ENDIF

;@int_output####################################################################
;INITIALIZE THE "OUTPUT" STRUCTURES AND FILLING WITH DEFAULT VALUES
;###############################################################################
IF KEYWORD_SET(full_gauss) THEN BEGIN
    print, 'Returning full Gaussian results'
    nearest_pixel = 0
    ;'full gaussian' output including the fitted widths
    ; Replicates the output of a previously seperate code, 'locate_things_fg.pro'
    located = {peaks:fltarr(nx,nt,6), errs:fltarr(nx,nt,6), allpeaks:intarr(nx,nt), $
               grad_left:fltarr(nx, nt), grad_right:fltarr(nx,nt), $
               td_img:fltarr(nx,nt), inverted:invert, despiked:despike, $
               dx:dx, units_dx:set_dx_units, $
               dt:dt, units_dt:set_dt_units, $
               res:res, cad:cad, km_per_arcsec:km_per_arcsec}
    end_pk_ind = 5
    end_err_ind = 5
    nearest_end_pk_ind = 2 ;default for nearest whole-pixel results
    nearest_end_err_ind = 2
ENDIF ELSE BEGIN
    ;[DEFAULT] basic output
    located = {peaks:fltarr(nx,nt,2), errs:fltarr(nx,nt), allpeaks:intarr(nx,nt), $
               grad_left:fltarr(nx, nt), grad_right:fltarr(nx,nt), $
               td_img:fltarr(nx,nt), inverted:invert, despiked:despike, $
               dx:dx, units_dx:set_dx_units, $
               dt:dt, units_dt:set_dt_units, $
               res:res, cad:cad, km_per_arcsec:km_per_arcsec}
    end_pk_ind = 1
    end_err_ind = 0
    nearest_end_pk_ind = 1 ;default for nearest whole-pixel results
    nearest_end_err_ind = 0
ENDELSE

located.td_img = td_diagram ;useful for debugging or plotting

;Filling the 'peaks' array with invalid negative values
;doing this should avoid confusion with real results (which in theory could be negitive valued)
located.peaks[0:nx-1,0:nt-1,0:end_pk_ind] = (min(td_diagram)-10.0) < (-10)

mini = min(located.peaks)

;IF NO ERRORS GIVEN USE WEIGHTS EQUAL TO 1
IF N_ELEMENTS(errors) LE 1 THEN BEGIN
   errors = fltarr(nx,nt)
   errors[0:nx-1,0:nt-1] = 1.0
ENDIF

IF KEYWORD_SET(debug) THEN BEGIN
    COMMON loc_debug_dat, loc_debug
    ;Array for debugging
    loc_debug = {shift:FLTARR(nx,nt), sigma:FLTARR(nx,nt), width:FLTARR(nx,nt), $
                 chisq:FLTARR(nx,nt), coeff:FLTARR(nx, nt), $
                 cut_shift:shift_cut, cut_sigma:1.5, cut_width:meas_size, $
                 cut_chisq:cut_chisq, cut_coeff:mini}
ENDIF

;@find_peaks####################################################################
;FIND ALL POSSIBLE LOCAL PEAKS
;###############################################################################
;New, faster loop for finding all peaks at once
allpeaks = intarr(nx,nt)
allpeaks[*,*] = 1 ;all pixels are potentially peaks to start
FOR bin=1, num_search_bins DO BEGIN
    ;Finds peaks but iteratively filtering out pixels that are NOT local maxima
    ;note: if two pixels have the same value (i.e. flat-topped), this algorithm
    ;will select the left-most pixel as the local max (IMPROVE IN THE FUTURE?)
    compare_right = td_diagram - SHIFT(td_diagram, -bin, 0)
    compare_left = td_diagram - SHIFT(td_diagram, bin, 0)
    not_pk_locs = WHERE(allpeaks GT 0 and (compare_left LE 0 or compare_right LT 0), /NULL)
    allpeaks[not_pk_locs] = 0 ;removes incorrectly flagged peaks
ENDFOR

compare_right = 0 ;free up memory
compare_left = 0

;Cleans out false peaks too close to the edges of the array
allpeaks[0:num_search_bins, *] = 0
allpeaks[-num_search_bins:-1, *] = 0

;Finds the i,j indices of each potential peak
pk_locs = WHERE(allpeaks GT 0, /NULL)
num_peaks_to_test = N_ELEMENTS(pk_locs)
pk_indices = ARRAY_INDICES(allpeaks, pk_locs)

;@filter_peaks##################################################################
;FILTER PEAKS BASED ON GRADIENT AND PIXEL SHIFT
;###############################################################################
;Loops over each potential peak and tests for the required gradient
grad_right = FLTARR(nx,nt)
grad_left = FLTARR(nx,nt)
ind_grad = [-2,-1,0,1,2]
ind_fit = -half_meas + INDGEN(meas_size)

; First, find the gradients
j = -1
FOR p=0L, num_peaks_to_test-1 DO BEGIN
    i = pk_indices[0, p]
    IF pk_indices[1, p] GT j THEN BEGIN
        j = pk_indices[1, p]
        img_slice = REFORM(td_diagram[0:nx-1,j,0])
        err_slice = REFORM(errors[0:nx-1,j,0])
    ENDIF

    ;Finds gradients on either side of the maximum
    IF NOT KEYWORD_SET(simp_grad) THEN BEGIN
        ;[DEFAULT] fit lines and find slope on either side
        left_line = POLY_FIT(ind_grad,img_slice[i-4:i],1,yfit=fit,measure_errors=err_slice[i-4:i])
        right_line = POLY_FIT(ind_grad,img_slice[i:i+4],1,yfit=fit2,measure_errors=err_slice[i:i+4])
        grad_left[i,j] = left_line[1]
        grad_right[i,j] = right_line[1]
    ENDIF ELSE BEGIN
        ;Analytic gradient (faster)
        grad_left[i,j] = TOTAL(ind_grad*(img_slice[i-4:i]-TOTAL(img_slice[i-4:i])/5.0))/10.0
        grad_right[i,j] = TOTAL(ind_grad*(img_slice[i:i+4]-TOTAL(img_slice[i:i+4])/5.0))/10.0
    ENDELSE
ENDFOR

; Secondly, if a percent grad is being used, auto-scale the gradient threashold
IF N_ELEMENTS(percent_grad) GT 0 THEN BEGIN
    IF percent_grad LT 0 THEN percent_grad = 0.0
    IF percent_grad GT 100 THEN percent_grad = 100.0

    grad_left_arr = grad_left[pk_locs]
    grad_right_arr = grad_right[pk_locs]

    abs_max_left = abs(max(grad_left_arr))
    abs_min_right = abs(min(grad_right_arr))
    abs_largest_grad = MAX([abs_max_left, abs_min_right])

    num_test_grads = FIX(FLOAT(abs_largest_grad)/0.01)
    test_grad_arr = FINDGEN(num_test_grads)*0.01

    temp_percent_rejected = 0.0
    pi = 1
    WHILE temp_percent_rejected LE percent_grad DO BEGIN
        pi = pi + 1 ;never bother with pi = 0 since that is the zero grad case
        temp_loc_rejected = WHERE(grad_left_arr LE test_grad_arr[pi] OR grad_right_arr GE -1.0*test_grad_arr[pi])
        temp_percent_rejected = 100*FLOAT(N_ELEMENTS(temp_loc_rejected))/FLOAT(num_peaks_to_test)
    ENDWHILE

    grad = test_grad_arr[pi-1]
    print, '### A %grad of', percent_grad, ' corresponds to a intensity grad of ', grad
ENDIF

; Thirdly, filter by grad and then find sub-pixel location using gaussian fit
j = -1
FOR p=0L, num_peaks_to_test-1 DO BEGIN
    i = pk_indices[0, p]
    IF pk_indices[1, p] GT j THEN BEGIN
        j = pk_indices[1, p]
        img_slice = REFORM(td_diagram[0:nx-1,j,0])
        err_slice = REFORM(errors[0:nx-1,j,0])
    ENDIF

    ;If gradients greater than a certain value begin
    ;Gradient of quadratic evaluated at x=0
    IF (grad_left[i,j] GT grad) AND (grad_right[i,j] LT (-1.0)*grad) THEN BEGIN

        IF KEYWORD_SET(nearest_pixel) THEN BEGIN
            ;replicates the results of "locate_things_min"
            ;creating fake fit results (for quick and dirty saving of values)
            chisq = -1
            chisq_red = -1
            coeff = [mini, 0.0, 0.0, 0.0, 0.0]
            sigma = [0.0, 0.0, 0.0, 0.0, 0.0]
        ENDIF ELSE IF KEYWORD_SET(weighted_mean) THEN BEGIN
            ;[EXPERIMENTAL!] find sub-pixel location using weighted mean method
            sum_itensities = FLOAT(TOTAL(img_slice[i-half_meas:i+half_meas]))
            weights = img_slice[i-half_meas:i+half_meas]/sum_itensities
            wg_mean = TOTAL(weights*ind_fit)/TOTAL(weights)
            wg_stddev = TOTAL(weights*(ind_fit-wg_mean)^2)/TOTAL(weights)
            ;creating fake fit results (for quick and dirty saving of values)
            chisq = 0.0
            chisq_red = 0.0
            coeff = [1.0, wg_mean, wg_stddev, 0.0, 0.0]
            sigma = [1.0, 1.0, 1.0, 0.0, 0.0]
        ENDIF ELSE BEGIN
            ;[DEFAULT] find gaussian fit to surrounding points
            estimates = [img_slice[i], 0.0, 2.0, min(img_slice[(i-half_meas):(i+half_meas)]), 0.1]
            coeff = mpfitfun('mygauss_plus_linear', ind_fit, $
                             img_slice[(i-half_meas):(i+half_meas)], err_slice[(i-half_meas):(i+half_meas)], $
                             estimates, perror=sigma, bestnorm=bestnorm, dof=dof, /quiet)
            chisq = bestnorm
            chisq_red = bestnorm/dof
        ENDELSE

        ;Plots useful stuff to a window (USE WITH CARE!)
        IF KEYWORD_SET(check) THEN BEGIN
            x=i-half_meas+findgen(2*half_meas+1)
            toplo=img_slice[(i-half_meas):(i+half_meas)]
            erplo=err_slice[(i-half_meas):(i+half_meas)]
            yran=[min(toplo)-0.01*abs(mean(toplo)),max(toplo)+0.01*abs(mean(toplo))]
            xran=[min(x)-1,max(x)+1]
            plot,x,toplo,thick=3,yst=1,title='Time frame '+strtrim(j,2),yrange=yran,xrange=xran,xst=1,$
                                       xtitle='Pixels',ytitle='Intensity units',position=[.1,.15,.8,.9],/norm,psym=1
            oploterror,x,toplo,erplo,thick=3,psym=1
            oplot,x,mygauss_plus_linear(ind_fit,coeff),linestyle=2
            oplot,x,coeff[3]+coeff[4]*ind_fit
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
            AND (chisq lt cut_chisq) AND (sigma[1] GT 0.) AND (coeff[0] gt mini) THEN BEGIN
                xyouts,[0.82,0.82],[0.20,0.20],'Meets criteria',/norm,charsize=1.2
            ENDIF

            clearline=FIFTEENB()
            form="($,'pause',a,a)"
            PRINT, form=form, '         ',clearline
            PAUSE, /quiet
        ENDIF

        ;For Gaussian fit results to be used the coefficients have to
        ;be less than one pixel from maximum and with an error less than
        ;one pixel. Otherwise position of maximum is used with 0.5 pixel error.
        IF (ABS(coeff[1]) LT shift_cut) AND (chisq LT cut_chisq) $
            AND (sigma[1] GT 0.0) AND (sigma[1] LT 1.5) $
            AND (coeff[2] LT meas_size) AND (coeff[0] GT mini) THEN BEGIN

            IF KEYWORD_SET(weighted_mean) THEN BEGIN
                ;Copy value at peak
                peak_val = img_slice[i]
                sig_peak = err_slice[i]
            ENDIF ELSE BEGIN
                ;Calculate value at peak
                peak_val = mygauss_plus_linear(coeff[1],coeff)
                errpeak = mygauss_plus_linear(coeff[1]+sigma[1],coeff)
                sig_peak = peak_val - errpeak ;estimate of error
            ENDELSE

            new_i = ROUND(i+coeff[1])
            located.peaks[new_i,j,0:end_pk_ind] = ([i+coeff[1], peak_val, ABS(coeff[2]), $
                                                    coeff[3], coeff[4], i])[0:end_pk_ind]
            located.errs[new_i,j,0:end_err_ind] = ([sigma[1], sig_peak, sigma[2], $
                                                    sigma[3], sigma[4], 1.0])[0:end_err_ind]
            allpeaks[i,j] = 0 ;clear out the old pre-shifted peak
            allpeaks[new_i,j] = 3

            IF KEYWORD_SET(debug) THEN BEGIN
                ;save values for debugging
                loc_debug.shift[new_i,j] = ABS(coeff[1])
                loc_debug.sigma[new_i,j] = sigma[1]
                loc_debug.width[new_i,j] = coeff[2]
                loc_debug.chisq[new_i,j] = chisq
                loc_debug.coeff[new_i,j] = coeff[0]
            ENDIF
        ENDIF ELSE BEGIN
            ;defaulting to nearest pixel (either by choice or because the gaussian fit failed)
            located.peaks[i,j,0:nearest_end_pk_ind] = ([i, img_slice[i], 0.0])[0:nearest_end_pk_ind]
            located.errs[i,j,0:nearest_end_err_ind] = ([0.5, err_slice[i], 0.0])[0:nearest_end_err_ind]
            allpeaks[i,j] = 2

            IF KEYWORD_SET(debug) THEN BEGIN
                ;save values for debugging
                loc_debug.shift[i,j] = ABS(coeff[1])
                loc_debug.sigma[i,j] = sigma[1]
                loc_debug.width[i,j] = coeff[2]
                loc_debug.chisq[i,j] = chisq
                loc_debug.coeff[i,j] = coeff[0]
            ENDIF
        ENDELSE
    ENDIF
ENDFOR

;transferring more values to the 'output' structure - 'located'
out_grad = grad
located.allpeaks = allpeaks
located.grad_left = grad_left
located.grad_right = grad_right
END
