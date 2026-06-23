# SPDX-License-Identifier: MIT
module CubeFitter
###============###
export calculate_moments, fit_cube, fit_spectrum_from_subcube
export make_spectrum_from_cutout, make_lines_mask, toggle_fnu_flam
export write_maps_to_fits, write_spectral_cube_to_fits, quickload_neblines
export quicklook_slice, quicklook_model, quicklook_fit_result_dict
###============================================================###
using Measurements: result, uncertainty, value
# Basic computing functionality, misc.
using Base: NullLogger
using Base.Threads, Printf, Logging, LoggingExtras
# Handle NaN's more gracefully than standard Julia
using NaNStatistics
# Physical units and unit conversion, and uncertainties
using Unitful 
using UnitfulAstro 
using UnitfulEquivalences 
using Measurements
using Uncertain
using Base.Threads
# Is the next one really necessary, or should I just hardcode the constants?
using PhysicalConstants.CODATA2018
# Data formats, table, the works
using FITSIO
using DataFrames
using DataStructures 
using CSV
# Interpolation, modeling and fitting
using Polynomials 
using FastRunningMedian
using GModelFit  # TODO: Clean up GMF imports/usings!
using GModelFit: AbstractComponent, Parameter, evaluate!
import GModelFit as gmf
import Interpolations: LinearInterpolation
# Numerical integration
using NumericalIntegration
import NumericalIntegration as nui
import Plots
# # Pretty but doesn't run in notebooks and also not good with threads:
# using Term.Progress  
using ProgressMeter
ProgressMeter.ijulia_behavior(:clear)  # TODO: Check up whether this is still necessary
# Spectroscopic helper functions
include("./SpecHelpers.jl")
using .SpecHelpers
export load_neblines, load_fits_cube, load_slices_dict

# Physical constants
const ckms = SpeedOfLightInVacuum |> u"km/s"
const caps = SpeedOfLightInVacuum |> u"Å/s"
datapath = joinpath(dirname(pathof(CubeFitter)), "..", "static_data")
include("./ContSubt.jl")
export cont_subt, cont_subt_spectrum

# Instrument specific types and functions:
include("./SpectralCubes.jl")
include("./NIRSpec.jl")
include("./MUSE.jl")
export AbstractSpectralCube, NIRSpecCube, MUSECube, MIRICube


"""    quickload_neblines()
Load the default spectral line table.
"""
function quickload_neblines()
    df = load_neblines(joinpath(datapath, "neblines.dat"))
    return df
end


"""    set_logging_level(level)
Set how much feedback you want to see.
## Arguments
- `level::Symbol`: The desired level of logging to be activaged.
Valid values are `:none`, `:info`, `:warning`, or `:debug`.
"""
function set_logging_level(
    level::Symbol; out=:console::Symbol, filename="CubeFitter.log"::String)
    if out == :file
        logger = FileLogger
    else
        logger = ConsoleLogger
        filename = stderr
    end 
    if level == :info
        global_logger(ConsoleLogger(stderr, Logging.Info))
    elseif level == :none
        global_logger(NullLogger());
    elseif level == :warning
        global_logger(ConsoleLogger(stderr, Logging.Warn))
    elseif level == :debug
        global_logger(ConsoleLogger(stderr, Logging.Debug))
    else
        println("""Please chose a valid level.""")
    end 
end 

# global_logger(ConsoleLogger(stderr, Logging.Info))  # Default
global_logger(NullLogger())
# global_logger(FileLogger("CubeFitter.log"))  # Default



"""    
    estimate_redshift_reference_line_flux(cube; kwargs...)
Estimate the redshift and line flux of a line given data and line ID.
This requires a reasonably accurate initial guess to be set, but 
gives a more accurate estimate. 

# Arguments
- `cube::AbstractSpectralCube`: Datacube in question.

# Keywords
- `xrange::UnitRange`: the range of x coordinates to extract spectrum from; 
  defaults to `nothing`, which includes everything.
- `yrange::UnitRange`: Same as `xrange`, but in the y direction.
- `mask::Array`: Mask determining spaxels to extract spectrum from; takes 
  precedence over `xrange` and `yrange`. 

# Returns
- Tuple{Float64, Float64}: Tuple containing the estimated values of line flux
  (as fraction of Hα) and redshift.
"""
function estimate_redshift_reference_line_flux(spec_cube; xrange=nothing, yrange=nothing, mask=nothing)
    wave = spec_cube.wave
    if !(mask isa Nothing); xrange=nothing; yrange=nothing; end
    spec_guess, _ = make_spectrum_from_cutout(spec_cube; xrange=xrange, yrange=yrange, mask=mask)
    neblines = spec_cube.linelist
    λ_obs = neblines[neblines.name .== string(spec_cube.ref_line), :lamvac][1] * (1 +spec_cube.z_init)
    iha_g = λ_obs - 15 .< spec_cube.wave .< λ_obs + 15
    fha_g = abs(nui.integrate(wave[iha_g], spec_guess[iha_g]))
    reds_g = (nansum(wave[iha_g] .* spec_guess[iha_g]) ./ nansum(spec_guess[iha_g])
        ./  neblines[neblines.name .== string(spec_cube.ref_line), :lamvac][1] - 1)
    @debug "Flux guess" fha_g reds_g
    return fha_g, reds_g
end 


"""    get_resolving_power(instrument[, setting=nothing])
This function figures out which instrument we are working with, and calls the relevant
instrument-specific function.
"""
function get_resolving_power(instrument::String; setting=nothing::Any)
    if lowercase(instrument) == "nirspec"
        return(get_resolving_power_nirspec(setting))
    elseif lowercase(instrument) == "muse"
        return get_resolving_power_muse(order=setting)
    else
    end
end


# """     get_resolving_power_nirspec(grating[, calib_path="../static_data/"])
# This function assumes that the FITS files containing the NIRSpec grating resolving power
# are named "jwst_nirspec_[grating]_disp.fits", as fetched from the official JDox.
# # Arguments
# - `grating::String`: The name of the grating used.
# # Keywords
# - `calib_path::String`: Path to the location of the calibration data, defaults to  the data folder installed with the package.
# # Returns
# - A callable `Interpolations.Extrapolation` objects which returns the linearly 
#   interpolated resolving power at the wavelength passed to it.
# """
# function get_resolving_power_nirspec(grating::String; calib_path=datapath) 
#     lsffilename = "jwst_nirspec_$(grating)_disp.fits"
#     @debug "Path to LSF file: " calib_path * lsffilename
#     fitsfile = FITS(joinpath(calib_path, lsffilename))
#     angwave = read(fitsfile[2], "WAVELENGTH") * u"μm" .|> u"angstrom"
#     R = read(fitsfile[2], "R")
#     @debug extrema(angwave) extrema(R)
#     itp = LinearInterpolation(ustrip.(angwave), R, extrapolation_bc=Flat())
#     return itp
# end


# """     get_resolving_power_muse(; order=3)
# ## Optional arguments
# - `order::Int`: The order of the polynomial to fit to the datapoints. Default is 3.
# ## Returns
# - A callable `Polynomials.Polynomial` object which returns the linearly interpolated
#   resolving power at the wavelength passed to it.
# """
# function get_resolving_power_muse(;order=3::Int)
#     ll = Array{Float64}(
#         [4650.0, 5000.0, 5500.0, 6000.0, 6500.0, 7000.0,
#          7500.0, 8000.0, 8500.0, 9000.0, 9350.0])
#     rr = Array{Float64}(
#         [1609, 1750, 1978, 2227, 2484, 2737,
#          2975, 3183, 3350, 3465, 3506])
#     drr = Array{Float64}([6, 4, 6, 6, 5, 4, 4, 4, 4, 5, 10])
#     p = Polynomials.fit(ll, rr, order, weights=1. ./ drr)
#     return p
# end
#


"""    convert_ergscms_Å_units(datadict::Dict; instrument="NIRSpec")
Makes a best possible attempt at finding out what the data units are of the cube,
and to convert those data units into Ångström and ergs/s/cm²/Å.
"""
function convert_ergscms_Å_units(datadict::Dict; instrument="NIRSpec")
    data, errs, wave, header = (
        datadict[:Data], datadict[:Errs], datadict[:Wave], datadict[:Header])
    wave = Array{Float64}(wave)
    if lowercase(instrument) in ["muse"]
        data .*= header["PHOTFLAM"]
        data .*= u"erg/s/cm^2/Å"
        errs .*= header["PHOTFLAM"]
        errs .*= u"erg/s/cm^2/Å"
        wave .*= u"Å"
    elseif lowercase(instrument) in ["nirspec"]
        data .*= header["PIXAR_SR"]
        data *= u"MJy"
        errs .*= header["PIXAR_SR"]
        errs *= u"MJy"
        wave = wave .* u"μm" .|> u"Å"
        @debug "Data dimension: " unique(dimension.(data))
        @debug "Data unit: " unique(unit.(data))
        @debug "wave unit: " unit(wave[1])
    end
    # Make sure data and waves have correct units
    flam_unit = u"erg/s/cm^2/Å"
    if unique(dimension.(data)) != dimension(flam_unit)
        #@info "Data unit dimension not equivalent to $flam_unit, converting."
        data = toggle_fnu_flam(data, wave)
        errs = toggle_fnu_flam(errs, wave)
    else
        #@info "Data unit was equal to or equivalent to $flam_unit"  #("Smooth dims bro!")
    end

    data = uconvert.(flam_unit, data)
    errs = uconvert.(flam_unit, errs)

    # Update the datacict to have better formatted output
    datadict[:Data] = data
    datadict[:Errs] = errs
    datadict[:Wave] = wave
    return datadict  #, data, errs, wave, header
end


"""    toggle_fnu_flam(flux, wave, verbose=false::Bool)
Converts fnu data to flam units, and vice versa.
Input must be Unitful™ Quantities.
"""
function toggle_fnu_flam(influx, wave; verbose=false::Bool)::Array
    # Broadcast spectrum or cube. A bit clumsy but works for now.
    if ndims(influx) > ndims(wave)
        wave = reshape(Array(wave), 1, 1, length(wave))
    end 
    # Establish base units and dimensions
    fnu_base_unit = u"erg/s/cm^2/Hz"
    flam_base_unit = u"erg/s/cm^2/Å"
    fnudim = dimension(fnu_base_unit)
    flamdim = dimension(flam_base_unit)
    # Convert wave to Å units.
    # `if` -clause because `Spectral()` cannot handle
    # two quantities of same dimension, for some reason.
    # println(dimension(wave[1]))
    if dimension(wave[1]) == dimension(u"Å")
        wave = uconvert.(u"Å", wave)
    else 
        wave = uconvert.(u"Å", wave, Spectral())
    end 
    conversion = ustrip.(wave).^2 ./ ustrip(caps)
    if verbose == true
        println(nanminimum(conversion), ",  ", nanmaximum(conversion))
        println(dimension(influx[1]), fnudim, flamdim)
        println(dimension(influx[1]) == fnudim)
        println(dimension(influx[1]) == flamdim)
    end 
    if dimension(influx[1]) == fnudim
        fnu = uconvert.(fnu_base_unit, influx)
        outflux = ustrip.(fnu) ./ conversion .* flam_base_unit
    elseif dimension(influx[1]) == flamdim
        flux = uconvert.(flam_base_unit, influx)
        outflux = ustrip.(flux) .* conversion .* fnu_base_unit
    end 
    return outflux
end


"""    
    make_spectrum_from_cutout(cube; kwargs...)

Extracts a spectrum either from rectangular pixel ranges or from a segment
defined by a boolean/bitwise mask. If the mask is instead an array, the
normalized number will be used as weights. This function is user-facing, but also used 
internally by `CubeFitter`; passing `format_nice=true` will give the nicer output for 
the user.

# Arguments
- `cube::AbstractSpectralCube`: The datacube to work on.

# Keywords
- `xrange`, `yrange`: Pixel coordinate ranges defining a rectangle from which to 
  extract the spectrum; both default to `nothing`.
- `mask::BitMatrix=nothing`: Alternative to `xrange` and `yrange`, a boolean mask 
  determining spaxels to include when extracting spectrum. Must match the spatial 
  dimensions of the cube. If given, yrange and xrange are ignored.
- `weights::Array{Float}=nothing`: Weights to give the spaxels when extracting spectrum. 
  If not given, the spaxels in the mask/range will all be weighed equally.
- `format_nice::Bool=false`

# Returns
- `Tuple{Vector{Float64}, Vector{Float64}}`: If `format_nice` is set to `false`, returns s  Tuple of two vectors containing the fluxes and errors of the extracted spectrum. 
- `DataFrames.DataFrame`: Pandas- or R-like DataFrame containing wave, flux and errors of
  the spectrum if `format_nice` is set to `true`
"""
function make_spectrum_from_cutout(cube; xrange=nothing, yrange=nothing, mask=nothing, weights=nothing, format_nice=false)
    if xrange isa Nothing; xrange = (1:size(cube.fluxcube)[1]); end
    if yrange isa Nothing; yrange = (1:size(cube.fluxcube)[2]); end
    if !(mask isa BitMatrix)
      if !(weights isa Nothing)
        @warn "Weights are given explicitly but also implicitly as a mask with numerical values.
        Using the explicitly given weights, using the mask just as a mask."
      end
      if !(mask isa Nothing)
        weights = mask |> copy
        mask = ((mask .!= 0) .& .!(mask .|> isnan))
      end
    end
    if !(weights isa Nothing)
      cube.fluxcube .*= weights
      cube.errscube .*= weights
    end
    if !(mask isa Nothing)
      cutout_spec = cube.fluxcube[mask, :]
      cutout_errs = cube.errscube[mask, :]
      slice = cutout_spec[:, 1]
      num_spaxels = length(slice)
      spec_init = nanmean(cutout_spec, dim=(1))
      errs_init = sqrt.(nansum(cutout_errs.^2, dim=(1))) / num_spaxels
    else
      cutout_spec = cube.fluxcube[xrange, yrange, :]
      cutout_errs = cube.errscube[xrange, yrange, :]
      slice = cutout_spec[:,:,1]
      num_spaxels = length(slice)
      if num_spaxels == 1  #1
        spec_init = vec(cutout_spec)
        errs_init = vec(cutout_errs)
      else
        spec_init = nanmean(cutout_spec, dim=(1,2))
        errs_init = sqrt.(nansum(cutout_errs.^2, dim=(1,2))) / num_spaxels
      end
    end
    return if format_nice
      OrderedDict(
        :wave => cube.wave, 
        :flam => spec_init,
        :dflam => errs_init,
      ) |> DataFrame
    else
      spec_init, errs_init
    end
    # if format_nice:
    # end
    # return spec_init, errs_init
end


"""    fit_spectrum_from_subcube(cube; xrange=nothing, yrange=nothing)

Extract a spatially rectangular subcube defined from a range in x- and y coordinates, 
coadd it as a spectrum, and pass it to the 1D spectrum fitter. 

This is the main workhorse of the user-facing functions in this package. Given any two
pixel ranges, it extrancts the spectrum and error vector from the cube in these ranges. If
not passed any range in a direction, it uses the entire available pixel range in that
direction.

To only use one row, column or spaxel, just pass in its pixel value as both start and end
of the range.

For a convenient wrapper around this function which fits an entire cube spaxel-by-spaxel,
see `fit_cube`.

# Arguments
Required arguments:
- `cube::AbstractSpectralCube`: The spectral cube object to fit. Must be one of the
  instrument-specific types of struct defined by this package.
Keyword arguments:
- `xrange::UnitRange{Int}`: A range (e.g. `2:3`) of x-axis pixel coordinates to include in
  the extracted spectrum.
- `yrange::UnitRange{Int}`: Same as `xrange`, but for y-axis coordinates.
- `broad_component::Bool`: Whether to include a second, broad, component.
- `line_selection::AbstractArray{Union{String,Symbol}}`: Lines to fit.
- `kinematics_from::AbstractArray{Union{String,Symbol}}`: Lines on which to base kinematic fit.
- `min_snr::Number`: SNR limit for for numeric flux estimate to discard chunk.

# Returns
- `Dict`: A dict containing wave, spectrum, errors, fit result and fit statistics.

# Examples
Extract and fit spectrum for the 10th to 15th columns, 12th to 17th row of a cube:
```julia-repl
julia> fit_spectrum_from_subcube(cube, xrange=10:15, yrange=12:17)
```
Extract the spectrum from the spaxel with coordinates 33, 44:
```julia-repl
julia> fit_spectrum_from_subcube(cube, xrange=33:33, yrange=44:44)
```
"""
function fit_spectrum_from_subcube(cube; xrange=nothing, yrange=nothing, broad_component=false,
    line_selection=nothing, kinematics_from=nothing, min_snr=0.5)::Dict
    if isa(xrange, Nothing); xrange=range(1, size(cube.fluxcube)[1]); end
    if isa(yrange, Nothing); yrange=range(1, size(cube.fluxcube)[2]); end
    wave_init = cube.wave
    spec, errs = make_spectrum_from_cutout(cube; xrange=xrange, yrange=yrange)
    @debug "Before NaNmask10'ing: " count(isnan.(spec)) count(isnan.(errs))
    notnan = nanmask10(spec, errs)
    goodwave, goodspec, gooderrs = wave_init[notnan], spec[notnan], errs[notnan]
    # @debug xrange yrange goodwave, count(notnan)
    # @debug "Goodwave: " goodwave
    mod = build_model(cube, xrange=xrange, yrange=yrange, dom_init=Domain(goodwave),
        broad_component=broad_component, line_selection=line_selection, kinematics_from=kinematics_from, min_snr=min_snr)
    @debug "First pass at making model done!" mod
    idx = make_lines_mask(mod, Domain(goodwave))
    # @debug "Made lines mask for $xrange, $yrange !"
    lwave, lspec, lerrs = goodwave[idx], goodspec[idx], gooderrs[idx]
    # @debug "Before and after line masking: " length(goodwave), length(lwave)
    if length(lwave) == 0; println("`fit_spectrum_from_subcube`: No wavelengths not masked") ;return NaN; end
    lmod = build_model(cube, xrange=xrange, yrange=yrange, dom_init=Domain(lwave),
        broad_component=broad_component, line_selection=line_selection, kinematics_from=kinematics_from, min_snr=min_snr)
    lmeas = Measures(Domain(lwave), lspec, lerrs)
    # return lmod
    # return lmod, lmeas
    @debug "Model before fitting: $xrange, $yrange:  " lmod lmeas lwave
    if !(cube.ref_line in keys(lmod)); return NaN; end
    # DEBUG:
    @debug lmod, lmeas
    @debug "Now trying for a fit" xrange yrange
    @debug lmod
    # results, stats = gmf.fit(lmod, lmeas, gmf.lsqfit())
    # results, stats = gmf.fit(lmod, lmeas, minimizer=gmf.lsqfit())
    #@info "Successfully fitted (col. $xrange, row $yrange)"
    # @debug "Results and stats!" results stats
    # END: 
    try  # Too many things can go wrong to catch eatch one separately.
        results, stats = gmf.fit(lmod, lmeas, gmf.lsqfit())
        # results, stats = gmf.fit(lmod, lmeas, minimizer=gmf.lsqfit())
        #@info "Successfully fitted (col. $xrange, row $yrange)"
        @debug "Results and stats!" results stats
        output = Dict(
            :wave => goodwave,
            :spec => goodspec,
            :errs => gooderrs,
            :measure => lmeas,
            :results => results,
            :fitstats => stats)
        return output
    catch
        return lmod
    end 
end


"""    fit_spectrum_from_mask(cube mask; kwargs...)

Extract a subcube defined by a spaxel mask (optionally with weights) as a spectrum, 
and pass it to the 1D spectrum fitter. 

This is the other version of the main workhorse of the user-facing functions in 
this package. Instead of spaxel ranges, it uses a user-provided mask to define , from 
which spaxels the spectrum is extracted; but otherwise works exactly like 
`fit_spectrum_from_subcube`.

For a convenient wrapper around this function which fits an entire cube spaxel-by-spaxel,
see `fit_cube`.

# Arguments
Required arguments:
- `cube::AbstractSpectralCube`: The spectral cube object to fit. Must be one of the
  instrument-specific types of struct defined by this package.
- `mask::BitMatrix`: Mask to use to extract spectrum.
Keyword arguments:
- `broad_component::Bool`: Whether to include a second, broad, component.
- `line_selection::AbstractArray{Union{String,Symbol}}`: Lines to fit.
- `kinematics_from::AbstractArray{Union{String,Symbol}}`: Lines on which to base kinematic fit.
- `min_snr::Number`: SNR limit for for numeric flux estimate to discard chunk.

# Returns
- `Dict`: A dict containing wave, spectrum, errors, fit result and fit statistics.

# Examples
Extract and fit spectrum from spaxels defined by the mask `mask`:
```julia-repl
julia> fitdict = fit_spectrum_from_mask(cube, mask)
```
"""
function fit_spectrum_from_mask(cube, mask; broad_component=false, line_selection=nothing, 
    kinematics_from=nothing, min_snr=0.5)
    # if isa(xrange, Nothing); xrange=range(1, size(cube.fluxcube)[1]); end
    # if isa(yrange, Nothing); yrange=range(1, size(cube.fluxcube)[2]); end
    wave_init = cube.wave
    spec, errs = make_spectrum_from_cutout(cube; mask=mask)
    # @debug "Before NaNmask10'ing: " count(isnan.(spec)) count(isnan.(errs))
    @debug "Before NaNmask10'ing: " count(!isnan.(spec)) count(!isnan.(errs)) spec errs
    notnan = nanmask10(spec, errs)
    goodwave, goodspec, gooderrs = wave_init[notnan], spec[notnan], errs[notnan]
    # @debug xrange yrange goodwave, count(notnan)
    # @debug "Goodwave: " goodwave
    mod = build_model(cube, mask=mask, dom_init=Domain(goodwave),
        broad_component=broad_component, line_selection=line_selection, kinematics_from=kinematics_from, min_snr=min_snr)
    @debug "First pass at making model done!" mod
    idx = make_lines_mask(mod, Domain(goodwave))
    # @debug "Made lines mask for $xrange, $yrange !"
    lwave, lspec, lerrs = goodwave[idx], goodspec[idx], gooderrs[idx]
    # @debug "Before and after line masking: " length(goodwave), length(lwave)
    if length(lwave) == 0; return NaN; end
    lmod = build_model(cube, mask=mask, dom_init=Domain(lwave),
        broad_component=broad_component, line_selection=line_selection, kinematics_from=kinematics_from, min_snr=min_snr)
    lmeas = Measures(Domain(lwave), lspec, lerrs)
    # return lmod
    # return lmod, lmeas
    @debug "Model before fitting: $xrange, $yrange:  " lmod lmeas lwave
    if !(cube.ref_line in keys(lmod)); return NaN; end
    # DEBUG:
    @debug lmod, lmeas
    @debug "Now trying for a fit" xrange yrange
    @debug lmod
    # results, stats = gmf.fit(lmod, lmeas, gmf.lsqfit())
    # results, stats = gmf.fit(lmod, lmeas, minimizer=gmf.lsqfit())
    #@info "Successfully fitted (col. $xrange, row $yrange)"
    # @debug "Results and stats!" results stats
    # END: 
    try  # Too many things can go wrong to catch eatch one separately.
        results, stats = gmf.fit(lmod, lmeas, gmf.lsqfit())
        # results, stats = gmf.fit(lmod, lmeas, minimizer=gmf.lsqfit())
        #@info "Successfully fitted (col. $xrange, row $yrange)"
        @debug "Results and stats!" results stats
        output = Dict(
            :wave => goodwave,
            :spec => goodspec,
            :errs => gooderrs,
            :measure => lmeas,
            :results => results,
            :fitstats => stats)
        return output
    catch
        return lmod
    end 
end


function makeresultdict(cube; broad_component=false, lines_select=nothing)
    @debug "Making output dict"
    xsize, ysize = size(cube.fluxcube)[1], size(cube.fluxcube)[2]
    neblines = cube.linelist
    wave = cube.wave
    if lines_select isa Nothing; lines_select=Symbol.(neblines[:, :name]); end
    # Create dictionary to hold the fit results.
    slices_dict = Dict()
    slices_dict[:refline] = cube.ref_line
    slices_dict[:header] = cube.header
    for scomp in neblines.name
        if ! (Symbol(scomp) in lines_select); continue; end
        lamvac = neblines[neblines.name .== scomp, :lamvac][1]
        obsvac = lamvac * (1 + cube.z_init)
        if (minimum(ustrip.(wave)) > obsvac) | (maximum(ustrip.(wave)) < obsvac)
            # @debug "Line $scomp outside wavelength range, skipping"
            continue
        end 
        comp = Symbol(scomp)
        if comp == :main; continue; end
        slice = zeros(xsize, ysize, 5) .* NaN
        slices_dict[comp] = slice
        if broad_component
           slices_dict[Symbol(String(comp) * "_broad")] = copy(slice)
        end
    end 
    slices_dict[:redshift] = zeros(xsize, ysize, 3) .* NaN
    slices_dict[:fwhm] = zeros(xsize, ysize, 3) .* NaN
    slices_dict[:mom0] = zeros(xsize, ysize, 3) .* NaN
    slices_dict[:mom1] = zeros(xsize, ysize, 3) .* NaN
    slices_dict[:mom2] = zeros(xsize, ysize, 3) .* NaN
    slices_dict[:fitstats] = zeros(xsize, ysize) .* NaN
    if broad_component
        slices_dict[:refline_broad] = Symbol(String(cube.ref_line) * "_broad")
        slices_dict[:redshift_broad] = zeros(xsize, ysize, 3) .* NaN
        slices_dict[:fwhm_broad] = zeros(xsize, ysize, 3) .* NaN
    end
    @debug "Output dict successfully done!"
    return slices_dict
end


"""    
    fit_cube(cube; kwargs...)
Fit the entire cube spaxel-by-spaxel or segment-by-segment.

This is a convenience function to perform that one action we want to perform 90% of the time.
Still, there are a few options to tweak:

# Arguments
- `cube::AbstractSpectralCube`: The SpectralCube object to fit.

# Keywords 
- `broad_component::Bool`: Wether or not to include a second kinematic component in the fit.
- `min_snr::Float`: Where the numerically estimated S/N ratio per spaxel is below this 
  value, the line will not be included in the fit for this line, and the value set to NaN.
- `line_selection::List{Symbol}`: A list of lines to fit. Defaults to all lines in the 
  loaded line list which fall within the wavelength range of the dataset.
- `kinematics_from_lines::List{Symbol}`: List of names of lines on which to base the 
  kinematc fits. Other lines will be fit with locked kinematics, assuming they share the 
  same properties as the ones in this list, allowing only the flux to vary.
  *NB*: Blended lines (e.g. [O II]λλ3727,29 or Hα + [N II]λλ6548,84) should be included 
  in this selection to ensure they are fitted simultaneously.
- `segmap::Array{Int}`: Optional segmentation map. If passed, an average will be made per 
  segment instead of a measurement for each spaxel.
- `window_kms::Float`: The width (in km/s) of the window to integrate over to find numerical 
  fluxes. Can be important for neighboring lines.

# Returns
- Dict: Dictionary containing maps of fluxes, kinematic properties, and fit statistics.
"""
function fit_cube(cube; broad_component=false, min_snr=1.0, line_selection=nothing, kinematics_from_lines=nothing, segmap=nothing, window_kms=1000.0)
    # xsize, ysize = size(cube.fluxcube)[1], size(cube.fluxcube)[2]
    slices_dict = makeresultdict(cube, broad_component=broad_component)
    user_segmap = segmap

    if kinematics_from_lines isa Nothing
        model_lines = line_selection
    else
        model_lines = kinematics_from_lines
    end

    if user_segmap isa Nothing
        # Make every spaxel a "segment"
        segmap = reshape(collect(1:length(cube.fluxcube[:,:,1])), size(cube.fluxcube)[1:end-1])
    else 
        # Use user provided segmentation map
        segmap = user_segmap
    end

    @debug "Made slices_dict"

    @showprogress @threads for f in (segmap |> unique |> sort)
        thismask = segmap .== f
        # print("Fitting segment No. $f with $(thismask |> sum) spaxels \r")
        numnan = count(isnan.(cube.fluxcube[thismask, :]))
        testspec, testerr = make_spectrum_from_cutout(cube; mask=thismask)
        nanratio = count(isnan.(testspec)) / length(testspec)
        if nanratio > 0.5
            # println("ratio of nans  $nanratio  in segment $f                 ")
            continue
        end
        # println("There were $badcount empty segments out of $(segmap |> unique |> count) in total")
        fitdict = fit_spectrum_from_mask(
            cube, thismask; broad_component=broad_component, line_selection=model_lines)
        @debug "What's a fitdict anyway?" fitdict typeof(fitdict)
        if !isa(fitdict, Dict)
            @debug "Fitdict $f was not a dict!"
            continue
        end
        linesdict = Dict(l => fitdict[:results][l] for l in fitdict[:results] |> keys if l != :main)
        if !(kinematics_from_lines isa Nothing)
            minwave, maxwave = (cube.wave |> minimum, cube.wave |> maximum) ./ (1 + cube.z_init)
            if line_selection isa Nothing
                line_idx = minwave .< cube.linelist[:, :lamvac] .< maxwave
                line_selection = cube.linelist[line_idx, :name] .|> Symbol
            end
            forced_lines = setdiff(line_selection, kinematics_from_lines)
            for l in forced_lines
                outdict = fit_spectrum_from_mask(
                    cube, thismask; broad_component=broad_component, line_selection=[l],
                    kinematics_from=fitdict[:results], min_snr=min_snr)
                # Do not fail gracefully, just flip the table and move on.
                # if !(outdict isa Dict); println("Failed segment $f with $(count(thismask)) spaxels"); continue; end
                if !(outdict isa Dict); continue; end
                if !(l in outdict[:results] |> keys); continue; end
                oo = outdict[:results][l]
                linesdict[l] = oo
            end
        end
        fitdict[:linesdict] = linesdict
        fill_in_fit_values!(slices_dict, fitdict, mask=thismask)
    end 
    mom0, mom1, mom2 = calculate_moments(cube; window_kms=window_kms |> float)
    if !(user_segmap isa Nothing)
        mom0 = mean_by_segmap(mom0 .|> value, user_segmap, errs=mom0 .|> uncertainty)
        mom1 = mean_by_segmap(mom1 .|> value, user_segmap, errs=mom1 .|> uncertainty)
        mom2 = mean_by_segmap(mom2 .|> value, user_segmap, errs=mom2 .|> uncertainty)
    end
    # Moments come with both signal and noise now
    slices_dict[:mom0] = stack([value.(mom0), uncertainty.(mom0)], dims=3)
    slices_dict[:mom1] = stack([value.(mom1), uncertainty.(mom1)], dims=3)
    slices_dict[:mom2] = stack([value.(mom2), uncertainty.(mom2)], dims=3)
    # Remove empty slices
    # println(slices_dict |> keys)
    # println("Entries in the slices dict before pruning: $(slices_dict |> keys |> length)")
    @info slices_dict |> keys
    for s in slices_dict |> keys
        if !(slices_dict[s] isa Array); continue; end
        slice = slices_dict[s][:,:,1]
        # println("for slice $s, there were $(count(isnan.(slice))) nans out of $(length(slice))")
        if count(isnan.(slice)) == length(slice)
            delete!(slices_dict, s)
            continue
        end
        if !(String(s) in cube.linelist[!, "name"]); continue; end
        m0, m1, m2 = calculate_moments(cube, refline=s, window_kms=window_kms)
        if !(user_segmap isa Nothing)
            m0 = mean_by_segmap(m0 .|> value, user_segmap, errs=m0 .|> uncertainty)
        end
        m0v, m0e = m0 .|> value, m0 .|> uncertainty
        m0s = stack([m0v,m0e], dims=3)
        slices_dict[s][:,:,4:5] = m0s
    end
    return slices_dict
end 


"""    
    build_model(cube; kwargs...)

Builds the `GModelFit.Model` object necessary for fitting the lines.
The model assumes that each line is described as a single Gaussian profile, with a shared
`redshift` and velocity width `fwhm`, while their individual fluxes are left as free parameters
(including lines which actually have fixed line ratios, which might not have that in the ). 

# Arguments
- `cube::AbstractSpectralCube`: A Domain object containing the wavelength range on which
  the model should be defined.

# Keywords
- `xrange::UnitRange{Int}`: Range of x-axis pixels to include in the fit; defaults to all.
- `yrange::UnitRange{Int}`: Range of y-axis pixels to include in the fit; defaults to all.
- `mask::Array`: Mask to extract spectrum from, alternative to `xrange` and `yrange`. 
  The `mask` takes precedence over `xrange` and `yrange`.
- `min_snr::Float64`: Minimum estimated SNR for a line to be included in the fit. Defaults 
  to 1.5. This is the SNR that is estimated numerically; set it to 0 to make sure no lines 
  are dropped when building the model.
- `fwhm_int::Float64`: The initial guess of FWHM of the line in km/s. Defaults to 100.
- `dom_init::GModelFit.Domain`: For finer control of which wavelength ranges can be included
  in the fit, pass a domain (will default to the full wl range in the cube otherwise). This 
  is the maximum allowed domain to include, parts of it might still be excluded by the function.
- `broad_component::Bool`: Whether to include a second, broad component in the model.
  So far, only up to 2 components is supported.
- `line_selection::Vector{Symbol}`: List of lines to fit; defaults to all the lines in the 
  loaded line list which fall inside the wavelength range of the data.

# Returns
- `model::GModelFit.Model`: The model containing the prescribed lines. The Reference line has
  norm=1., and the others are scaled according to the values of `foverha` given in the line list
  input file.
"""
function build_model(cube; xrange=nothing, yrange=nothing, mask=nothing, min_snr=0.5, fwhm_int=100,
    dom_init=nothing, broad_component=false, line_selection=nothing, kinematics_from=nothing, lock_kinematics=false)
    redss = (1. + cube.z_init)  # Just for convenience
    wave = cube.wave
    neblines = cube.linelist
    if line_selection isa Nothing; line_selection = Symbol.(neblines[:, "name"]); end
    if !(kinematics_from isa Nothing); lock_kinematics=true; end
    if broad_component & !(kinematics_from isa Nothing)
        println("Building kinematics from existing model is not implemented for two components yet")
        return NaN
    end
    union!(line_selection, [cube.ref_line])
    if isa(dom_init, Nothing); dom_init = Domain(wave); end
    @debug "build_model() domain: " xrange yrange dom_init
    spec_init, errs_init = make_spectrum_from_cutout(cube; xrange=xrange, yrange=yrange, mask=mask)
    components = OrderedDict()
    complist = Vector{Symbol}()
    for l in neblines.name[1:end]
        if !(Symbol(l) in line_selection); continue; end
        lamvac  = neblines[neblines.name.==l, :lamvac][1]
        lam_obs = lamvac * redss
        @debug "Build model: domain min max" (dom_init |> coords |> extrema)
        if !(minimum(coords(dom_init)) + 10 < lam_obs < maximum(coords(dom_init)) - 10)
            @debug "`build_model()` for $xrange,$yrange: Line $l was outside the wavelength domain"
            continue
        end 
        fitwindow = v_to_deltawl(1000, lam_obs) # 1000
        idx = lam_obs - fitwindow .< wave .< lam_obs + fitwindow
        if length(wave[idx]) < 10; continue; end
        snr = estimate_line_snr(wave[idx], spec_init[idx], err=errs_init[idx])
        line_too_faint = (snr < min_snr) | isnan(snr)
        if line_too_faint
            if l != cube.ref_line
                @debug   "$l was too faint with SNR = $snr"
                continue
            end 
        end
        fha_g, reds_g = estimate_redshift_reference_line_flux(
            cube, xrange=xrange, yrange=yrange, mask=mask)
        if isnan(fha_g) || isnan(reds_g)
            @debug "Numerical integration of line $l failed due to lack of data." l
            continue
        end
        norm = estimate_line_snr(wave[idx], spec_init[idx], err=errs_init[idx], only_flux=true)
        fwhm_inst = cube.lsf_fitter(lam_obs)
        # @debug "Resolving power at line $l: " fwhm_inst
        if kinematics_from isa Nothing
            fwhm_tot = sqrt(fwhm_int^2)  # + fwhm_int^2)
            zz = cube.z_init
        else
            fwhm_tot = sqrt(kinematics_from[cube.ref_line].fwhm_kms.val^2)  # ^2 + fwhm_inst^2)
            zz = kinematics_from[cube.ref_line].redshift.val
        end
        g = GModelFit.FComp(_gauss_line, [:waves];
            lab_wave=lamvac, redshift=zz,
            fwhm_kms=fwhm_tot, norm=norm, lsf=fwhm_inst)
        components[Symbol(l)] = g # = ref_line
        append!(complist, [Symbol(l)])
        
        if broad_component
            g2 = GModelFit.FComp(
                _gauss_line, [:waves], lab_wave=lamvac, redshift=zz,
                fwhm_kms=fwhm_tot*1.5, norm=norm*0.3, lsf=fwhm_inst)
            components[Symbol(l * "_broad")] = g2 # = ref_line
            append!(complist, [Symbol(l * "_broad")])
        end
        @debug "Added the line $l to the model"
    end 
    # First, define model with mutually independent components
    # @debug "What if broad model?" dom_init components
    components[:main] = SumReducer(complist)
    @debug "Components before building model" components
    mod = Model(components)
    # mod[:main] = SumReducer(complist)
    @debug "SumReducer succeeded!" components
    # Now: We set initial values, locked variables and interdependencies.
    for comp in keys(mod)
        @debug "Setting up $comp"
        # INFO: If a component is patched to itself, it in essence gets fixed!!
        # So definitely do not let the reference line get `patch`'ed! 
        # TODO: Find a more elegant solution for the `if` statement below.
        if comp == :main
            @debug "Comp was `:main`, moving on"
            continue
        elseif !(comp in [Symbol(cube.ref_line), Symbol(String(cube.ref_line) * "_broad")])
            @debug "Fixing parameters for $comp"
            if endswith(String(comp), "_broad")
                refline = Symbol(String(cube.ref_line) * "_broad")
            else
                refline = Symbol(String(cube.ref_line))
            end 
            mod[comp].lab_wave.fixed = true
            mod[comp].redshift.low = cube.z_init - cube.z_init * 0.1
            mod[comp].redshift.high = cube.z_init + cube.z_init * 0.1

            if kinematics_from isa Nothing
                mod[comp].redshift.patch = refline
                mod[comp].fwhm_kms.patch = refline
            else
                mod[comp].fwhm_kms.fixed = true
                mod[comp].redshift.fixed = true
            end
            mod[comp].norm.low = 0.
            mod[comp].lsf.fixed = true
        else
            @debug "Fixing parameters for $comp"
            mod[comp].fwhm_kms.low = 1.
            mod[comp].fwhm_kms.high = 2000.
            mod[comp].redshift.low = cube.z_init - cube.z_init * 0.1
            mod[comp].redshift.high = cube.z_init + cube.z_init * 0.1
            mod[comp].lab_wave.fixed = true
            mod[comp].norm.low = 0.
            mod[comp].lsf.fixed = true
            if !(kinematics_from isa Nothing)
                mod[comp].redshift.fixed = true
                mod[comp].fwhm_kms.fixed = true
            end
            if endswith(String(comp), "_broad")
                mod[comp].redshift.val = cube.z_init - cube.z_init * 0.001
                mod[comp].fwhm_kms.patch = @fd (m, v) -> v + m[Symbol(cube.ref_line)].fwhm_kms
                mod[comp].fwhm_kms.low = 0.
            end 
        end 
    end 
    @debug "`build_model()` successfully done!" mod
    return mod
end 


"""    _gauss_line(waves, lab_wave, redshift, fwhm_kms, norm, lsf)
"""
function _gauss_line(
  waves::Array, lab_wave::Float64, redshift::Float64, fwhm_kms::Float64, norm::Float64, lsf::Float64)
  cen_wave = lab_wave * (1. + redshift)
  fwhm_inst = cen_wave / lsf
  fwhm_aa = v_to_deltawl(fwhm_kms, cen_wave)
  fwhm_tot = sqrt(fwhm_inst^2 + fwhm_aa^2)
  sigma_aa = fwhm_to_sigma(fwhm_tot)
  flux = gauss(waves, ustrip(sigma_aa), ustrip(cen_wave)) .* norm
  return flux
end

 
"""
    calculate_moments(cube; refline=nothing, window_kms=1000)
Computes moments of a flux slab with associated wavelengths.

Computes the 0th, 1st and 2nd moments of a spectral chunk (a line map) in one or more
pixels. Errors/uncertainties are not included in this computation.
# Arguments
- `cube::AbstractSpectralCube`: Datacube struct to be computed.

# Keywords
- `refline::Symbol`: The line to integrate, if not the default one set in the `cube`.
- `window_kms::Float64`: Width of the spectroscopic range included in the calculation,
  given in km/s.

# Returns
- A tuple of three 2D arrays (images) of the moments mapped to each spatial pixel.
"""
function calculate_moments(cube; refline=nothing, window_kms::Float64=1000.)
  if isa(refline, Nothing)
    refline = cube.ref_line
  end 
  # println("Reference line: $refline")
  # println("")
  lam_rest = cube.linelist[cube.linelist.name .|> String .== refline |> String, :lamvac][1]
  lam_obs = lam_rest * (1 + cube.z_init)  
  # cube.linelist[cube.linelist.name .|> String .== refline |> String, :lamvac][1] * (1 + cube.z_init)
  # println("$lam_rest, $lam_obs")
  window_ang = v_to_deltawl(window_kms/2, lam_obs)
  idx = lam_obs - window_ang .< cube.wave .< lam_obs + window_ang
  waves, flux = cube.wave[idx], cube.fluxcube[:,:,idx] .± cube.errscube[:,:,idx]
  # println("Calculate moments: Size of flux chunk is: ", size(flux))
  # Make sure things are well formatted now before they get passed along.
  @assert size(waves)[1] == size(flux)[3]
  # Empty arrays to fill in:
  @debug waves
  mom0 = similar(flux[:, :, 1])
  mom1 = similar(flux[:, :, 1])
  mom2 = similar(flux[:, :, 1])
  ncols, nrows = size(flux)
  for i in 1:ncols
    for j in 1:nrows
      theflux = flux[i, j, :]  # Just for convenience
      deltalam = waves[2:end] .- waves[1:end-1] 
      append!(deltalam, deltalam[end])
      mom0b = (deltalam .* theflux) |> nansum
      mom0[i, j] = abs(mom0b)
      mom1[i, j] = nansum(waves .* theflux) ./ nansum(theflux)
      @debug i, j theflux mom0[i, j] mom1[i, j]
      mom2[i, j] = sqrt(
        abs.(nansum((waves .- mom1[i, j]).^2 .* theflux)) /
        abs.(nansum(theflux))
      )
    end
  end 
  return mom0, mom1, mom2
end 


# Uh, these two look like they need to be cleaned out, but don't have time now. 
# My suspicion is they aren't actually used anywhere.
function moments(wave::AbstractFloat, flux::AbstractFloat, err::AbstractFloat)
    flux = flux .± err
    mom0, mom1, mom2 = moments(wave, flux)
    return mom0, mom2, mom3
end

function moments(wave::AbstractFloat, flux::Measurement{AbstractFloat})
    flux = abs.(flux)
    deltawave = wave[2:end] - wave[1:end-1] 
    append!(deltawave, deltawave[end])
    mom0 = (deltawave .* flux) |> nansum
    mom1 = ((waves .* deltawave .* flux)  |> nansum) ./ (flux |> nansum)
    mom2 = sqrt(nansum((waves .- mom1).^2 .* flux) / nansum(flux))# |> nansum
end


"""    make_lines_mask(mod; window_width_kms=1000, plot_it=false)
Creates a mask which keeps data in a window around each line in the model, and discards
the rest. The window width is customizable. Mainly used internally in the package.
## Input
- `mod::GModelFit.Model`: The model with the lines already defined.
## Optional input
- `window_width::Float64`: Half width in km/s of the fitting window around each line.
"""
function make_lines_mask(mod::Model, dom::Domain; window_width_kms::Float64=1000., plot_it::Bool=false, plot_comps::Bool=false)
    wave = coords(dom)
    mask = falses(size(wave))
    for comp in keys(mod)
        if comp == :main
            continue
        end
        cenwave = mod[comp].lab_wave.val * (1 + mod[comp].redshift.val)
        window_width_ang = v_to_deltawl(window_width_kms, cenwave)
        mask[cenwave-window_width_ang .< wave .< cenwave+window_width_ang] .= 1
    end 
    if plot_it == true
        theexponent = log10.(mod()) |> maximum |> floor
        thefactor = 10^theexponent
        Plots.plot(coords(mod.domain), mask .* 0.5 .* maximum(mod())/thefactor, label="Mask")
        if plot_comps
            for comp in keys(mod.buffers)
                comp == :main && continue
                Plots.plot!(coords(mod.domain), mod.buffers[comp]./thefactor, label=string(comp))
            end 
        end
        Plots.plot!(coords(mod.domain), mod()./thefactor, title="Model and mask", label="Model")
        Plots.xlabel!("Wavelength")
        Plots.gui()
    end 
    return mask
end 


"""    nanmask10(flux, dflux)
Removes wave/velocity bins in which wither spec or errs are NaN, or where either is == 0.
The `flux` and `dflux` (error) arrays must have the same size. Internal function, not meant
to be user facing.
# Arguments:
- `flux::Array{Float64}` - the flux (flux density) array.
- `dflux::Array{Float64}` - the standard errors to the flux
# Returns:
- A Tuple of the same two arrays, with the `NaN` values replaced as described above.
"""
function nanmask10(
    flux::Array, dflux::Array)::BitVector
    idx1 = (.!isnan.(flux)) .& (.!isnan.(dflux)) .& (dflux .!= 0.)
    idx = idx1
    return idx
end 


"""    gauss(λ::Array, σ::Float64, μ::Float64)
Return a simple normalized Gaussian with an integral of 1.
"""
function gauss(λ::Array,  σ::Float64, μ::Float64)
    g = @. 1/σ/sqrt(2π) * exp(-(λ - μ)^2/2/σ^2)
    return g
end 


"""    estimate_line_snr(wave, flux; err=nothing)

Quick and coarse estimate of the integrated flux within a wavelength window
without actual fitting.
"""
function estimate_line_snr(wave, flux; err=nothing, only_flux=false)
    if err isa Nothing
        median_fit = running_median(flux, 3, :sym, nan=:ignore)
        nuisance_be_gone = flux .- median_fit
        stddev = nanstd(nuisance_be_gone) .* ones(size(flux))
    else
        stddev = err
    end 
    stddev = abs.(stddev)
    uflux = flux .± stddev
    usum = nansum(uflux)
    snr = Measurements.value(usum)/Measurements.uncertainty(usum)
    if isnan(snr); snr = 0.1; flux=[0]; end
    if flux |> maximum |> isnan; outflux = 0; else; outflux = flux |> maximum; end
    if only_flux; return Measurements.value(usum); end
    return snr
end 


"""    fill_in_fit_values(dict::Dict, fitresults::Dict, pix_coords::Tuple{Int,Int}; fitstats=NaN)
"""
function fill_in_fit_values!(
    dict::Dict{Any,Any}, fitdict::Dict; fitstats=NaN, pix_coords=nothing, mask=nothing)

    if mask isa Nothing
        mask = falses(dict[:fwhm] |> size)
    end

    fitstats = fitdict[:fitstats].fitstat
    fitresults = fitdict[:results]
    @debug dict[:refline] fitresults[Symbol(dict[:refline])] 
    dict[:fitstats][mask] .= fitstats
    dict[:redshift][mask, 1] .= fitresults[Symbol(dict[:refline])].redshift.val
    dict[:redshift][mask, 2] .= fitresults[Symbol(dict[:refline])].redshift.unc
    dict[:redshift][mask, 3] .= fitresults[Symbol(dict[:refline])].redshift.val / 
        fitresults[Symbol(dict[:refline])].redshift.unc
    dict[:fwhm][mask, 1] .= fitresults[Symbol(dict[:refline])].fwhm_kms.val
    dict[:fwhm][mask, 2] .= fitresults[Symbol(dict[:refline])].fwhm_kms.unc
    dict[:fwhm][mask, 3] .= fitresults[Symbol(dict[:refline])].fwhm_kms.val /
        fitresults[Symbol(dict[:refline])].fwhm_kms.unc
    if :redshift_broad in keys(dict)
        dict[:redshift_broad][mask, 1] .= fitresults[Symbol(dict[:refline_broad])].redshift.val
        dict[:redshift_broad][mask, 2] .= fitresults[Symbol(dict[:refline_broad])].redshift.unc
        dict[:redshift_broad][mask, 3] .= fitresults[Symbol(dict[:refline_broad])].redshift.val /
            fitresults[Symbol(dict[:refline_broad])].redshift.unc
        dict[:fwhm_broad][mask, 1] .= fitresults[Symbol(dict[:refline_broad])].fwhm_kms.val
        dict[:fwhm_broad][mask, 2] .= fitresults[Symbol(dict[:refline_broad])].fwhm_kms.unc
        dict[:fwhm_broad][mask, 3] .= fitresults[Symbol(dict[:refline_broad])].fwhm_kms.val /
            fitresults[Symbol(dict[:refline_broad])].fwhm_kms.unc
    end 
    fwhmu = dict[:fwhm][mask, 1] .± dict[:fwhm][mask, 2]
    sigmau = fwhm_to_sigma(fwhmu)
    ## Fill in the fitted flux of each line in each spaxel
    for l in keys(dict)
        if !(l in keys(fitdict[:linesdict]))
            continue
        end
        component = fitdict[:linesdict][l]
        normu = component[:norm].val ± component[:norm].unc
        fluxu = sqrt(2π) * sigmau * normu
        # println("$pix_coords.val, $fluxu.val, $normu.val, $(fluxu.val/normu.val)")
        # fluxu = normu
        # dict[l][mask, 1] = fluxu.val
        # dict[l][mask, 2] = fluxu.err
        # dict[l][mask, 3] = fluxu.val/fluxu.err
        dict[l][mask, 1] .= normu.val
        dict[l][mask, 2] .= normu.err
        dict[l][mask, 3] .= normu.val / normu.err
    end
end

###=============================================================================
#    File I/O functions

"""    write_spectral_cube_to_fits(filepath, spectralcube)
Write a subtype of AbstractSpectralCube to a FITS file.
# Arguments
- `filepath::String`: Path to save the cube to.
- `spectralcube::AbstractSpectralCube`: The datacube to save.
"""
function write_spectral_cube_to_fits(filepath::String, spectralcube::NIRSpecCube)
    # TODO: Make possible to choose units. 
    # TODO: Make sure units-related FITS keywords are correct.
    f = FITS(filepath, "w")
    primhead = spectralcube.primheader; primhead["CONTINUE"] = ""
    header = spectralcube.header; header["CONTINUE"] = ""
    flux = spectralcube.fluxcube .* spectralcube.flux_convert
    errs = spectralcube.errscube .* spectralcube.flux_convert
    write(f, Float64[], header=primhead)
    write(f, flux, header=header, name="SCI")
    write(f, errs, header=header, name="ERR")
    close(f)
end


# """    write_maps_to_fits(filepath, mapsdict)
# Write a dictionary of line maps as output by the cube fitter function to a multi-
# extension FITS file.
# """
# function write_maps_to_fits(filepath::String, mapsdict::Dict)
#     f = FITS(filepath, "w")
#     keylist = sort(collect(keys(mapsdict)))
#     primhead = mapsdict[:header]; primhead["CONTINUE"] = ""  #push!(primhead, "CONTINUE"="")
#     sechead = mapsdict[:header]; sechead["CONTINUE"] = ""
#     write(f, Float64[], header=primhead)
#     # for k in keys(mapsdict)
#     for k in keylist
#         if k in [:header, :primhead, :refline]
#             continue
#         end 
#         write(f, mapsdict[k], header=sechead, name=String(k))
#         # write(f, mapsdict[k],  name=String(k))
#     end 
#     close(f)
#     return nothing
# end 

"""    write_maps_to_fits(filepath, mapsdict)
Write a dictionary of line maps as output by the cube fitter function to a multi-
extension FITS file. NEW VERSION
# Arguments
- `filepath::String`: Path to save the file to.
- `mapsdict::Dict`: Dictionary of the type that is output from `fit_cube`.
# Returns
- `nothing`
"""
function write_maps_to_fits(filepath::String, mapsdict::Dict)
  f = FITS(filepath, "w")
  keylist = sort(collect(keys(mapsdict)))
  primhead = mapsdict[:header]; primhead["CONTINUE"] = ""
  sechead = mapsdict[:header]; sechead["CONTINUE"] = ""
  write(f, collect(0:0.01:1), header=primhead)
  for k in keylist
    if k in [:header, :primhead, :refline, :refline_broad]
        continue
    end 
    @debug "What is this nonsense?" String(k) mapsdict[k]
    write(f, mapsdict[k], header=sechead, name=String(k))
  end 
  close(f)
  return nothing
end 


"""    load_slices_dict(filepath::String)
Load in a FITS file of "slices" saved by the `write_maps_to_fits` function.
Convenient e.g. for inspection.
# Arguments
- `filepath::String`: Path to the file to open.
"""
function load_slices_dict(filepath::String)
  f = FITS(filepath, "r", extendedparser=true)
  outdict = Dict()
  for (i, hdu) in enumerate(f)
    try
      extname = read_key(hdu, "EXTNAME")[1]
      println(i-1, " ", Symbol(extname))
      outdict[Symbol(extname)] = read(hdu)
    catch
      outdict[:PRIMARY] = read_header(hdu)
      println("No EXT name for ", i-1)
    end
  end
  return outdict
end


"""    quicklook_slice(slicedict, name; what="data", norm="sqrt", colorlimits=nothing)

Quick and convenient visualization of slices output fromthe `fit_cube()` function.

# Parameters
- `slicedict::Dict`: The dictionary containing the data slices
- `name::Symbol`: Key for the slice to visualize
# Optional parameters
- `norm`: Function to normalize/color scale the data. Must be callable with a single 
  argument. E.g., `identity`, `sqrt`, `log10`, or similar.
- `what::Symbol`: One of wither `:data`, `:errs`, `:snr`, `:numdata`, or `:numerrs`.
- `cmap::Symbol`: The colormap to use with `Plots.heatmap()`.
- `colorlimits::Tuple`: Tuple of (min, max) color cut values.
"""
function quicklook_slice(slicedict, name; what=:data, norm=identity, colorlimits=nothing, cmap=:cubehelix, comment="")
    theslice = slicedict[Symbol(name)]
    layers = Dict(:data => 1, :errs => 2, :snr => 3, :numdata => 4, :numerrs => 5)
    # thelayer = layers[Symbol(what)]
    data = transpose(theslice[:,:,layers[Symbol(what)]])
    if colorlimits isa Nothing
        colorlimits = (nanpctile(data, 2), nanpctile(data, 99))
    end
    data = clamp.(data, colorlimits[1], colorlimits[2])
    Plots.heatmap(data .|> norm, aspect_ratio=:equal, color=cmap)
    Plots.title!("Showing: $name; kind: $what, $comment")
end 


"""    quicklook_fit_result_dict(indict)
Plots of the results of a fit of a single spectrum. 
# Parameters
- `indict::Dict`: The output dictionary of `fit_spectrum_from_subcube`
"""
function quicklook_fit_result_dict(indict)
    Plots.plot(indict[:wave], indict[:spec], color=:gray)
    Plots.plot!(indict[:wave], indict[:errs])
    inmodel = indict[:results]
    Plots.plot!(inmodel.domain.axis[1], inmodel(), linewidth=2)
    for b in inmodel.buffers |> keys
        Plots.plot!(inmodel.domain.axis[1], inmodel.buffers[b])
    end
    Plots.hline!([0], color=:black, linewidth=1,)
end


"""    quicklook_model(inmodel)
# Parameters
- `inmodel::`: GModelFit.Model
"""
function quicklook_model(inmodel, cube)
    dom = Domain(cube.wave)
    Plots.plot(dom.axis[1], inmodel(dom), linewidth=2)
    # for b in inmodel.buffers |> keys
    # for b in inmodel.comps|> keys
    #     Plots.plot!(dom.axis[1], inmodel.comps[b](dom))
    # end
    Plots.hline!([0], color=:black, linewidth=1,)
end


###=============================================================================
#   Segmentation map and Voronoi binning related helper functions
#
"""    mean_by_segmap(data, segmap; errs=nothing)

Given a 2D array and a map of spatial bins of the same size, outputs a new array 
with each bin filled with the mean of the data witin it. Is NaN/missing-aware, and 
optionally takes an array of standard errors which is propagated in the standard way.

For bins of 1 pixel size, this is just the identity operator. 
"""
function mean_by_segmap(data, segmap; errs=nothing)
    if !(errs isa Nothing)
        data = data .±ᵤ errs
        # data = data .± errs
    end
    outmap = similar(data)
    for i in (segmap |> unique |> sort)
        segmask = segmap .== i
        spaxs = count(segmask)
        print("`mean_by_segmap` processing segment # $i with $spaxs spaxels \r")
        dpoints = vec(data[segmask])
        outmap[segmap.==i] .= nanmean(dpoints)
    end
    return outmap .|> Measurements.Measurement
end


###=============================================================================
#   Placeholder functions for extensions
#
""" Please load the `VoronoiBinning` package to activate this functionality.
"""
function voronoi_bin_slice(::Any)
  println("This function requires the `VoronoiBinning` module to be installed and loaded.")
  return nothing
end

export voronoi_bin_slice, mean_by_segmap

end  # End module
