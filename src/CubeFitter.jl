module CubeFitter
###============###

export AbstractSpectralCube, NIRSpecCube, MUSECube
export calculate_moments, fit_cube, fit_spectrum_from_subcube
export make_spectrum_from_cutout, make_lines_mask, toggle_fnu_flam

using Measurements: result
###===========================================###

# Basic computing functionality, misc.
using Base: available_text_colors_docstring, NullLogger
using Base.Threads, Printf, Logging, LoggingExtras
# Handle NaN's more gracefully than standard Julia
using NaNStatistics
# Physical units and unit conversion, and uncertainties
# using Unitful, UnitfulAstro, UnitfulEquivalences, Measures, PhysicalConstants.CODATA2018
using Unitful, UnitfulAstro, UnitfulEquivalences, Measurements, PhysicalConstants.CODATA2018
# Data formats, table, the works
using FITSIO, DataFrames, DataStructures, CSV
# Interpolation, modeling and fitting
using Polynomials, GModelFit, FastRunningMedian
import GModelFit as gmf
import Interpolations: LinearInterpolation
# Numerical integration
using NumericalIntegration
import NumericalIntegration as nui
using Term.Progress
include("./SpecHelpers.jl")
using .SpecHelpers

# Physical constants
const ckms = SpeedOfLightInVacuum |> u"km/s"
const caps = SpeedOfLightInVacuum |> u"Å/s"

# datapath = "./static_data/"
# datapath = joinpath(@__DIR__, "..", "static_data")
# println(@__MODULE__)
datapath = joinpath(dirname(pathof(CubeFitter)), "..", "static_data")
# include("./SpectralCubes.jl")
# using .SpectralCubes


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
    # if level == :info
    #     #@info "Setting logging level to info-only"
    #     global_logger(logger(filename, Logging.Info))
    #     # global_logger(ConsoleLogger(stderr, Logging.Info))
    # elseif level == :none
    #     #@info "Turning off logging"
    #     global_logger(NullLogger());
    # elseif level == :warning
    #     #@info "Including info and warnings in feedback"
    #     global_logger(logger(filename, Logging.Warn))
    # elseif level == :debug
    #     #@info "Turning on full, debuggin-level logging"
    #     global_logger(logger(filename, Logging.Debug))
    # else
    #     println("""Please chose a valid level.""")
    # end 
    if level == :info
        #@info "Setting logging level to info-only"
        global_logger(ConsoleLogger(stderr, Logging.Info))
    elseif level == :none
        #@info "Turning off logging"
        global_logger(NullLogger());
    elseif level == :warning
        #@info "Including info and warnings in feedback"
        global_logger(ConsoleLogger(stderr, Logging.Warn))
    elseif level == :debug
        #@info "Turning on full, debuggin-level logging"
        global_logger(ConsoleLogger(stderr, Logging.Debug))
    else
        println("""Please chose a valid level.""")
    end 
end 

# global_logger(ConsoleLogger(stderr, Logging.Info))  # Default
global_logger(NullLogger())
# global_logger(FileLogger("CubeFitter.log"))  # Default


# include("SpectralCubes.jl")

# Declaring an abstract SpectralCube type will allow us to make SpectralCube objects
# inheriting from it for various ls instruments, which require different treatments and
# different kind of, and still be able to write functions which can operate on all of them. 
abstract type AbstractSpectralCube end

mutable struct GenericSpectralCube <: AbstractSpectralCube
    instrument::String
    wave::Array
    fluxcube::Array
    errscube::Array
end

mutable struct NIRSpecCube <: AbstractSpectralCube
    """ A struct representing a NIRSpec IFU datacube"""
    grating::String
    wave::Array
    fluxcube::Array
    errscube::Array
    header::FITSHeader
    primheader::FITSHeader
    linelist::DataFrame
    lsf_fitter
    z_init
    ref_line
    function NIRSpecCube(filepath::String, grating::String;
        # linelist_path="./static_data/neblines.dat",
        linelist_path=joinpath(datapath, "neblines.dat"),
        # lsf_file_path="./static_data/jwst_nirspec_$(grating)_disp.fits",
        lsf_file_path=joinpath(datapath, "jwst_nirspec_$(grating)_disp.fits"),
        z_init=0, reference_line=:OIII_5007)

        linelist = load_neblines(linelist_path)
        ddict = load_fits_cube(filepath)
        ddict = convert_ergscms_Å_units(ddict, "NIRSpec")
        wave = ustrip.(ddict[:Wave])
        fluxcube = ustrip.(ddict[:Data])
        errscube = ustrip.(ddict[:Errs])
        header = ddict[:Header]
        primheader = ddict[:Primheader]
        itp = get_resolving_power("NIRSpec", setting=grating)
        new(grating, wave, fluxcube, errscube, header, primheader, linelist, itp, z_init, reference_line)
    end
end


"""    MUSECube(filepath[, linelist_path="])
"""
mutable struct MUSECube <: AbstractSpectralCube
    setting::String
    resol_fit_degree::Int
    header::FITSIO.FITSHeader
    primheader::FITSIO.FITSHeader
    linelist::DataFrame
    lsf_fitter
    function MUSECube(filepath::String;
                      linelist_path=:)
        ddict = load_fits_cube(filepath, "r")
    end
end


function estimate_redshift_reference_line_flux(spec_cube; xrange=nothing, yrange=nothing)
    wave = spec_cube.wave
    spec_guess, _ = make_spectrum_from_cutout(spec_cube, xrange, yrange)
    neblines = spec_cube.linelist
    λ_obs = neblines[neblines.name .== string(spec_cube.ref_line), :lamvac][1] * (1 +spec_cube.z_init)
    iha_g = λ_obs - 15 .< spec_cube.wave .< λ_obs + 15
    fha_g = abs(nui.integrate(wave[iha_g], spec_guess[iha_g]))
    reds_g = (nansum(wave[iha_g] .* spec_guess[iha_g]) ./ nansum(spec_guess[iha_g])
        ./  neblines[neblines.name .== string(spec_cube.ref_line), :lamvac][1] - 1)
    @debug "Flux guess" fha_g reds_g
    return fha_g, reds_g
end 



"""    load_fits_cube(inpath)
Convenience function to load a FITS cube and return a Dict of various, sometimes useful,
quantities read or derived from the cube.
## Returns
- Dict containing wave array, data and error cube, and headers from the primary HDU and
  first extension from the data file.
"""
function load_fits_cube(inpath::String)::Dict{Symbol, Any}
    fitsfile = FITS(inpath, "r")
    header = read_header(fitsfile[2])
    primary = read_header(fitsfile[1])
    data = read(fitsfile[2])
    errs = read(fitsfile[3])
    naxis3 = header["NAXIS3"]
    crval3 = header["CRVAL3"]
    if haskey(header, "CDELT3")
        cdelt3 = header["CDELT3"]
    elseif haskey(header, "CD3_3")
        cdelt3 = header["CD3_3"]
    else
        println("Could not find spectral axis keywords in FITS header")
    end
    waves = crval3 .+ cdelt3 .* (0:naxis3-1)
    out = Dict(:Wave => waves, :Data => data, :Errs => errs,
               :Header => header, :Primheader => primary)
    return out
end


"""    get_resolving_power(instrument[, setting=nothing])
This function figures out which instrument we are working with, and calls the relevant
instrument-specific function
"""
function get_resolving_power(instrument::String; setting=nothing::Any)
    if lowercase(instrument) == "nirspec"
        return(get_resolving_power_nirspec(setting))
    elseif lowercase(instrument) == "muse"
        return get_resolving_power_muse(degree=setting)
    else
    end
end


"""     get_resolving_power_nirspec(grating[, calib_path="../static_data/"])
This function assumes that the FITS files containing the NIRSpec grating resolving power
are named "jwst_nirspec_[grating]_disp.fits", as fetched from the official JDox.
## Returns
- A callable `Interpolations.Extrapolation` objects which returns the linearly interpolated
  resolving power at the wavelength passed to it.
"""
function get_resolving_power_nirspec(grating::String; calib_path=datapath)  #"./static_data/")
    lsffilename = "jwst_nirspec_$(grating)_disp.fits"
    @debug "Path to LSF file: " calib_path * lsffilename
    fitsfile = FITS(joinpath(calib_path, lsffilename))
    angwave = read(fitsfile[2], "WAVELENGTH") * u"μm" .|> u"angstrom"
    R = read(fitsfile[2], "R")
    @debug extrema(angwave) extrema(R)
    itp = LinearInterpolation(ustrip.(angwave), R)
    return itp
end


function convert_ergscms_Å_units(datadict::Dict, instrument="NIRSPec")
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



"""    toggle_fnu_flam()
Converts fnu data to flam units, and vice versa.
Input must be Unitful™ (i.e., Quantities).
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


"""    make_spectrum_from_cutout()
Converts fnu data to flam units, and vice versa.
Input must be Unitful™ (i.e., Quantities).
"""
function make_spectrum_from_cutout(cube, xrange=nothing, yrange=nothing)
    if isa(xrange, Nothing); xrange = (1:size(cube.fluxcube)[1]); end
    if isa(yrange, Nothing); yrange = (1:size(cube.fluxcube)[2]); end
    cutout_spec = cube.fluxcube[xrange, yrange, :]
    cutout_errs = cube.errscube[xrange, yrange, :]
    slice = cutout_spec[:,:,1]
    if length(slice) == 1
        spec_init = vec(cutout_spec)
        errs_init = vec(cutout_errs)
    else
        spec_init = nanmean(cutout_spec, dim=(1,2))
        errs_init = sqrt.(nansum(cutout_errs.^2, dim=(1,2))) / length(slice)
    end
    # @debug "`make_spectrum_from_cutout()`: Number of spatial pix in cutout" size(slice)
    return spec_init, errs_init
end



"""    fit_subcube(cube[; xrange=nothing, yrange=nothing])
Only returns one fit result; if a region of many pixels are given, they are first combined
to one spectrum; then fitted.
## Returns
- `Dict`: A dict containing wave, spectrum, errors, fit result and fit statistics.
"""
function fit_spectrum_from_subcube(cube; xrange=nothing, yrange=nothing)
    wave_init = cube.wave
    spec, errs = make_spectrum_from_cutout(cube, xrange, yrange)
    @debug "Before NaNmask10'ing: " count(isnan.(spec)) count(isnan.(errs))
    notnan = nanmask10(spec, errs)
    goodwave, goodspec, gooderrs = wave_init[notnan], spec[notnan], errs[notnan]
    # @debug xrange yrange goodwave, count(notnan)
    @debug "Goodwave: " goodwave
    mod = build_model(cube, xrange=xrange, yrange=yrange, dom_init=Domain(goodwave))
    @debug "First pass at making model done!"
    idx = make_lines_mask(mod)
    @debug "Made lines mask for $xrange, $yrange !"
    lwave, lspec, lerrs = goodwave[idx], goodspec[idx], gooderrs[idx]
    @debug "Before and after line masking: " length(goodwave), length(lwave)
    if length(lwave) == 0; return NaN; end
    lmod = build_model(cube, xrange=xrange, yrange=yrange, dom_init=Domain(lwave))
    lmeas = Measures(Domain(lwave), lspec, lerrs)
    @debug "Model before fitting: $xrange, $yrange:  " lmod lmeas lwave
    if !(cube.ref_line in keys(lmod)); return NaN; end
    try  # Too many things can go wrong to catch eatch one separately.
        results, stats = gmf.fit(lmod, lmeas)
        #@info "Successfully fitted (col. $xrange, row $yrange)"
        @debug results stats
        output = Dict(
            :wave => goodwave,
            :spec => goodspec,
            :errs => gooderrs,
            :measure => lmeas,
            :results => results,
            :fitstats => stats)
        return output 
    catch
        #@warn "Fitting (col. $xrange, row $yrange) failed! Moving on"
        return NaN
    end 
end



function makeresultdict(cube)
    @debug "Making output dict"
    xsize, ysize = size(cube.fluxcube)[1], size(cube.fluxcube)[2]
    neblines = cube.linelist
    wave = cube.wave
    # Create dictionary to hold the fit results.
    slices_dict = Dict()
    slices_dict[:refline] = cube.ref_line
    slices_dict[:header] = cube.header
    # slices_dict[:primhead] = cube.primheader
    # slices_dict[:primhead]["REFLINE"] = cube.ref_line
    # set_comment!(slices_dict[:primhead], "REFLINE", "Reference line used for line fitting.")
    for scomp in neblines.name
        lamvac = neblines[neblines.name .== scomp, :lamvac][1]
        obsvac = lamvac * (1 + cube.z_init)
        if (minimum(ustrip.(wave)) > obsvac) | (maximum(ustrip.(wave)) < obsvac)
            @debug "Line $scomp outside  wavelength range, skipping"
            continue
        end 
        comp = Symbol(scomp)
        if comp == :main; continue; end
        slice = zeros(xsize, ysize, 3) .* NaN
        slices_dict[comp] = slice
    end 
    slices_dict[:redshift] = zeros(xsize, ysize, 3) .* NaN
    slices_dict[:fwhm] = zeros(xsize, ysize, 3) .* NaN
    slices_dict[:mom0] = zeros(xsize, ysize, 3) .* NaN
    slices_dict[:mom1] = zeros(xsize, ysize, 3) .* NaN
    slices_dict[:mom2] = zeros(xsize, ysize, 3) .* NaN
    slices_dict[:fitstats] = zeros(xsize, ysize) .* NaN
    @debug "Output dict successfully done!"
return slices_dict
end


function fit_cube(cube)
    xsize, ysize = size(cube.fluxcube)[1], size(cube.fluxcube)[2]
    slices_dict = makeresultdict(cube)
    @debug "Made slices_dict"
    @track for x in 1:xsize
        for y in 1:ysize
            print("Fitting col. $x, row $y \r")
            # numnan = size((collect(isnan.(cube.fluxcube[x, y, :]))))
            numnan = count(isnan.(cube.fluxcube[x, y, :]))
            nanratio = numnan/length(cube.wave)
            @debug "Number of nans in x=$x, y=$y: " numnan size(cube.fluxcube) nanratio
            if nanratio > 0.5
                #@warn "Pixel ($x, $y) has > 50% NaN; moving on"
                continue
            end
            fitdict = fit_spectrum_from_subcube(cube, xrange=(x:x), yrange=(y:y))
            if !isa(fitdict, Dict)
                @debug "Fitdict $x, $y was not a dict!"
                continue
            end
            # @debug x, y, typeof(x), typeof(y) typeof(fitdict[:results]) typeof(fitdict[:fitstats])
            @debug fitdict[:results] fitdict[:fitstats]
            fill_in_fit_values!(slices_dict, fitdict, (x, y))
            sleep(0.001)
        end 
    end 
    mom0, mom1, mom2 = calculate_moments(cube)
    slices_dict[:mom0] = mom0
    slices_dict[:mom0] = mom1
    slices_dict[:mom0] = mom2
return slices_dict
end 



# TODO: Make sure docstring is up to date! 
# TODO: Could perhaps be `setting` for generality? 🤔
"""    build_model(cube[; xrange=nothing, yrange=nothing, min_snr=1.5, fwhm_int=100, dom_init=nothing])
Builds the `GModelFit.Model` object necessary for fitting the lines.
The model assumes that each line is described as a single Gaussian profile, with a shared
`redshift` and velocity width `fwhm`, while their individual fluxes are left as free parameters
(including lines which actually have fixed line ratios, which might not have that in the ). 
## Input
- `domain::GModelFit.Domain{1}`: A Domain object containing the wavelength range on which
  the model should be defined.
- `neblines::DataFrame`: The dataframe containing results from `load_neblines`.
- `redsh::Float64`: An initial guess for the redshift of the object we're looking at.
## Optional input
- `fwhm_int::Float64`: The initial guess for the physical (*not* observed) FWHM of the
  line, in units of km/s. Default = 100.
- `instrument::String`: Either 'NIRSpec' or 'MUSE' for the time being. Default is 'NIRSpec'.
- `order::Int`: Order of the polynomial to fit. Only relevant if `instrument`='MUSE'.
  Defaults to 3.
- `grism::String`: Should actually be `grating`. Only relevant for 'NIRSpec', at the time
  being. Default value is 'g140h'. 
## Output
- `model::GModelFit.Model`: The model containing the prescribed lines. The Reference line has
  norm=1., and the others are scaled according to the values of `foverha` given in the line list
  input file.
"""
function build_model(cube; xrange=nothing, yrange=nothing, min_snr=1.5, fwhm_int=100, dom_init=nothing)
    redss = (1. + cube.z_init)  # Just for convenience
    wave = cube.wave
    neblines = cube.linelist
    if isa(dom_init, Nothing); dom_init = Domain(wave); end
    @debug "buld_model() domain: " xrange yrange dom_init
    spec_init, errs_init = make_spectrum_from_cutout(cube, xrange, yrange)
    components = OrderedDict()
    complist = Vector{Symbol}()
    λobs_ref = neblines[neblines.name .== string(cube.ref_line), :lamvac][1] * (redss)
    for l in neblines.name[1:end]
        lamvac  = neblines[neblines.name.==l, :lamvac][1]
        lam_obs = lamvac * redss
        # if (lam_obs > maximum(coords(dom_init)) - 10) | (lam_obs < minimum(coords(dom_init)) + 10)
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
            cube, xrange=xrange, yrange=yrange)
        if isnan(fha_g) || isnan(reds_g)
            @debug "Numerical integration of line $l failed due to lack of data." l
            continue
        end
        @debug "Results of numerical integration" fha_g
        foverha = neblines[neblines.name.==l, :foverha][1]
        @debug l neblines.name
        norm = fha_g * foverha
        fwhm_inst = cube.lsf_fitter(λobs_ref)
        # @debug "Resolving power at line $l: " fwhm_inst
        fwhm_tot = sqrt(fwhm_int^2 + fwhm_int^2)
        @assert isa(fwhm_tot, Number)
        @assert isa(lamvac, Number)
        @assert isa(cube.z_init, Number)
        @assert isa(foverha, Number)
        @assert isa(norm, Number) "$foverha"
        @assert isa(lamvac, Number)
        g = GModelFit.FComp(
            _gauss_line, [:waves], lab_wave=lamvac, redshift=cube.z_init,
            fwhm_kms=fwhm_tot, norm=norm, lsf=fwhm_inst)
        components[Symbol(l)] = g # = ref_line
        append!(complist, [Symbol(l)])
        @debug "Added the line $l to the model"
    end 

    # First, define model with mutually independent components
    components[:main] = SumReducer(complist)
    mod = Model(dom_init, components)

    for comp in keys(mod)
        @debug "Setting up $comp"
        # INFO: If a component is patched to itself, it in essence gets fixed!!
        # So definitely do not let the reference line get `patch`'ed! 
        # TODO: Find a more elegant solution for the `if` statement below.
        if comp == :main
            @debug "Comp was `:main`, moving on"
            continue
        elseif comp != Symbol(cube.ref_line)
            @debug "Fixing parameters for $comp"
            mod[comp].lab_wave.fixed = true
            mod[comp].redshift.low = cube.z_init - cube.z_init * 0.1
            mod[comp].redshift.high = cube.z_init + cube.z_init * 0.1
            mod[comp].redshift.patch = Symbol(cube.ref_line)
            mod[comp].fwhm_kms.patch = Symbol(cube.ref_line)
            mod[comp].fwhm_kms.low = 1.
            mod[comp].fwhm_kms.high = 2000.
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
        end 
    end 
    @debug "`build_model()` successfully done!"
    return mod  #, components
end 



"""    _gauss_line(waves, lab_wave, redshift, fwhm_kms, norm, lsf)
"""
function _gauss_line(
    waves::Array, lab_wave::Float64, redshift::Float64, fwhm_kms::Float64, norm::Float64, lsf::Float64)  #, grating="g140h")
    cen_wave = lab_wave * (1. + redshift)
    fwhm_inst = cen_wave / lsf  # get_resolving_power(
        # cen_wave, instrument=instrument, order=order, grism=grating)
    fwhm_aa = v_to_deltawl(fwhm_kms, cen_wave)
    #= fwhm_tot = sqrt(fwhm_inst^2 + fwhm_kms^2) =#
    fwhm_tot = sqrt(fwhm_inst^2 + fwhm_aa^2)
    sigma_aa = fwhm_to_sigma(fwhm_tot)
    #= sigma_aa = fwhm_tot / ckms * cen_wave / (2*sqrt(2*log10(2))) =#
    flux = gauss(waves, ustrip(sigma_aa), ustrip(cen_wave)) .* norm
    return flux
end
#= function _gauss_line( =#
#=     waves::Array, lab_wave::Float64, redshift::Float64, fwhm_kms::Float64, norm::Float64, lsf::Float64)  #, grating="g140h") =#
#=     cen_wave = lab_wave * (1. + redshift) =#
#=     fwhm_inst = ustrip(ckms) / lsf  # get_resolving_power( =#
#=         # cen_wave, instrument=instrument, order=order, grism=grating) =#
#=     fwhm_tot = sqrt(fwhm_inst^2 + fwhm_kms^2) =#
#=     sigma_aa = fwhm_tot / ckms * cen_wave / (2*sqrt(2*log10(2))) =#
#=     flux = gauss(waves, ustrip(sigma_aa), ustrip(cen_wave)) .* norm =#
#=     return flux =#
#= end =#


 
"""
Computes moments of a flux slab with associated wavelengths.
Computes the 0th, 1st and 2nd moments of a spectral chunk (a line map) in one or more
pixels.
## Input:
- `waves::Array{Float64, 1}`: Wavelength array, must be 1-dimensionsl (a vector)
- `flux::Array{Float64, 2}`: The flux (or rather flux density) data, must be 3-dimensional.
Errors/uncertainties are not included in this computation
## Returns
- A tuple of three 2D arrays (images) of the moments mapped to each spatial pixel.
"""
# function calculate_moments(waves::Array{Float64, 1}, flux::Array{Float64, 3})
function calculate_moments(cube; refline=nothing, window_kms=1000)
    if isa(refline, Nothing)
        refline = cube.ref_line
    end 
    λ_obs = cube.linelist[cube.linelist.name .|> String .== refline |> String, :lamvac][1] * (1 + cube.z_init)
    window_ang = v_to_deltawl(window_kms, λ_obs)
    idx = λ_obs - window_ang .< cube.wave .< λ_obs + window_ang
    waves, flux = cube.wave[idx], cube.fluxcube[:,:,idx]
    println(size(flux))
    # Make sure things are well formatted now before they get passed along.
    @assert size(waves)[1] == size(flux)[3]
    # Empty arrays to fill in:
    @debug waves
    mom0 = similar(flux[:, :, 1])
    mom1 = similar(flux[:, :, 1])
    mom2 = similar(flux[:, :, 1])
    ncols, nrows = size(flux)
    @track for i in 1:ncols
        for j in 1:nrows
            theflux = flux[i, j, :]  # Just for convenience
            mom0a = NumericalIntegration.integrate(theflux, waves)
            @debug "Moment0: " mom0a
            if mom0a <= 0; mom0a = NaN; end
            mom0[i, j] = mom0a
            mom1[i, j] = nansum(waves .* theflux) ./ nansum(theflux)
            @debug i, j theflux mom0[i, j] mom1[i, j]
            if isnan.([mom0[i, j], mom1[i, j]]) |> any
                mom2[i, j] = NaN
            else 
    # mom2 = np.sqrt(np.nansum((wave - mom1)**2 * flux, axis=0) / np.absolute(np.nansum(flux, axis=0)))
                mom2[i, j] = sqrt(abs.(
                sum((waves .- mom1[i, j]).^2 .* theflux))
                / abs.(nansum(theflux)))# - λ_obs
            end 
        end
    end 
    return mom0, mom1, mom2
end 



"""    make_lines_mask(mod[; window_width_kms, plot_it])
Creates a mask which keeps data in a window around each line in the model, and discards
the rest. The window width is customizable.
## Input
- `mod::GModelFit.Model`: The model with the lines already defined.
## Optional input
- `window_width::Float64`: Half width in km/s of the fitting window around each line.
"""
function make_lines_mask(mod::Model; window_width_kms::Float64=1000., plot_it::Bool=false)
    wave = coords(mod.domain)
    mask = falses(size(wave))
    for comp in keys(mod)
        if comp == :main
            continue
        end
        cenwave = mod[comp].lab_wave.val * (1 + mod[comp].redshift.val)
        window_width_ang = v_to_deltawl(window_width_kms, cenwave)
        @debug "Make line mask: " wave window_width_ang
        mask[cenwave-window_width_ang .< wave .< cenwave+window_width_ang] .= 1
    end 
    if plot_it == true
        plt.plot(coords(mod.domain), mod())
        for comp in keys(mod.buffers)
            comp == :main && continue
            plt.plot(coords(mod.domain), mod.buffers[comp])
        end 
        plt.plot(coords(mod.domain), mask .* 0.5 .* maximum(mod()))
    end 
    return mask
end 



"""    nanmask10(flux, dflux)
Removes wave/velocity bins in which wither spec or errs are NaN, or where either is == 0.
The `flux` and `dflux` (error) arrays must have the same size.
## Input:
- `flux::Array{Float64}` - the flux (flux density) array.
- `dflux::Array{Float64}` - the standard errors to the flux
## Returns:
- A Tuple of the same two arrays, with the `NaN` values replaced as described above.
"""
function nanmask10(
    flux::Array{Float64}, dflux::Array{Float64})::BitVector
    idx1 = (.!isnan.(flux)) .& (.!isnan.(dflux)) .& (dflux .!= 0.)
    idx = idx1# .& idx2
    return idx  #flux[idx], dflux[idx]
end 



"""    gauss(λ::Array, σ::Float64, μ::Float64)
A simple normalized Gaussian with an integral of 1.
Takes as input wavelength array, and values of line center μ and line width σ.
"""
function gauss(λ::Array,  σ::Float64, μ::Float64)
    g = @. 1/σ/sqrt(2π) * exp(-(λ - μ)^2/2/σ^2)
    return g
end 



""" Quick and coarse estimate of the integrated flux within a wavelength window
without actual fitting.
"""
function estimate_line_snr(wave, flux; err=nothing)
    if isa(err, Nothing)
        median_fit = running_median(flux, 3, :sym, nan=:ignore)
        nuisance_be_gone = flux .- median_fit
        stddev = nanstd(nuisance_be_gone) .* ones(size(flux))
    else
        stddev = err
    end 
    uflux = flux .± stddev
    usum = nansum(uflux)
    snr = Measurements.value(usum)/Measurements.uncertainty(usum)
    if isnan(snr); snr = 0.1; end
    return snr
end 



"""     get_resolving_power_muse([, order=3])
## Optional arguments
- `order::Int`: The order of the polynomial to fit to the datapoints. Default is 3.
## Returns
- A callable `Polynomials.Polynomial` object which returns the linearly interpolated
  resolving power at the wavelength passed to it.
"""
function get_resolving_power_muse(;order=3::Int)
    ll = Array{Float64}(
        [4650.0, 5000.0, 5500.0, 6000.0, 6500.0, 7000.0,
         7500.0, 8000.0, 8500.0, 9000.0, 9350.0])
    rr = Array{Float64}(
        [1609, 1750, 1978, 2227, 2484, 2737,
         2975, 3183, 3350, 3465, 3506])
    drr = Array{Float64}([6, 4, 6, 6, 5, 4, 4, 4, 4, 5, 10])
    p = Polynomials.fit(ll, rr, order, weights=1. ./ drr)
    return p
end


""" load_neblines(infile)
## Input:
- `filepath::String`: The file to use. Must be a whitespace separated file,
  containing at least the columns 'name', 'lamvac', 'lamair', and 'foverha'.
"""
function load_neblines(infile="./static_data/neblines.dat")
    lines = CSV.File(infile, delim=" ", ignorerepeated=true, comment="#") |> DataFrame
    select!(lines, [:name, :lamvac, :lamair, :foverha])
    #@info "Loaded list of nebular lines"
    return lines
end 


"""    fill_in_fit_values(dict::Dict, fitresults::Dict, pix_coords::Tuple{Int,Int}; fitstats=NaN)
"""
function fill_in_fit_values!(
    dict::Dict{Any,Any}, fitresults::Dict, pix_coords::Tuple{Int,Int}; fitstats=NaN)
    row, column = pix_coords
    fitstats = fitresults[:fitstats].fitstat
    fitresults = fitresults[:results]
    @debug dict[:refline] fitresults[Symbol(dict[:refline])] 
    dict[:fitstats][row, column] = fitstats
    dict[:redshift][row, column, 1] = fitresults[Symbol(dict[:refline])].redshift.val
    dict[:redshift][row, column, 2] = fitresults[Symbol(dict[:refline])].redshift.unc
    dict[:fwhm][row, column, 1] = fitresults[Symbol(dict[:refline])].fwhm_kms.val
    dict[:fwhm][row, column, 2] = fitresults[Symbol(dict[:refline])].fwhm_kms.unc
    for l in keys(dict)
        if !(l in keys(fitresults))
            continue
        end
        dict[l][row, column, 1] = fitresults[l][:norm].val
        dict[l][row, column, 2] = fitresults[l][:norm].unc
        dict[l][row, column, 3] = fitresults[l][:norm].val / fitresults[l][:norm].unc
    end
end



""" `write_to_fits(filepath, mapsdict)`
"""
function write_to_fits(filepath::String, mapsdict::Dict)
    f = FITS(filepath, "w")
    keylist = sort(collect(keys(mapsdict)))
    primhead = mapsdict[:header]; primhead["CONTINUE"] = ""  #push!(primhead, "CONTINUE"="")
    sechead = mapsdict[:header]; sechead["CONTINUE"] = ""
    write(f, Float64[], header=primhead)
    # for k in keys(mapsdict)
    for k in keylist
        if k in [:header, :primhead, :refline]
            continue
        end 
        write(f, mapsdict[k], header=sechead, name=String(k))
        # write(f, mapsdict[k],  name=String(k))
    end 
    close(f)
    return nothing
end 


end  # End module
