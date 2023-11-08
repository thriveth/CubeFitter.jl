#Threads, because it realy can speed things up
using Base.Threads
using FITSIO
using Printf
using FastRunningMedian
# Because normal statistics are not NaN-friendly
using NaNStatistics
# For sigma clipping:
using Photometry
using Interpolations

#= include("./SpecHelpers.jl") =#
### This part is only relevant for my Sunburst NIRSpec project
# Set up the paths and filenames necessary to do the thing.
in_path = "../Products/NIRSpec/science-ready/"
out_path = "../Products/ContSub/"

pathsdict = Dict(
    "pos1_blu" => in_path * "sunburst-P1-sigmaclipped-s3d-g140h.fits", 
    "pos1_red" => in_path * "sunburst-P1-sigmaclipped-s3d-g235h.fits", 
    "pos2_blu" => in_path * "sunburst-P2-sigmaclipped-s3d-g140h.fits", 
    "pos2_red" => in_path * "sunburst-P2-sigmaclipped-s3d-g235h.fits", 
    "pos3_blu" => in_path * "sunburst-P3-sigmaclipped-s3d-g140h.fits", 
    "pos3_red" => in_path * "sunburst-P3-sigmaclipped-s3d-g235h.fits")

outpaths = Dict(
    "pos1_blu" => out_path * "P1-contsub-s3d-g140h.fits", 
    "pos1_red" => out_path * "P1-contsub-s3d-g235h.fits", 
    "pos2_blu" => out_path * "P2-contsub-s3d-g140h.fits", 
    "pos2_red" => out_path * "P2-contsub-s3d-g235h.fits", 
    "pos3_blu" => out_path * "P3-contsub-s3d-g140h.fits", 
    "pos3_red" => out_path * "P3-contsub-s3d-g235h.fits")
### End of Sunburst-NIRSpec-only part


"""    load_fits_cube_contsub()
Convenience function for this context only. Loads an IFU datacube from a FITS file,
returns a dict of the essentials.
## Input:
- `inpath::String`: path to the file to open.
## Returns:
- `out::Dict`: A dictionary of datacube, errorcube, wavelength vector,
  FITS primary and first extension header.
"""
function load_fits_cube_contsub(inpath::String)::Dict
    # TODO: Phase this out! 
    fitsfile = FITS(inpath, "r")
    header = read_header(fitsfile[2])
    primary = read_header(fitsfile[1])
    data = read(fitsfile[2])
    errs = read(fitsfile[3])

    naxis1 = header["NAXIS1"]
    naxis2 = header["NAXIS2"]
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

    out = Dict(
        "Wave" => waves,
        "Data" => data,
        "Errs" => errs,
        "Header" => header,
        "Primary header" => primary)

    return out
end


"""
Leaner version of the cont_subt function. This functions does nothing else than
continuum subtract the data cube and return the result. 
## Input
- `incube::AbstractSpectralCube`: The 3D array containing the non-subtracted datacube. We are assuming that
  it is... BLA BLA
"""
function cont_subt(incube; 
        #= TODO: Finish rewriting this thingie =#
        clip_sigmas=3., window_pix=51, 
        neblines=joinpath(datapath, "neblines.dat"),
        #= neblines="./static_data/neblines.dat", =#
        velwidth=1000,  # km/s
        redshift=0,
        mask_lines=false,
        min_strength=0.01,  # Fraction of Hα
    )::Array{Float64, 3}
    inarray = incube[:Data]
    nlines = load_neblines(neblines)
    contsub_cube = zeros(eltype(inarray), size(inarray))
    if mask_lines
        mask = falses(size(inarray)[3])
        for r in eachrow(nlines)
            if r[:foverha] < min_strength
                continue
            end
            cenwave = mod[comp].lab_wave.val * (1 + mod[comp].redshift.val)
            window_width_ang = v_to_deltawl(window_width_kms, cenwave)
            @debug "Make line mask: " wave window_width_ang
            mask[cenwave-window_width_ang .< wave .< cenwave+window_width_ang] .= 1
        end 
    end
    # println(size(contsub_cube))
    naxis1, naxis2 = size(inarray)
    # Iterate over the spatial dimensions and populate the spectral datacube
    @threads for i in 1:naxis1
        for j in 1:naxis2
            @printf("Now continuum subtracting row %03d, col. %03d. \r", i, j)
            flush(stdout)
            spec = inarray[i, j, :] 
            scspec = copy(spec)
            scspec[abs.(spec) .> nanmedian(spec) + clip_sigmas .* nanstd(spec)] .= NaN
            contsub_cube[i, j, :] = spec - running_median(scspec, 51, nan=:ignore)
        end
    end 
    return contsub_cube
end


""" Continuum subtract a FITS Datacube file, and optionally write the result to a new file.

This function contains the plumbing, or bookkeeping, or whatever you call it. For the actual
subraction, it calls the function `cont_subt()`. This function is in turn called by
`contsub_many()`.

This function requires the datadicts to be populated with reasonable names and file paths.
## Input
- `inpath::String`: Path to the non-subtracted datacube file.
## Optional input
- `saveit::Bool`: Whether or not to save the generated output. For safety reasons, it is
  set to `false` by default.
- `outpath::String`: Path to the output file.
  Default: `"continuum_subtracted.fits"`.
- `clip_sigmas::Float64`: Continuum subtraction works on sigma clipped data. This argument
  allows to customize how strong the clipping should be. Default value is `3`.
- `window_pix::Int`: Pixel width of the rolling window in which the median is taken.
  Default value is `51`. Odd numbers give better results.
## Returns
- `out::Dict`: Mainly for convenience, the function returns a `Dict` containing a continuum
  subtracted data cube, an unaltered error cube, as well as FITS headers from both the primary
  and the data extension of the input FITS file.
"""
function continuum_subtract_and_save(
    inpath::String;
    saveit=false::Bool,
    outpath="continuum_subtracted.fits"::String,
    clip_sigmas=3.::Float64,
    window_pix=51::Int,
    )::Dict
    # This is just for convenience, if I want to define some shortcuts for some of the file names,
    # making it easier to call the function interactively.
    if inpath in keys(pathsdict)
        datadict = load_fits_cube_contsub(pathsdict[inpath])
    else 
        datadict = load_fits_cube_contsub(inpath)
    end 

    contsubt = cont_subt(datadict["Data"], clip_sigmas=clip_sigmas, window_pix=window_pix)

    if saveit == true
        outfile = FITS(outpath, "w")
        write(outfile, Float64[], header=datadict["Primary header"])
        write(outfile, contsubt; header=datadict["Header"], name="SCI")
        write(outfile, datadict["Errs"], header=datadict["Header"], name="ERR")
        close(outfile)
    end
    datadict["Cont.subt."] = contsubt
    return datadict
end


"""
# Continuum subtract one or more files
Continuum subtracts a a single file and save the result to a new file.
"""
function contsub_many(which="all"::String; saveit=false::Bool)
    for (thing, inpath) in pathsdict
        if which != "all"
            if which != thing
                continue
            end
        end
        println(" Input file:  $(inpath)\n Output file:  $(outpaths[thing])")
        data, wave, contsub = cont_subt(inpath, saveit=saveit, outpath=outpaths[thing], setting=thing)
    end
    return nothing
end
