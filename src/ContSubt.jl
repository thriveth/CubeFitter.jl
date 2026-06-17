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
# using PyPlot

"""    cont_subt(incube::AbstractSpectralCube; ..)

Leaner version of the cont_subt function. This functions does nothing else than
continuum subtract the data cube and return the result. 
## Input
- `incube::AbstractSpectralCube`: The datacube (type `CubeFitter.AbstractSpectralCube`) containing the 
  data to subtract.
- `window_pix::Number`: Number of pixels (preferably odd number) to indluce in the running median kernel.
- `velwidth::Number`: Velocity width to assume for each emission line, given in km/s.
- `mask_lines::Bool`: Whether to mask out known emission lines before computing the running median.
- `min_strength::Number`: The minimum strength of an emission line (relative to Hα in the line list) 
  to be masked out before computing running median. This is ignored if `mask_lines` is set to `false`.
- `fancy::Bool`: Whether or not to use experimental "fancy" continuum modeling method. 
  At the moment, the "fancy" version is worse than the "simple" one and also slower, so don't use it.
"""
function cont_subt(incube; 
        # TODO: Rewrite to use cont_subt_spectrum for internals 
        # to avoid duplication of code
        window_pix=51, 
        sigma_clip=true, clip_sigmas=3.,
        velwidth=350,  # km/s
        mask_lines=false,
        min_strength=0.01,  # Fraction of Hα
        fancy=false,
    )
    outcube = deepcopy(incube)
    inarray = outcube.fluxcube
    wave = outcube.wave
    redshift = outcube.z_init
    nlines = outcube.linelist
    contsub_cube = zeros(eltype(inarray), size(inarray))
    mask = trues(size(inarray)[3])
    if mask_lines
        for r in eachrow(nlines)
            if r[:foverha] < min_strength
                continue
            end
            #= cenwave =  mod[comp].lab_wave.val * (1 + mod[comp].redshift.val) =#
            cenwave = r[:lamvac] * (1 + redshift)
            window_width_ang = v_to_deltawl(velwidth, cenwave)
            @debug "Make line mask: " wave window_width_ang
            mask[cenwave-window_width_ang .< wave .< cenwave+window_width_ang] .= 0
            # So we can interpolate properly:
            mask[1] = mask[end] = 1
        end 
    end
    naxis1, naxis2 = size(inarray)
    # Iterate over the spatial dimensions and populate the spectral datacube
    @track for i in 1:naxis1
        for j in 1:naxis2
            @printf("Now continuum subtracting row %03d, col. %03d. \r", i, j)
            flush(stdout)
            if fancy
                # Despite being fancy, this method is actually worse!
                # But I keep it because it is good to have infrastructure that
                # calls an external function, then can just rewrite the 
                # external function.
                specvec = inarray[i, j, :] |> deepcopy
                inspec = DataFrame(Dict(:wave => wave, :flux => specvec))
                mask = .~mask
                mask .|= (specvec .|> isnan)
                mask[1] = mask[end] = 0
                outspec = subtract_continuum(inspec, mask; pix_window=window_pix,  sig_lim=clip_sigmas)
                ipl = LinearInterpolation(wave,outspec[!, :continuum])(wave)
            else
                # Despite its simplicity, this method seems unreasonably good!
                specvec = inarray[i, j, :] |> deepcopy
                spec = inarray[i, j, :][mask]
                #= Sigma clip the (possibly line-masked) spectrum =#
                scwave = wave[mask]
                scspec = deepcopy(spec)
                scmask = trues(size(scwave))
                scmask[abs.(scspec) .> nanmedian(scspec) + clip_sigmas .* nanstd(scspec)] .= 0
                scmask[1] = scmask[end] = 1  # Must be included to allow interpolationc
                finalspec, finalwave = scspec[scmask], scwave[scmask]
                continuum = running_median(finalspec, window_pix, nan=:ignore)
                ipl = LinearInterpolation(finalwave, continuum)(wave)
            end
            contsub_cube[i, j, :] = inarray[i, j, :] - ipl#(wave)
            # contsub_cube[i, j, :] = ipl#(wave)  # DEBUG
        end
    end 
    outcube.fluxcube = contsub_cube
    return outcube
end

# function make_lines_mask(linelist, min_strength=0.01)
# end
#
#

"""    cont_subt_spectrum(input_spectrum; kwargs...)

`kwargs` are the same as input to `cont_subt`
Expects the input spectrum to be a DataFram with columns :wave, :flux/:fnu/:flam, and :dflux/:dfnu/:dflam
"""
function cont_subt_spectrum(
    inspec; 
    linelist=nothing, window_pix=41, sigma_clip=true, clip_sigmas=3., velwidth=350,  # km/s
    redshift=0.0, mask_lines=false, fancy=false, min_strength=0.02,  # Fraction of Hα 
    )
    if "flux" in inspec |> names 
        flab, dflab, clab = :flux, :dflux, :cont
    elseif "flam" in inspec |> names
        flab, dflab, clab = :flam, :dflam, :contflam
    else
        flab, dflab, clab = :fnu, :dfnu, :contfnu
    end
    # inspec.names = Symbol.(inspec.names)
    outspec = inspec |> deepcopy
    # rename!(outspec, Symbol.(inspec |> names))
    rename!(string, outspec)
    # return outspec
    wave = inspec[!,:wave]
    mask = trues(inspec[!, 1] |> length)
    if mask_lines
        for r in eachrow(linelist)
            # if r[!,:foverha] < min_strength
            if r.foverha < min_strength
                continue
            end
            #= cenwave =  mod[comp].lab_wave.val * (1 + mod[comp].redshift.val) =#
            cenwave = r.lamvac * (1 + redshift)
            window_width_ang = v_to_deltawl(velwidth, cenwave)
            @debug "Make line mask: " wave window_width_ang
            mask[cenwave-window_width_ang .< wave .< cenwave+window_width_ang] .= 0
            # So we can interpolate properly:
            mask[1] = mask[end] = 1
        end 
    end

    if fancy
        # Despite being fancy, this method is actually worse!
        # But I keep it because it is good to have infrastructure that
        # calls an external function, then can just rewrite the 
        # external function.
        mask = .~mask
        mask .|= (inspec[!,:flux] .|> isnan)
        mask[1] = mask[end] = 0
        outspec = subtract_continuum(inspec, mask; pix_window=window_pix,  sig_lim=clip_sigmas)
        ipl = LinearInterpolation(wave,outspec[!, :continuum])(wave)
    else
        # Despite its simplicity, this method seems unreasonably good!
        # specvec = inarray[i, j, :] |> deepcopy
        # spec = inarray[i, j, :][mask]
        spec = inspec[!, flab][mask]
        #= Sigma clip the (possibly line-masked) spectrum =#
        scwave = wave[mask]
        scspec = deepcopy(spec)
        scmask = trues(size(scwave))
        scmask[abs.(scspec) .> nanmedian(scspec) + clip_sigmas .* nanstd(scspec)] .= 0
        scmask[1] = scmask[end] = 1  # Must be included to allow interpolationc
        finalspec, finalwave = scspec[scmask], scwave[scmask]
        continuum = running_median(finalspec, window_pix, nan=:ignore)
        ipl = LinearInterpolation(finalwave, continuum)(wave)
    end
    outspec[!, clab] = ipl
    outspec[!, flab] = inspec[!, flab] .- ipl
    return outspec
end


# function subtract_continuum(spec, mask=nothing; pix_window=31, sig_lim=3, plotit=false)
#     # STRATEGY: 
#     # 1. Mask around lines
#     # 2. Smooth the remaining pixels, subtract the smoothed.
#     # 3. Take Std. of diff, add to massk all where diff is more than
#     #    sig_lim std's. 
#     # if mask isa Nothing; mask = zeros(length(spec)); end
#     origmask = mask |> copy
#     specm = spec[.~mask, :]
#     med1 = running_median(specm[!, :flux], pix_window, nan=:ignore)
#     # Back to orig wave grid
#     # println(size(specm), "  ,  ", size(spec))
#     med1ipl = LinearInterpolation(specm[!, :wave], med1)(spec[!, :wave])
#     med2ipl = running_median(med1ipl, (pix_window * 3) + 0, nan=:ignore)
#     meddiff = med1ipl - med2ipl
#     std2 = nanstd(meddiff)
#     mask2 = abs.(meddiff) .> sig_lim * std2
#     mask .|= mask2
#     spec[!, :contmask] = mask
#
#     specm = spec[.~mask, :]
#     med3 = running_median(specm[!, :flux], pix_window, nan=:ignore)
#     med3ipl = LinearInterpolation(specm[!, :wave], med3)(spec[!, :wave])
#     ### NOW, 
#     #   - subtract continuum-thus-far from original data, 
#     #   - Mask out known mask thingies, 
#     #   - Sigma clip the rest, add to mask.
#     meddiff3 = spec[!, :flux] .- med3ipl
#     std3 = nanstd(meddiff3[.~mask])
#     mask3 = abs.(meddiff3) .> sig_lim * std3
#     mask .|= mask3
#     spec[!, :contmask] = mask
#     specm = spec[.~mask, :]
#     med4 = running_median(specm[!, :flux], pix_window, nan=:ignore)
#     med4ipl = LinearInterpolation(specm[!, :wave], med4)(spec[!, :wave])
#     spec[!, :continuum] = med4ipl
#     if plotit
#         figure()
#         plot(spec[!, :wave]./zz, spec[!, :flux])
#         plot(spec[!, :wave]./zz, med3ipl)
#         plot(spec[!, :wave]./zz, mask)
#         plot(spec[!, :wave]./zz, origmask .* 0.9)
#         plot(spec[!, :wave]./zz, med4ipl)
#         plot(spec[!, :wave]./zz, mask .* 1.1)
#     end
#     return spec
# end
#
#
# """    continuum_subtract_and_save()
#
# Continuum subtract a FITS Datacube file, and optionally write the result to a new file.
#
# This function contains the plumbing, or bookkeeping, or whatever you call it. For the actual
# subraction, it calls the function `cont_subt()`. This function is in turn called by
# `contsub_many()`.
#
# This function requires the datadicts to be populated with reasonable names and file paths.
# ## Input
# - `inpath::String`: Path to the non-subtracted datacube file.
# ## Optional input
# - `saveit::Bool`: Whether or not to save the generated output. For safety reasons, it is
#   set to `false` by default.
# - `outpath::String`: Path to the output file.
#   Default: `"continuum_subtracted.fits"`.
# - `clip_sigmas::Float64`: Continuum subtraction works on sigma clipped data. This argument
#   allows to customize how strong the clipping should be. Default value is `3`.
# - `window_pix::Int`: Pixel width of the rolling window in which the median is taken.
#   Default value is `51`. Odd numbers give better results.
# ## Returns
# - `out::Dict`: Mainly for convenience, the function returns a `Dict` containing a continuum
#   subtracted data cube, an unaltered error cube, as well as FITS headers from both the primary
#   and the data extension of the input FITS file.
# """
# function continuum_subtract_and_save(
#     inpath::String;
#     saveit=false::Bool,
#     outpath="continuum_subtracted.fits"::String,
#     clip_sigmas=3.::Float64,
#     window_pix=51::Int,
#     )::Dict
#     # This is just for convenience, if I want to define some shortcuts for some of the file names,
#     # making it easier to call the function interactively.
#     if inpath in keys(pathsdict)
#         datadict = load_fits_cube_contsub(pathsdict[inpath])
#     else 
#         datadict = load_fits_cube_contsub(inpath)
#     end 
#
#     contsubt = cont_subt(datadict["Data"], clip_sigmas=clip_sigmas, window_pix=window_pix)
#
#     if saveit == true
#         outfile = FITS(outpath, "w")
#         write(outfile, Float64[], header=datadict["Primary header"])
#         write(outfile, contsubt; header=datadict["Header"], name="SCI")
#         write(outfile, datadict["Errs"], header=datadict["Header"], name="ERR")
#         close(outfile)
#     end
#     datadict["Cont.subt."] = contsubt
#     return datadict
# end
#

# """
# # Continuum subtract one or more files
# Continuum subtracts a a single file and save the result to a new file.
# """
# function contsub_many(which="all"::String; saveit=false::Bool)
#     for (thing, inpath) in pathsdict
#         if which != "all"
#             if which != thing
#                 continue
#             end
#         end
#         println(" Input file:  $(inpath)\n Output file:  $(outpaths[thing])")
#         data, wave, contsub = cont_subt(inpath, saveit=saveit, outpath=outpaths[thing], setting=thing)
#     end
#     return nothing
# end
