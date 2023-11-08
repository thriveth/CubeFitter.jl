module SpecHelpers

export fwhm_to_sigma, sigma_to_fwhm, wl_to_v, v_to_wl, deltawl_to_v, v_to_deltawl
export vactoair, airtovac, load_neblines, load_fits_cube

using Unitful, PhysicalConstants.CODATA2018, CSV, DataFrames, FITSIO

const ckms = SpeedOfLightInVacuum |> u"km/s"
const caps = SpeedOfLightInVacuum |> u"Å/s"

#= staticdatapath = joinpath(dirname(pathof(SpecHelpers)), "..", "static_data") =#


###=========================================================
#           Convenience functions used here and there.
"""    
    load_neblines(infile)
Loads a table of spectral lines, good for e.g. masking.
# Arguments
- `filepath::String`: The file to use. Must be a whitespace separated file,
  containing at least the columns 'name', 'lamvac', 'lamair', and 'foverha'.
"""
function load_neblines(infile)
    lines = CSV.File(infile, delim=" ", ignorerepeated=true, comment="#") |> DataFrame
    select!(lines, [:name, :lamvac, :lamair, :foverha])
    #@info "Loaded list of nebular lines"
    return lines
end 


"""    load_fits_cube(inpath::String)
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



###==========================================
#           Often-used computations

function fwhm_to_sigma(fwhm)
    σ = fwhm / (2. * sqrt(2. * log(2.)))
    return σ
end


function sigma_to_fwhm(σ)
    fwhm = σ * 2. * sqrt(2. * log(2.))
    return fwhm
end 


""" Converts a wavelength range and a reference wavelength to a velocity range.
Uses the relation:
>     λ / λ0 = 1 + v / c =>
>          v = (1 + λ / λ0) * c
Velocities are in km/s.
"""
function wl_to_v(wave::Float64, wl0::Float64)::Float64
    v = (wave / wl0 - 1.) * ckms / 1000.
    return v 
end 


""" Converts a velocity range and a reference wavelength to a wavelength range.
Uses the relation:
>       λ / λ0 = 1 + v / c =>
>            λ = λ0 * (1 + v / c)
Velocities are in km/s. """
function v_to_wl(v, wl0)
    wave = wl0 * (1. + v / ustrip(ckms))  #(con.c / 1000.))
    return wave
end 


""" Gives Δλ as function of v and λ_0.
>      v = c Δλ / λ0 =>
>     Δλ = v λ0 / c
## Input:
- `v::Float64`: Velocity in km/s 
- `wl0::Float64`: Reference wavelegth in Å.
"""
function v_to_deltawl(v, λ0)
    Δλ = v * λ0 / ustrip(ckms)
    return Δλ
end 


""" Gives delta-lambda as function of v and lambda-0.
>      v = c Δλ / λ0 =>
>     Δλ = v λ0 / c
## Input:
- `Δλ::Float64`: Wavelength difference in Å
- `λ0::Float64`: Reference wavelegth in Å.
"""
function deltawl_to_v(Δλ, λ0)
    v = ustrip(ckms) * Δλ / λ0
    return v
end 



"""    airtovac(airwl[, nouvconv=true])
Translated from the Astropysics
implementation: http://packages.python.org/Astropysics/ To not have
this module as dependency as that is the only one necessary.

Returns vacuum wavelength of the provided air wavelength array or scalar.
Good to ~ .0005 angstroms.

If nouvconv is True, does nothing for air wavelength < 2000 angstroms.

Input must be in angstroms.

Adapted from idlutils airtovac.pro, based on the IAU standard
for conversion in Morton (1991 Ap.J. Suppl. 77, 119)
"""
function airtovac(airwl::AbstractFloat, nouvconv=true::Bool)
    if airwl < 2000 && nouvconv
        return airwl
    else 
        sig2 = (1e4 / airwl)^2
        convfact = (1. + 6.4328e-5 + 2.94981e-2 / (146. - sig2)
            + 2.5540e-4 / (41. - sig2))
        return airwl * convfact
    end
end


""" Also copied from Astropysics, see air_to_vacuum() above. Returns
air wavelength of the provided vacuum wavelength array or scalar. Good
to ~ .0005 angstroms.

If nouvconv is True, does nothing for air wavelength < 2000 angstroms.

Input must be in angstroms.

Adapted from idlutils vactoair.pro.
"""
function vactoair(vacwl::AbstractFloat, nouvconv=true::Bool)
    vacwl < 2000 && nouvconv && return vacwl
    wave2 = vacwl^2
    convfact = 1.0 + 2.735182e-4 + 131.4182/wave2 + 2.76249e8/(wave2*wave2)
    return vacwl/convfact
end

end  # End module
