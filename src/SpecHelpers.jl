module SpecHelpers

export fwhm_to_sigma, sigma_to_fwhm, wl_to_v, v_to_wl, deltawl_to_v, v_to_deltawl
export vactoair, airtovac

using Unitful, PhysicalConstants.CODATA2018

const ckms = SpeedOfLightInVacuum |> u"km/s"
const caps = SpeedOfLightInVacuum |> u"Å/s"

###==========================================
#           Functions

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
