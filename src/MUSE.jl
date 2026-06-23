# INSTRUMENT KIT. 
# ===============
# For each instrument, this file should contain: 
# - A an AbstractSpectralCube -derived mutable struct.
# - A function that returns the resolving power R for a given observed wavelength.
# - OPTIONALLY A Base.show() method taylored to said instrument, if the default 
#   Base.show(AbstractSpectralCube) doesn't get the job done.

"""    MUSECube(filepath; kwargs...) -> MUSECube
A MUSE datacube with metadata. The MUSE pipeline stores data in ergs/cm2/s/Å units, so 
`CubeFitter.MUSECube` assumes that the input data are in those units. 

# Arguments
- `filepath::String`: Path to data file in FITS datacube format.
# Keywords
- `linelist_path::String`: Path to user-provided line list. Defaults to using the built-in 
  line list. 
- `lsf_polynomium_degree::Int=3`: Degree of the polynomial used to fit/interpolate the MUSE
  Resolving power curve. 
- `z_init::Float`: An initial guess of the cosmological redshift of the observed object. 
  This should be within ~5% of the true value, otherwise the code doesn't know what to do 
  with the data. 
- `reference_line::Symbol`: The (strong) emission line to use for initial redshift and flux
  estimates. Defaults to :OIII_5007. NB! This parameter is important to set if the default 
  is not appropriate!
- `data_ext::String`: Name of the FITS extension with the flux data cube. 
  Defaults to "DATA". Can also be an extension number.
- `error_ext::String`: Name of the FITS extension with the variance data cube.
  Defaults to "STAT". Can also be an extension number. This is converted to standard 
  errors, since this is what `CubeFitter` uses internally. NB! This also means that the 
  output errors from `CubeFitter` are std errors, _not_ variances. 
"""
mutable struct MUSECube <: AbstractSpectralCube
    const instrument::String
    # setting::String
    # resol_fit_degree::Int
    wave::Array
    fluxcube::Array
    errscube::Array
    header::FITSHeader
    primheader::FITSHeader
    linelist::DataFrame
    lsf_fitter
    z_init
    ref_line
    flux_convert
    function MUSECube(filepath::String;
                      linelist_path::String=joinpath(datapath, "neblines.dat"),
                      lsf_polynomium_degree::Int=3, lsf_file_path=nothing,
                      z_init=0.0, reference_line=:OIII_5007, data_ext="DATA", 
                      error_ext="STAT", 
                      )
        linelist = load_neblines(linelist_path)
        ddict = load_fits_cube(filepath, data_ext=data_ext, error_ext=error_ext)
        origflux = deepcopy(ddict[:Data])
        wave = ustrip.(ddict[:Wave])
        fluxcube = ustrip.(ddict[:Data])
        errscube = sqrt.(ustrip.(ddict[:Errs]))
        header = ddict[:Header]
        primheader = ddict[:Primheader]
        itp = get_resolving_power("MUSE", setting=lsf_polynomium_degree)
        fluxconverter = nanmean(origflux ./ fluxcube, dims=(1,2))
        println("MUSE cube successfully loaded")
        new("MUSE", wave, fluxcube, errscube, header, primheader, linelist, itp, z_init, reference_line, fluxconverter)
    end
end

"""     get_resolving_power_muse(; order=3)
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


# # How to show the type in the REPL 
# # (equivalent to the .__repr__() function in Python)
# function Base.show(io::IO, cube::MUSECube)
#     cubesize = cube.fluxcube |> size
#     outstr = "IFU datacube w. metadata. \nInstrument: $(cube.instrument) \n"
#     outstr *= "Data dimensions: x: $(cubesize[1]), y: $(cubesize[2]), λ: $(cubesize[3])\n"
#     outstr *= "Initial redshift: $(cube.z_init)\nReference transition: :$(cube.ref_line)"
#     print( io, outstr)
# end
