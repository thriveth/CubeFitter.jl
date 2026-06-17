datapath = joinpath(dirname(pathof(CubeFitter)), "..", "static_data")


# Declaring an abstract SpectralCube type will allow us to make SpectralCube objects
# inheriting from it for various ls instruments, which require different treatments and
# different kind of, and still be able to write functions which can operate on all of them. 
abstract type AbstractSpectralCube end


"""    GenericSpectralCube(filepath::String; kwargs...)
Thus far, this is not actually implemented.
It is planned to be a barebones solution for cubes that are not from any of the
explicitly supported.
"""
mutable struct GenericSpectralCube <: AbstractSpectralCube
    instrument::String
    wave::Array
    fluxcube::Array
    errscube::Array
end


"""    NIRSpecCube(filepath, grating; kwargs...) -> NIRSpecCube
A NIRSpec IFU datacube with errors and metadata. Assumes that the input is in the format 
saved by the NIRSpec data pipeline.

# Arguments
- `filepath::String`: Path to data file in FITS datacube format.
- `grating::String`: The grating used to make the observations. 
# Keywords
- `linelist_path::String`: Path to user-provided line list. Defaults to using the built-in 
  line list. 
- `lsf_file_path::String`: Path to any non-default version of the LSF solution file. This 
  must at present time have the same format as the standard dispersion files included in 
  this package.
- `z_init::Float`: An initial guess of the cosmological redshift of the observed object. 
  This should be within ~5% of the true value, otherwise the code doesn't know what to do 
  with the data. 
- `reference_line::Symbol`: The (strong) emission line to use for initial redshift and flux
  estimates. Defaults to :OIII_5007. NB! This parameter is important to set if the default 
  is not appropriate!
- `data_ext::String`: Name of the FITS extension with the flux data cube. 
  Defaults to "SCI". Can also be an extension number.
- `error_ext::String`: Name of the FITS extension with the standard error data cube.
  Defaults to "ERR". Can also be an extension number.
"""
mutable struct NIRSpecCube <: AbstractSpectralCube
    const instrument::String
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
    flux_convert
    function NIRSpecCube(filepath::String, grating::String; 
            linelist_path=joinpath(datapath, "neblines.dat"), 
            lsf_file_path=joinpath(datapath, "jwst_nirspec_$(grating)_disp.fits"), 
            z_init=0.0, reference_line=:OIII_5007, data_ext="SCI", error_ext="ERR") 
        grating = lowercase(grating)
        linelist = load_neblines(linelist_path)
        ddict = load_fits_cube(filepath, data_ext=data_ext, error_ext=error_ext)
        origflux = deepcopy(ddict[:Data])
        ddict = convert_ergscms_Å_units(ddict; instrument="NIRSpec")
        wave = ustrip.(ddict[:Wave])
        fluxcube = ustrip.(ddict[:Data])
        errscube = ustrip.(ddict[:Errs])
        header = ddict[:Header]
        primheader = ddict[:Primheader]
        itp = get_resolving_power("NIRSpec", setting=grating)
        fluxconverter = nanmean(origflux ./ fluxcube, dims=(1,2))
        println("NIRSpecCube successfully loaded.")
        new("NIRSpec", grating, wave, fluxcube, errscube, header, primheader, linelist, 
            itp, z_init, reference_line, fluxconverter)
    end
end


"""    MIRICube(filepath, grating; kwargs...) -> MIRICube
_Not yet implemented_
"""
mutable struct MIRICube <: AbstractSpectralCube
    const instrument::String
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
    function MIRICube(filepath::String, grating::String;
        linelist_path=joinpath(datapath, "neblines.dat"),
        lsf_file_path=joinpath(datapath, "jwst_miri_$(grating)_disp.fits"),
        z_init=0, reference_line=:OIII_5007)
        linelist = load_neblines(linelist_path)
        ddict = load_fits_cube(filepath)
        ddict = convert_ergscms_Å_units(ddict, instrument="NIRSpec")
        wave = ustrip.(ddict[:Wave])
        fluxcube = ustrip.(ddict[:Data])
        errscube = ustrip.(ddict[:Errs])
        header = ddict[:Header]
        primheader = ddict[:Primheader]
        itp = get_resolving_power("NIRSpec", setting=grating)
        new("MIRI", grating, wave, fluxcube, errscube, header, primheader, 
            linelist, itp, z_init, reference_line)
    end
end


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



function Base.show(io::IO, cube::AbstractSpectralCube)
    cubesize = cube.fluxcube |> size
    outstr = "\nIFU datacube w. metadata. \nInstrument: $(cube.instrument) \n"
    if cube.instrument in ["NIRSpec", "MIRI"]
        outstr *= "Grating: $(uppercase(cube.grating))\n"
    end
    outstr *= "Data dimensions: x: $(cubesize[1]), y: $(cubesize[2]), λ: $(cubesize[3])\n"
    outstr *= "Initial redshift: $(cube.z_init)\nReference transition: :$(cube.ref_line)"

    print(
          io, 
          outstr
          # "\nIFU datacube w. metadata. \nInstrument: $(cube.instrument) \n" * 
          # "Grating: $(uppercase(cube.grating))\n" *
          # "Data dimensions: x: $(cubesize[1]), y: $(cubesize[2]), λ: $(cubesize[3])\n" *
          # "Initial redshift: $(cube.z_init)\nReference transition: :$(cube.ref_line)"
         )
end
