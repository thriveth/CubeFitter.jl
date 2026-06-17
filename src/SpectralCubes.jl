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
Load data from a JWST/NIRSpec IFU datacube into a Julia `struct`.
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


"""    MUSECube(filepath; linelist_path) -> MUSECube
_Not yet implemented_
"""
mutable struct MUSECube <: AbstractSpectralCube
    setting::String
    resol_fit_degree::Int
    header::FITSIO.FITSHeader
    primheader::FITSIO.FITSHeader
    linelist::DataFrame
    lsf_fitter
    function MUSECube(filepath::String;
                      linelist_path::String=joinpath(datapath, "neblines.dat"),
                      lsf_polynomium_degree::Int=3,
                      )

        ddict = load_fits_cube(filepath, "r")
        linelist = load_neblines(linelist_path)
        ddict = convert_ergscms_Å_units(ddict, "MUSE")
        wave = ustrip.(ddict[:Wave])
        fluxcube = ustrip.(ddict[:Data])
        errscube = ustrip.(ddict[:Errs])
        header = ddict[:Header]
        primheader = ddict[:Primheader]
        itp = get_resolving_power("MUSE", setting=lsf_polynomium_degree)
    end
end


function Base.show(io::IO, cube::AbstractSpectralCube)
    cubesize = cube.fluxcube |> size
    print(
          io, 
          "\nIFU datacube w. metadata. \nInstrument: $(cube.instrument) \n" * 
          "Grating: $(uppercase(cube.grating))\n" *
          "Data dimensions: x: $(cubesize[1]), y: $(cubesize[2]), λ: $(cubesize[3])\n" *
          "Initial redshift: $(cube.z_init)\nReference transition: :$(cube.ref_line)"
         )
end
