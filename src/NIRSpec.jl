# INSTRUMENT KIT. 
# ===============
# For each instrument, this file should contain: 
# - A an AbstractSpectralCube -derived mutable struct.
# - A function that returns the resolving power R for a given observed wavelength.
# - A Base.show() method taylored to said instrument.

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
        # itp = get_resolving_power("NIRSpec", setting=grating)
        itp = get_resolving_power_nirspec(grating;)
        fluxconverter = nanmean(origflux ./ fluxcube, dims=(1,2))
        # println("NIRSpecCube successfully loaded.")
        @info "NIRSpecCube successfully loaded."
        new("NIRSpec", grating, wave, fluxcube, errscube, header, primheader, linelist, 
            itp, z_init, reference_line, fluxconverter)
    end
end



"""     get_resolving_power_nirspec(grating[, calib_path="../static_data/"])
This function assumes that the FITS files containing the NIRSpec grating resolving power
are named "jwst_nirspec_[grating]_disp.fits", as fetched from the official JDox.
# Arguments
- `grating::String`: The name of the grating used.
# Keywords
- `calib_path::String`: Path to the location of the calibration data, defaults to  the data folder installed with the package.
# Returns
- A callable `Interpolations.Extrapolation` objects which returns the linearly 
  interpolated resolving power at the wavelength passed to it.
"""
function get_resolving_power_nirspec(grating::String; calib_path=datapath) 
    lsffilename = "jwst_nirspec_$(grating)_disp.fits"
    @debug "Path to LSF file: " calib_path * lsffilename
    fitsfile = FITS(joinpath(calib_path, lsffilename))
    angwave = read(fitsfile[2], "WAVELENGTH") * u"μm" .|> u"angstrom"
    R = read(fitsfile[2], "R")
    @debug extrema(angwave) extrema(R)
    itp = LinearInterpolation(ustrip.(angwave), R, extrapolation_bc=Flat())
    return itp
end


# How to show the type in the REPL 
# (equivalent to the .__repr__() function in Python)
function Base.show(io::IO, cube::NIRSpecCube)
    cubesize = cube.fluxcube |> size
    outstr = "IFU datacube w. metadata. \nInstrument: $(cube.instrument) \n"
    outstr *= "Grating: $(uppercase(cube.grating))\n"
    outstr *= "Data dimensions: x: $(cubesize[1]), y: $(cubesize[2]), λ: $(cubesize[3])\n"
    outstr *= "Initial redshift: $(cube.z_init)\nReference transition: :$(cube.ref_line)"
    print(
          io, 
          outstr
         )
end


