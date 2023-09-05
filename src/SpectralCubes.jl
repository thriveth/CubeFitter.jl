using FITSIO
using DataFrames

# const datapath = joinpath(dirname(pathof(CubeFitter)), "..", "static_data")
datapath = "../static_data/"

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
        lsf_file_path="./static_data/jwst_nirspec_$(grating)_disp.fits",
        # lsf_file_path=joinpath(datapath, "jwst_nirspec_$(grating)_disp.fits"),
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
