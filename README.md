# CubeFitter

<!-- [![Build Status](https://github.com/thriveth/CubeFitter.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/thriveth/CubeFitter.jl/actions/workflows/CI.yml?query=branch%3Amain) -->

## Introduction

A package for automatically fitting emission lines in astronomical spectral cubes.
Hopefully, this could develop into a more general package for handling spectral cubes in
Julia.

Development is still in early days, so everything may break - but currently, the NIRSpec
related code _Works For Me™_. I am in the process of testing it for MIRI and MUSE, and
most other instruments should be quite easy to add as well.


## What does it do (and how)?

This code contains three main components (which may get separated into
their own packages later): 

- A family of `SpectralCube` types which handle the data and metadata.
- A continuum subtraction function (see description below)
- The functionality for fitting emission lines.

The workhorses of this package are the various types `NIRSpecCube`, `MUSECube`,
etc.; as well as the functions `cont_subt()`, `fit_cube()`, and
`fit_spectrum_from_subcube()`.


## Usage

### Enabling the package

This package is not yet registered to install (maybe later). To run it now, clone the
repository to your preferred location; then:

```julia
julia> import Pkg
julia> Pkg.activate("/path/to/CubeFitter.jl")
julia> Pkg.instantiate()  # Install dependencies, needs only be done once.
julia> using CubeFitter
```

Alternatively, you can enter the `Pkg>` prompt, run `activate /path/to/CubeFitter.jl`,
then (first time) `instantiate()`. Press Backspace to return to the normal `julia>`
prompt, and run `using CubeFitter`.

### Quick example

A simple example of a session using `CubeFitter.jl`:

```julia
julia> cube = NIRSpecCube("/path/to/datacube.fits", "g140m"; z_init=2.43)
julia> cscube = cont_subt(cube)
julia> write_spectral_cube_to_fits("mycube.fits", cscube)
julia> out = fit_cube(cscube)
julia> write_maps_to_fits("/desired/path/to/output/file.fits", out)
```

It is of course possible to just load the continuum subtracted cube and not go
through that step again. 

The function `fit_cube()` will return a dictionary of $N\times M \times 3$ arrays, one for
each fitted line, as well as one for the 0th to 2nd moment, the spaxel-wise redshift, and
line width. Each such array will have the measured values as the first slice, the
parameter error in the second slice, and the fit statistics in the last. Each of these
slices will be saved as an individual HDU in the output FITS file.

The function fit_cube is mainly a convenience function enabling one to cycle through every
spaxel in the cube, and running the fits with default settings. One can also use the
function `fit_spectrum_from_subcube()` directly. This function allows to extract a
spectrum from the cube based on an x- and y range passed to it (this must always be a
range; if only one spaxel is wanted, give the range for the spaxel (_i_, _j_) as
`xrange=i:i, yrange=j:j`). The function extracts a simple, un-weighted spectrum over the
given (_x, y_) range, runs the fit, and return a dictionary of the fit results and
statistics, along with the extracted spectrum, errors, and wavelength range for
convenience. See the function docstring to learn more.


### Data format

The function `fit_cube()` expects the data to already be continuum subtracted, but
otherwise saved in the same data format as the final pipeline products of a given
instrument; this can be done using the `cont_subt()` function as described above. Using
JWST/NIRSpec as an example, this means that the data cubes should be saved as a FITS file
with an empty primary HDU, and the flux and error cubes saved as the first and second
extensions, respectively. At the moment, this requirement is hardcoded, and any other
formats will either throw an error or, at worst, yield wrong and meaningless results.


## Screens

There really isn't much to look at but here it is:

![Screenshot of CubeFitter in action](./Screenshots/CubeFitter.png)


## Planned features / whishlist

In order of approximate priority: 

- [x] Include continuum subtraction functionality in the package.
- [x] Make it possible to add a second and perhaps third kinematic component. 
- [ ] Test and ensure the `MUSECube` and `MIRICube` structs actually works as advertised
- [ ] Allow for measuring flux or upper limits in lines that are currently excluded as
      having too low S/N ratio.
- [ ] Add support for more instruments. Suggestions welcome (especially if accompanied
      with a suitable test dataset).
- [ ] Write quicklook-functions allowing to quickly view the fit outputs with minimum
      input. But still be tinker-friendly, don't hide stuff from the user.
