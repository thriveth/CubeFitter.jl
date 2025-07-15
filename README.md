# CubeFitter

[![Build Status](https://github.com/thriveth/CubeFitter.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/thriveth/CubeFitter.jl/actions/workflows/CI.yml?query=branch%3Amain)

![](./Media/Logo.png)

## Introduction

A package for automatically fitting emission lines in astronomical spectral cubes.
Hopefully, this could develop into a more general package for handling spectral cubes in
Julia.

Development is still in early days, so everything may break - but currently, the NIRSpec
related code _Works For Me™_. I am in the process of testing it for MIRI and MUSE, and
most other instruments should be quite easy to add as well.


## What does it do (and how)?

This package contains of a few main components of functionality, plus some plumbing to make them play nice together. 

- A family of `SpectralCube` types which handle the data and metadata.
- A continuum subtraction function (see description below)
- The functionality for fitting emission lines, either -per-spaxel, or for a number of arbitrary spatial bins/fragments.
- A wrapper around the `VoronoiBinning` package, which is an optional dependency. 

The workhorses of this package are the various types `NIRSpecCube`, `MUSECube`,
etc.; as well as the functions `cont_subt()`, `fit_cube()`, and
`fit_spectrum_from_subcube()`.


## Usage

### Installing

This package is not yet registered to install through the Julia package manager
(maybe later). You can install the GitHub `main` branch by running:

```julia 
julia> import Pkg; Pkg.add("https://github.com/thriveth/CubeFitter.jl")
```


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

The function `fit_cube()` will return a dictionary of $N\times M \times 5$
arrays, one for each fitted line, as well as one for the 0th to 2nd moment, the
spaxel-wise redshift, the line FWHM, and the fit statistics. Each such array
will have the measured values as the first slice, the parameter error in the
second slice, and S/N in the third. Fourth and fifth slice in this array
contain the numerically coadded line flux and its associated standard error.
Each of these entries will be saved as an individual HDU in the output FITS
file.

The function `fit_cube` is mainly a convenience function enabling one to cycle
through either individual spaxels or a map of numbered bins, and extracting and
fitting a spectrum for every mask or spaxel, and filling in best-fit fitted
values of line fluxes and kinematic parameters and their associated standard
errors and fit statistcs into the corresponding spaxel or fragment.

In addition to the line fits, the routine also measures the flux with standard
errors by numerical integration in each spaxel; this is saved together with the
flux from fitting in the output. 

One can also use the function `fit_spectrum_from_subcube()` directly. This
function allows to extract a spectrum from the cube based on an x- and y range
passed to it (this must always be a range; if only one spaxel is wanted, give
the range for the spaxel (_i_, _j_) as `xrange=i:i, yrange=j:j`). The function
extracts a simple, un-weighted spectrum over the given (_x, y_) range, runs the
fit, and return a dictionary of the fit results and statistics, along with the
extracted spectrum, errors, and wavelength range for convenience.
Alternatively, one can pass a mask (`true` for included spaxels, `false` for
excluded ones. See the function docstring to learn more.


### Data format

The function `fit_cube()` expects the data to already be continuum subtracted,
but otherwise saved in the same data format as the final pipeline products of a
given instrument; this can be done using the `cont_subt()` function as
described above. Using JWST/NIRSpec as an example, this means that the data
cubes should be saved as a FITS file with an empty primary HDU, and the flux
and error cubes saved as the first and second extensions, respectively. 

When instantiating the cube, it is possible to pass the data and error
extension indices or names as keyword arguments. Examples: 

```julia
julia> thecube = NIRSpecCube("/path/to/fits/file.fits", data_ext=0, err_ext=3)
julia> thecube = NIRSpecCube("/path/to/fits/file.fits", data_ext="SCI", err_ext="ERR")
```


### Output format

The fitting function (`fit_cube`) returns a dictionary of maps plus a few other
properties. Those maps are the *redshift* (common to all lines), the *fwhm*
(also common to all lines), as well as an entry for each line measured. These
entries each contain an $N \times M \times 5$ array, where $N \times M$ are the
spatial dimensions of the cube, and each consists of five layers: 

1. Measured line flux in each spaxel from line fitting.
2. Standard errors of layer 1. 
3. S/N from line fitting (mainly for convenience), 
4. Numerical flux from each line, and 
5. Standard errors of layer 4. 

This output dictionary can be saved to a FITS file using the function
`write_maps_to_fits`, which writes each of the map entries to a separate FITS
extension HDU named for the dictionary key, and containing the output array in
its data part. The output data file follows the convention of having the first,
primary, extension be empty of data and containing the most elaborate header
(which at this point is basically just a copy of the header of the input
datacube). 


### Specify lines to fit

If nothing is specified in the call to `fit_cube` or
`fit_spectrum_from_subcube`, these functions will attempt to fit all lines in
the line list loaded along with the datacube, but will skip any that either is
outside the wavelength range of the cube, or has flux below the value of `min_snr`.

It is also possible to specify a list of lines to fit at the time of calling
these functions, with the optional argument `line_selection`, taking a list of
line names in String of Symbol format.

If no kinematic template lines are selected (see below), the lines in this list
will be fitted simultaneously with shared kinematic parameters.

If kinematic template lines are specified, then the remaining lines will be
modeled one by one.


### Use selection of lines as kinematics template

the `fit_cube` function can as an optional argument take a list of lines to use
for a kinematics template. In this case, the lines of this list will first be
fit simultaneously as above. Afterwards, the remaining lines will be fit one by
one with their kinematics fixed to the one found from the template lines, and
the flux left as the only free parameters. 

```julia
julia> cscube = NIRSpecCube("/path/to/cont_subtracted/datacube.fits", "g140m"; z_init=0.76)
julia> out = fit_cube(cscube, kinematics_from_lines=[:HI_4861, :OIII_4959, :OIII_5007])
julia> quicklook_slice(out, :OIII_5007, norm=sqrt, cmap=:inferno)
```

**NB!** Lines that are part of blended features should be included in this list
to be modeled properly, even if they are not particularly strong. I will
implement something more elegant if I find the time (pull requests are always welcome!) 


### Quicklook functions

`CubeFitter` includes three functions to quickly view fitting output for
quality control. These use `Plots.jl` and whichever backend is set up as the
default. See the docstrings for each function in the Julia REPL for
documentation. The functions are:

- **`quicklook_slice`**: Views the 2D line and kinematics maps output by
  `fit_cube`. 

- **`quicklook_fit_result_dict`**: Previews the output of a single spectrum fit
  generated by `fit_spectrum_from_subcube`. This is useful to see if the fit in
  a specific spaxel or set of spaxels makes sense.


- **`quicklook_model`**: Previews a model generated by the function
  `build_model`. See the docstrings of these functions for directions.


### Fitting from fragmentation map

As mentioned above, `fit_cube` and `fit_spectrum_from_subcube` can take an
arbitrary (set of) mask(s) as input, for which to extract spectrum and perform
fitting. This set of masks can be a fragmentation map, a map of Voronoi bins
(see under "extras" below), or created in any arbitrary way.

#### Map/mask format

- **`fit_cube`**: A 2D array of **integers** of the same dimensions as the
  spatial dimensions of the data cube in question. Every number is interpreted
  as a separate fragment. 
- **`fit_spectrum_from_subcube`**: A 2D BitArray (created as an array of type
  Bool), where excluded spaxels get the value `false`/0, and included spaxels are
  `true`/1. 


### Extras

#### Voronoi binning 

If the package
![VoronoiBinning](https://github.com/Michael-Reefe/VoronoiBinning.jl) is
installed and loaded, it activates the CubeFitter function `voronoi_bin_slice`. 

This function takes as argument a fit result dictionary (output from
`fit_cube`) and a slice key (any slice that contains both signal and noise is
valid), a target S/N value, and a few other arguments (see the docstring for
more information). The function outputs the resulting 2D map of bins as well as
the binned value of the selected quantity to bin for. The latter comes as the
type `Matrix{Measurement{Float64}}`, use the package `Measurements` and its
functions `value` and `uncertainty` to separate out these two quantities as
individual arrays. 

## Screenshots

There really isn't much to look at but here it is:

![Screenshot of CubeFitter in action](./Screenshots/CubeFitter.png)


# Development

To hack on your own branch, clone the repository to your preferred location;
then activate it. It is a good idea to use `Revise.jl` to have changes
impolemented and precompiled on-the-go: 

```julia
julia> using Revise   # Optional but a good idea
julia> import Pkg
julia> Pkg.activate("/path/to/CubeFitter.jl")
julia> Pkg.instantiate()  # Install dependencies, needs only be done once.
julia> using CubeFitter
```

Alternatively, you can enter the `Pkg>` prompt, run `activate /path/to/CubeFitter.jl`,
then (first time) `instantiate()`. Press Backspace to return to the normal `julia>`
prompt, and run `using CubeFitter`.


## Planned features / whishlist

In order of approximate priority: 

- [x] Include continuum subtraction functionality in the package.
- [x] Make it possible to add a second and perhaps third kinematic component. 
- [x] Write quicklook-functions allowing to quickly view the fit outputs with
  minimum input. But still be tinker-friendly, don't hide stuff from the user.
- [x] Allow to measure flux/upper limits numerically when S/N threshold is not
  met (now does this all the time, whether S/N threshold is met or not - it is
  computationally cheap and simpler this way. .
- [ ] Spectrum extraction and fitting from arbitrary masks, e.g. fragmentation
  maps.
- [ ] Implement adaptive binning.
- [ ] Allow user to fix ratio between lines of doublets with shared upper
  levels.
- [ ] Create interface to select lines to always fit together (useful for
  blended features).
- [ ] Test and ensure the `MUSECube` and `MIRICube` structs actually work as
  advertised
- [ ] Add support for more instruments. Suggestions welcome (especially if
  accompanied with a suitable test dataset).


# Other projects

- [`Loki.jl`](https://github.com/Michael-Reefe/Loki.jl) is another (and way
  more advanced and ambitious) Julia package for IFU data than this one. Loki
  seems to implement much more astrophysics and much more complex models than
  the simple non-parametric continuum subtraction and Gaussian peaks of
  CubeFitter. On the other hand, with the higher specialization also comes a
  more limited scope in terms of which targets you could apply it on; and I
  find that CubeFitter is simpler to use for its different scope. CubeFitter
  also is written in pure Julia and does not depend on any Python packages
  through `PyCall`. 

  If you want to implement physical per-spaxel dust models, PAH features and
  stellar kinematics in your model, then Loki is probably what you want. If on
  the other hand, you are simply looking for a 1- or 2-Gaussian component model
  of your lines with shared kinematics by default, CubeFitter might be right
  for you. 
