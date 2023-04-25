# NeRCA

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://tgal.pages.km3net.de/NeRCA.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://tgal.pages.km3net.de/NeRCA.jl/dev)
[![Build Status](https://git.km3net.de/tgal/NeRCA.jl/badges/main/pipeline.svg)](https://git.km3net.de/tgal/NeRCA.jl/pipelines)
[![Coverage](https://git.km3net.de/tgal/NeRCA.jl/badges/main/coverage.svg)](https://git.km3net.de/tgal/NeRCA.jl/commits/main)

NeRCA (**Ne**utrino **R**esearch with **C**osmics in the **A**byss) is a Julia
package to access and analyse KM3NeT related data.
It can read the official KM3NeT ROOT and formats, calibration information (DETX
files) and also live data from running detectors.
The ROyFit reconstruction is also included and can be found in `src/fit.jl`.
