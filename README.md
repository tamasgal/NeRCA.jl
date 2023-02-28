# NeRCA

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://tamasgal.github.io/NeRCA.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://tamasgal.github.io/NeRCA.jl/dev)
[![Build Status](https://github.com/JuliaPhysics/Corpuscles.jl/workflows/CI/badge.svg)](https://github.com/tamasgal/NeRCA.jl/actions)
[![Codecov](https://codecov.io/gh/tamasgal/NeRCA.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/tamasgal/NeRCA.jl)
[![Build Status](https://api.cirrus-ci.com/github/tamasgal/NeRCA.jl.svg)](https://cirrus-ci.com/github/tamasgal/NeRCA.jl)

NeRCA (**Ne**utrino **R**esearch with **C**osmics in the **A**byss) is a Julia
package to access and analyse KM3NeT related data.
It can read the official KM3NeT ROOT and formats, calibration information (DETX
files) and also live data from running detectors.
The ROyFit reconstruction is also included and can be found in `src/fit.jl`.
