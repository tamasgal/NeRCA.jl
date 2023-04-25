# NeRCA.jl Package

## Introduction

Welcome to the `NeRCA.jl` (**Ne**utrino **R**esearch with **C**osmics in the **A**byss) package documentation. This piece of software was written
for those who like to explore the KM3NeT world with a new, modern and ultra fast
programming language designed for scientific computing: [Julia](https://www.julialang.org).
The I/O for KM3NeT related data is provided by [KM3io.jl](https://git.km3net.de/common/KM3io.jl).

## Installation

`NeRCA.jl` was initially registered on the officially Julia package registry but since version 0.11.0, it resides in the [KM3NeT Julia Registry](https://git.km3net.de/common/julia-registry). To install the latest version, you first need to add the KM3NeT Julia Registry to your local registries. Follow the instructions in its README or simply do

    git clone https://git.km3net.de/common/julia-registry ~/.julia/registries/KM3NeT

After that, you can install `NeRCA.jl` just like any other Julia package with `Pkg`:

```julia-repl
using Pkg
Pkg.add("NeRCA")
```

## General information

`NeRCA.jl` started as a hobby project and contained all kinds of KM3NeT related stuff. In 2023, many of the functionalities were extracted to distinct packages. Here are some examples:

- [KM3io.jl](https://git.km3net.de/common/KM3io.jl): general I/O (reading ROOT, DETX and other KM3NeT files)
- [KM3DB.jl](https://git.km3net.de/tgal/KM3DB.jl): accessing the KM3NeT Oracle Database (StreamDS, run setups, detector calibrations, ...)
- [KM3Acoustics.jl](https://git.km3net.de/acoustics/KM3Acoustics.jl): Acoustics stuff, simulation of acoustic events, calibration procedures, ..
- [KM3NeTTestData.jl](https://git.km3net.de/km3py/km3net-testdata): A collection of KM3NeT related sample files (used for tests)
