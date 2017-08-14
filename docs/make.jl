using Documenter, KM3NeT

makedocs(modules=[KM3NeT],
         doctest=true,
         format=:html,
         sitename = "KM3NeT.jl",
         pages = [
             "index.md",
             "Calibration" => "calibration.md",
             #= "Subsection" => [ =#
             #=     ... =#
             #= ] =#
         ]
         )

deploydocs(deps   = Deps.pip("mkdocs", "python-markdown-math"),
    repo = "github.com/tamasgal/KM3NeT.jl.git",
    julia  = "0.6.0",
    osname = "linux")
