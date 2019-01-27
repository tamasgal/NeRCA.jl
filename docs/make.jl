using Documenter, KM3NeT

makedocs(modules=[KM3NeT],
         doctest=true,
         format = Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true"),
         authors = "Tamas Gal",
         sitename = "KM3NeT.jl",
         pages = [
             "index.md",
             "Calibration" => "calibration.md",
             #= "Subsection" => [ =#
             #=     ... =#
             #= ] =#
         ]
         )
#=  =#
#= deploydocs(deps   = Deps.pip("mkdocs", "python-markdown-math"), =#
#=     repo = "github.com/tamasgal/KM3NeT.jl.git", =#
#=     julia  = "0.7.0", =#
#=     osname = "linux") =#
