using Documenter, NeRCA

makedocs(modules=[NeRCA],
         doctest=true,
         format = Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true"),
         authors = "Tamas Gal",
         sitename = "NeRCA.jl",
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
#=     repo = "github.com/tamasgal/NeRCA.jl.git", =#
#=     julia  = "0.7.0", =#
#=     osname = "linux") =#
