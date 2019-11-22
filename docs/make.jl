using Documenter, NeRCA

makedocs(;
    modules=[NeRCA],
     format = Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true"),
    pages=[
        "Home" => "index.md",
        "I/O" => "io.md",
        "Hits" => "hits.md",
        "ROyFit" => "fit.md",
    ],
    repo="https://github.com/tamasgal/NeRCA.jl/blob/{commit}{path}#L{line}",
    sitename="NeRCA.jl",
    authors="Tamas Gal",
    assets=String[],
)

deploydocs(;
    repo="github.com/tamasgal/NeRCA.jl",
)
