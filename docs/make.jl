using Documenter, NeRCA

makedocs(;
    modules = [NeRCA],
    sitename = "NeRCA.jl",
    authors = "Tamas Gal",
    format = Documenter.HTML(;
        assets = ["assets/custom.css"],
        sidebar_sitename = false,
        collapselevel = 4,
        warn_outdated = true,
    ),
    pages = [
        "Home" => "index.md",
        "I/O" => "io.md",
        "Hits" => "hits.md",
        "ROyFit" => "fit.md",
    ],
    repo = "https://git.km3net.de/tgal/NeRCA.jl/blob/{commit}{path}#L{line}",
)

deploydocs(;
  repo = "git.km3net.de/tgal/NeRCA.jl",
  devbranch = "main",
  push_preview=true
)
