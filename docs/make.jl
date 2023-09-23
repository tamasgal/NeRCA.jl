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
    warnonly = [:missing_docs],
    pages = [
        "Home" => "index.md",
        "I/O" => "io.md",
        "Hits" => "hits.md",
        "ROyFit" => "fit.md",
    ],
    repo = Documenter.Remotes.URL(
        "https://git.km3net.de/tgal/NeRCA.jl/blob/{commit}{path}#L{line}",
        "https://git.km3net.de/tgal/NeRCA.jl"
    ),
)

deploydocs(;
  repo = "git.km3net.de/tgal/NeRCA.jl",
  devbranch = "master",
  push_preview=true
)
