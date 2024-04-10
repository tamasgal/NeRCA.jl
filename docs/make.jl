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
    checkdocs = :exports,
    pages = [
        "Home" => "index.md",
        "Hits" => "hits.md",
        "FibonacciFit" => "ffit.md",
    ],
    repo = Documenter.Remotes.URL(
        "https://git.km3net.de/tgal/NeRCA.jl/blob/{commit}{path}#L{line}",
        "https://git.km3net.de/tgal/NeRCA.jl"
    ),
)

deploydocs(;
  repo = "git.km3net.de/tgal/NeRCA.jl",
  devbranch = "main",
  push_preview=true
)
