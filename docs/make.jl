using PreCoIBM
using Documenter

makedocs(;
    modules=[PreCoIBM],
    authors="Hetao Z.",
    repo="https://github.com/HetaoZ/PreCoIBM.jl/blob/{commit}{path}#L{line}",
    sitename="PreCoIBM.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://HetaoZ.github.io/PreCoIBM.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/HetaoZ/PreCoIBM.jl",
)
