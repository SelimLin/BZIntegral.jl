using BZIntegral
using Documenter

DocMeta.setdocmeta!(BZIntegral, :DocTestSetup, :(using BZIntegral); recursive=true)

makedocs(;
    modules=[BZIntegral],
    authors="Yihao",
    repo="https://github.com/SelimLin/BZIntegral.jl/blob/{commit}{path}#{line}",
    sitename="BZIntegral.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://SelimLin.github.io/BZIntegral.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/SelimLin/BZIntegral.jl",
    devbranch="main",
)
