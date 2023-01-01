using EverySingleStreet
using Documenter

DocMeta.setdocmeta!(EverySingleStreet, :DocTestSetup, :(using EverySingleStreet); recursive=true)

makedocs(;
    modules=[EverySingleStreet],
    authors="Ole Kröger <o.kroeger@opensourc.es> and contributors",
    repo="https://github.com/Wikunia/EverySingleStreet.jl/blob/{commit}{path}#{line}",
    sitename="EverySingleStreet.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://Wikunia.github.io/EverySingleStreet.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/Wikunia/EverySingleStreet.jl",
    devbranch="main",
)
