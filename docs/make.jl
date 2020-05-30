using SurveyDataWeighting
using Documenter

makedocs(;
    modules=[SurveyDataWeighting],
    authors="Graham Stark",
    repo="https://github.com/grahamstark/SurveyDataWeighting.jl/blob/{commit}{path}#L{line}",
    sitename="SurveyDataWeighting.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://grahamstark.github.io/SurveyDataWeighting.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/grahamstark/SurveyDataWeighting.jl",
)
