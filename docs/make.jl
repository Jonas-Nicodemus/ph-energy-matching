using Documenter

using PortHamiltonianSystems, PortHamiltonianModelReduction, EnergyMatching, QuadraticOutputSystems

makedocs(;
    modules=[PortHamiltonianSystems, PortHamiltonianModelReduction, EnergyMatching, QuadraticOutputSystems],
    authors="Jonas Nicodemus <jonas.nicodemus@icloud.com> and contributors",
    repo="https://github.com/Jonas-Nicodemus/ph-energy-matching/blob/{commit}{path}#{line}",
    sitename="ph-energy-matching",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://Jonas-Nicodemus.github.io/ph-energy-matching",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "EnergyMatching" => "EnergyMatching.md",
        "PortHamiltonianModelReduction" => "PortHamiltonianModelReduction.md",
        "PortHamiltonianSystems" => "PortHamiltonianSystems.md",
        "QuadraticOutputSystems" => "QuadraticOutputSystems.md",
        "API" => "API.md",
    ],
)

deploydocs(;
    repo="github.com/Jonas-Nicodemus/ph-energy-matching",
    devbranch="main",
)
