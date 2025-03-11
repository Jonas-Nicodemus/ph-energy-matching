using DrWatson
@quickactivate "ph-energy-matching"
using Documenter

using PortHamiltonianSystems, PortHamiltonianModelReduction, EnergyMatching, QuadraticOutputSystems

@info "Building Documentation"
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
        "QuadraticOutputSystems" => "QuadraticOutputSystems.md",
        "PortHamiltonianSystems" => "PortHamiltonianSystems.md",
        "PortHamiltonianModelReduction" => "PortHamiltonianModelReduction.md",   
        "API" => "API.md",
    ],
)

@info "Deploying Documentation"
deploydocs(;
    repo="github.com/Jonas-Nicodemus/ph-energy-matching.git",
    target = "build",
    push_preview = true,
    devbranch="main",
)
