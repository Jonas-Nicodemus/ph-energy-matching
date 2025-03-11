# Energy matching in reduced passive and port-Hamiltonian systems
This repository contains the code for the paper [Energy matching in reduced passive and port-Hamiltonian systems](https://arxiv.org/abs/2309.05778).
The goal is to obtain low-dimensional port-Hamiltonian (pH) models that effectively approximate both the input-output dynamics and the energy (Hamiltonian) of a full-order model (FOM).

## Citing
If you use this project for academic work, please consider citing our
[publication](https://arxiv.org/abs/2309.05778):

    T. Holicki, J. Nicodemus, P. Schwerdtner, and B. Unger
    Energy matching in reduced passive and port-Hamiltonian systems
    ArXiv e-print 2309.05778, 2023.

## Installation
This code base is using the [Julia Language](https://julialang.org/) and
[DrWatson](https://juliadynamics.github.io/DrWatson.jl/stable/)
to make a reproducible scientific project named
> ph-energy-matching

To (locally) reproduce this project, do the following:

1. Download this code base. Notice that raw data are typically not included in the
   git-history and may need to be downloaded independently.
2. Open a Julia console and do:
   ```
   julia> using Pkg
   julia> Pkg.add("DrWatson") # install globally, for using `quickactivate`
   julia> Pkg.activate("path/to/this/project")
   julia> Pkg.instantiate()
   ```

This will install all necessary packages for you to be able to run the scripts and
everything should work out of the box, including correctly finding local paths.

You may notice that most scripts start with the commands:
```julia
using DrWatson
@quickactivate "ph-energy-matching"
```
which auto-activate the project and enable local path handling from DrWatson.

## Usage
The executable script `main.jl` is located in the `scripts` directory. 
It performs the following steps:
1. Set up the experiment.
2. Apply minimal realization (Kalman-like decomposition).
2. Declare and run the methods (`Reductors` and `EnergyMatcher`).
3. Evaluate the ROMs.
4. Analyze the results.

## Project structure 
This project contains four packages:
- [EnergyMatching.jl](https://jonas-nicodemus.github.io/ph-energy-matching/dev/EnergyMatching/): Contains the methods for solving the energy matching problem.
- [PortHamiltonianModelReduction.jl](https://jonas-nicodemus.github.io/ph-energy-matching/dev/PortHamiltonianModelReduction/): Contains two structure preserving model reduction algorithms for port-Hamiltonian systems, `phirka` and `prbt`.
- [PortHamiltonianSystems.jl](https://jonas-nicodemus.github.io/ph-energy-matching/dev/PortHamiltonianSystems/): Contains methods for the analysis of port-Hamiltonian systems, as well as the `PortHamiltonianStateSpace` data type.
- [QuadraticOutputSystems.jl](https://jonas-nicodemus.github.io/ph-energy-matching/dev/QuadraticOutputSystems/): Contains methods for the analysis of linear dynamical systems with quadratic output, as well as the `QuadraticOutputStateSpace` data type.

## License
Distributed under the MIT License. See `LICENSE` for more information.

## Contact
Jonas Nicodemus - jonas.nicodemus@simtech.uni-stuttgart.de

Tobias Holicki - tobias.holicki@imng.uni-stuttgart.de\
Paul Schwerdtner - paul.schwerdtner@nyu.edu\
Benjamin Unger - benjamin.unger@simtech.uni-stuttgart.de