# Energy matching in reduced passive and port-Hamiltonian systems

This repository contains the code for the paper Energy matching in reduced passive and port-Hamiltonian systems.
The goal is to find low-dimensional port-Hamiltonian (pH) models that not only match the input-output dynamic of a full order model (FOM), but also its energy (Hamiltonian) trajectory.

<!-- TABLE OF CONTENTS -->
<details open="open">
  <summary><h2 style="display: inline-block">Table of Contents</h2></summary>
  <ol>
    <!-- <li>
      <a href="#citing">Citing</a>
    </li> -->
    <li>
      <a href="#installation">Installation</a>
    </li>
    <li><a href="#usage">Usage</a></li>
    <li><a href="#license">License</a></li>
    <li><a href="#contact">Contact</a></li>
  </ol>
</details>

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


<!-- USAGE EXAMPLES -->
## Usage

The executable script `main.jl` is located in the `scripts` directory. 
It performs the following steps:
1. Set up a pH full order model and declare the methods (`Reductors` and (Energy-)`Matcher`) to run.
2. Run the methods.
3. Evaluate the ROMs.
4. Analyze the results.

Note that for the exact reproduction of the results in the paper for the poroelasticity model, 
the solution of the positive-real algebraic Riccati equation from MATLAB's `icare` is required.

This requires a running version of MATLAB and the package [MATLAB.jl](https://github.com/JuliaInterop/MATLAB.jl).
Then you need to uncomment
- the `using MATLAB` line in `src/PortHamiltonianSystems/PortHamiltonianSystems.jl`.
- in the `prgram` function in `src/PortHamiltonianSystems/gramians.jl`, the MATLAB related lines must be uncommented.

## License
Distributed under the MIT License. See `LICENSE` for more information.

## Contact
Jonas Nicodemus - jonas.nicodemus@simtech.uni-stuttgart.de

Tobias Holicki - tobias.holicki@imng.uni-stuttgart.de\
Paul Schwerdtner - paul.schwerdtner@nyu.edu\
Benjamin Unger - benjamin.unger@simtech.uni-stuttgart.de

Project Link: [https://github.com/Jonas-Nicodemus/ph-energy-matching][project-url]

[license-shield]: https://img.shields.io/github/license/Jonas-Nicodemus/ph-energy-matching.svg?style=for-the-badge
[license-url]: https://github.com/Jonas-Nicodemus/ph-energy-matching/blob/main/LICENSE
[project-url]:https://github.com/Jonas-Nicodemus/ph-energy-matching/
[docs-shield]:https://img.shields.io/badge/docs-online-blue.svg?style=for-the-badge
[docs-url]:https://jonas-nicodemus.github.io/ph-energy-matching/