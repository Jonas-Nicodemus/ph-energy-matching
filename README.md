[![arXiv][arxiv-shield]][arxiv-url]
[![DOI][doi-shield]][doi-url]
[![Documentation Status][docs-shield]][docs-url]
[![MIT License][license-shield]][license-url]

# [Energy matching in reduced passive and port-Hamiltonian systems][arxiv-url]
This repository contains the code for the paper [Energy matching in reduced passive and port-Hamiltonian systems][arxiv-url].
The goal is to obtain low-dimensional port-Hamiltonian (pH) models that effectively approximate both the input-output dynamics and the energy (Hamiltonian) of a full-order model (FOM).

<!-- TABLE OF CONTENTS -->
<details open="open">
  <summary><h2 style="display: inline-block">Table of Contents</h2></summary>
  <ol>
    <li>
      <a href="#citing">Citing</a>
    </li>
    <li>
      <a href="#installation">Installation</a>
    </li>
    <li><a href="#usage">Usage</a></li>
    <li><a href="#license">License</a></li>
    <li><a href="#contact">Contact</a></li>
  </ol>
</details>

## Citing
If you use this project for academic work, please consider citing our
[publication][arxiv-url]:

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

## License
Distributed under the MIT License. See `LICENSE` for more information.

## Contact
Jonas Nicodemus - jonas.nicodemus@simtech.uni-stuttgart.de

Tobias Holicki - tobias.holicki@imng.uni-stuttgart.de\
Paul Schwerdtner - paul.schwerdtner@nyu.edu\
Benjamin Unger - benjamin.unger@simtech.uni-stuttgart.de

Project Link: [https://github.com/Jonas-Nicodemus/ph-energy-matching][project-url]

[doi-shield]: https://img.shields.io/badge/DOI-10.5281%20%2F%20zenodo.8335231-blue.svg?style=for-the-badge
[doi-url]: https://doi.org/10.5281/zenodo.8335231
[arxiv-shield]: https://img.shields.io/badge/arXiv-2204.13474-b31b1b.svg?style=for-the-badge
[arxiv-url]: https://arxiv.org/abs/2309.05778
[license-shield]: https://img.shields.io/github/license/Jonas-Nicodemus/ph-energy-matching.svg?style=for-the-badge
[license-url]: https://github.com/Jonas-Nicodemus/ph-energy-matching/blob/main/LICENSE
[project-url]:https://github.com/Jonas-Nicodemus/ph-energy-matching/
[docs-shield]:https://img.shields.io/badge/docs-online-blue.svg?style=for-the-badge
[docs-url]:https://jonas-nicodemus.github.io/ph-energy-matching/