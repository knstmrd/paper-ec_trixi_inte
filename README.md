
[![License: MIT](https://img.shields.io/badge/License-MIT-success.svg)](https://opensource.org/licenses/MIT)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.10869717.svg)](https://doi.org/10.5281/zenodo.10869717)

This repository contains information and code to reproduce the results presented in the article
```bibtex
@online{oblapenko2024entropyconservative,
  title={Entropy-conservative high-order methods for high-enthalpy gas flows},
  author={Oblapenko, Georgii and Torrilhon, Manuel},
  year={2024},
  eprint={2403.16882},
  archivePrefix={arXiv},
  primaryClass={physics.flu-dyn},
  eprintclass={physics.flu-dyn},
  journal={arXiv preprint 2403.16882},
  doi={10.48550/arXiv.2403.16882}
}
```
If you find these results useful, please cite the article mentioned above. If you use the implementations provided here, please also cite this repository as
```bibtex
@misc{oblapenko2024entropyconservativeRepro,
  title={Reproducibility repository for "Entropy-conservative high-order methods for high-enthalpy gas flows"},
  author={Oblapenko, Georgii and Torrilhon, Manuel},
  year={2024},
  howpublished={\url{https://github.com/knstmrd/paper-ec_trixi_inte}},
  doi={10.5281/zenodo.10869716}
}
```

## Abstract

A framework for numerical evaluation of entropy-conservative volume fluxes in gas flows with internal energies is developed, for use with high-order discretization methods. The novelty of the approach lies in the ability to use arbitrary expressions for the internal degrees of freedom of the constituent gas species. The developed approach is implemented in an open-source discontinuous Galerkin code for solving hyperbolic equations. Numerical simulations are carried out for several model 2-D flows and the results are compared to those obtained with the finite volume-based solver DLR TAU.

## Reproducing the results

### Installation

To download the code using `git`, use 

```bash
git clone git@github.com:knstmrd/paper-ec_trixi_inte.git
``` 


To instantiate the environment execute the following two commands:
```bash
cd paper-ec_trixi_inte
julia --project=. -e 'using Pkg; Pkg.resolve(); Pkg.instantiate()'
```

Note that the results are obtained using Julia 1.10.5 and Trixi 0.8.8.
Thus, you might need to install the [old Julia 1.10.5 release](https://julialang.org/downloads/oldreleases/) first
and *replace* the `julia` calls from this README with
`/YOUR/PATH/TO/julia-1.10.5/bin/julia`


### The code

The implementation of the equations and necessary functions (pre-processing, computation of primitive variables,
evaluation of fluxes, etc.) is in the `src/compressible_euler_2d_intenergy.jl` file.
Some models for the rotational and vibrational energies are provided in the `src/internal_energy_models.jl` file.

### Running the examples

The scripts for the different test cases are located in the `examples` directory.

To execute them provide the path (if inside the `paper-ec_trixi_inte` directory):

```bash
julia --project=. examples/elixir_euler_vibr_aho_cylinder.jl
```

30x30 and 60x60 grids are possible, the grid size can be changed by changing the `Nx = ` line in the example file.
The polynomial degree can be changed by changing the `polydeg = ` line in the example file. 

## Authors

* [Georgii Oblapenko](https://acom.rwth-aachen.de/the-lab/team-people/name:georgii_oblapenko) (Corresponding Author)
* [Manuel Torrilhon](https://www.acom.rwth-aachen.de/the-lab/team-people/name:manuel_torrilhon)

## Disclaimer

Everything is provided as is and without warranty. Use at your own risk!

## Acknowledgments

This work has been supported by the German Research Foundation within the research unit DFG-FOR5409. 
