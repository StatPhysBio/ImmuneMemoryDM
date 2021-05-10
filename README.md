# ImmuneMemoryDM
This repository contains the code and analysis associated with the manuscript:

Schnaack, Nourmohammad, [Optimal evolutionary decision-making to store immune memory](https://elifesciences.org/articles/61346). eLife 2021

It allows reproduction and a simplified visualization of the numerical results reported in the manuscript.

## Dependencies

The code is written in [Julia](https://julialang.org) and depends on several packages. Below we give a list of version numbers of the packages for which the code is known to run.
- Julia 1.5.3
- Distributions 0.24.6
- SpecialFunctions 1.1.0
- QuadGK 2.4.1
- FileIO 1.4.5
- IJulia 1.23.1
- Plots 1.9.1

## Usage

The repository includes the source code to generate and analyze data as well as jupiter notebooks that can be used to produce all figures used in the manuscript. To generate your own figures, you can, for example, install [IJulia](https://github.com/JuliaLang/IJulia.jl) and run 
```bash
using IJulia 
notebook()
```
in your Julia terminal to launch the IJulia notebook in your browser. Then select one of the notebooks to reproduce the figures of the manuscript.
In order to reduce the computation time, the default simulation parameters are set to fewer repetitions and a smaller resolution compared to the manuscript, but the results capture the reported behavior. To get comparable results set the parameters to the values given in the manuscript. Note: As most of the simulations are stochastic you generally do not expect precisely equivalent plots.

## Contact

Any issues or questions should be addressed to [me](mailto:oskar.schnaack@ds.mpg.de).
