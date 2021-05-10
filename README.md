# ImmuneMemoryDM
This repository contains the code and analysis associated with the manuscript
Schnaack, Nourmohammad, [Optimal evolutionary decision-making to store immune memory](https://elifesciences.org/articles/61346). eLife 2021
It allows reproduction and a simplefied visulaization of the numerical results reported in the manuscript.

## Dependencies

The code is written in Julia and depends on a number of packeges. Below we give a list of version numbers of the packages for which the code is known to run.
- Julia 1.5.3
- Distributions 0.24.6
- SpecialFunctions 1.1.0
- QuadGK 2.4.1
- FileIO 1.4.5
- IJulia 1.23.1
- Plots 1.9.1

### Usage

The reprosetory includes to source code to generate and analyze data as well as jupiter notebooks that can be used to produce all figures used in the manuscript. To gernate your own figures you can for exsample install [IJulia](https://github.com/JuliaLang/IJulia.jl) and run 
```bash
using IJulia 
notebook()
```
in your Julia terminal to launch the IJulia notebook in your browser. Then select one of the notebooks to reproduce the figures of the manuscript. 

In order to reduce the computation time, the default simulation perameters are set less repetitions and to a smaller resolution compared to the mamuscript but the results capture the behavior reported in the manuscript. To get the compareble results set the perameters to the values given in the manuscript. Note: As most the simulations are stochastic you generally do not expect precisely equivalent plots.

