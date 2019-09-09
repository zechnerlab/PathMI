# PathMI.jl
Julia code for calculating the path mutual information between complete trajectories of two molecular components.

The provided code is customized to reproduce the case studies of the paper "Path mutual information for a class of biochemical reaction networks" (http://arxiv.org/abs/1904.01988) and is not meant to be a general tool. 

For further information, contact Lorenzo Duso (duso@mpi-cbg.de) & Christoph Zechner (zechner@mpi-cbg.de)

# Installation and requirements

Julia version greater than 1.0 and the package PyPlot are required.

Install PyPlot by running in the Julia REPL
```julia
import Pkg
Pkg.add("PyPlot")
using PyPlot
```

# Reproducing the figures of the paper

To get started, include the module PathMI.jl
```julia
include("PathMI.jl")
```
Each figure in the paper can be readily reproduced by a dedicated function: 
```julia
PathMI.figure1a()  # for figure 1a
PathMI.figure1b()  # for figure 1b
PathMI.figure1c()  # for figure 1c
PathMI.figure1d()  # for figure 1d
PathMI.figure2a()  # for figure 2a
PathMI.figure2b()  # for figure 2b
```
If it is desired to get precisely the identical figures as reported in the paper, it is possible to pass the argument ` paper=true `, for example `PathMI.figure1b(paper=true)` . This command sets the seed of the random number generator to the same value used when making the figures.

# Advanced usage

This section exaplains how to independently perform simulations and further play with the case studies provided in the paper. 
First, it is necessary to initialize the system. Run
```julia
S = PathMI.initialize1() # for the first case study (mRNA-Protein transcription)
S = PathMI.initialize2() # for the second case study (mRNA-Protein transcription, transient behavior)
S = PathMI.initialize3() # for the third case study (Lotka-Volterra model)
```
to store in the variable `S` the system settings of the respective case study. 

If desired, the default initial condition x(0)=(A(0),B(0)) for the two chemical species can be changed by modifying the entries of `S.x0` . For instance,
```julia
S.x0[1] = 100  # sets A(0)=100
S.x0[2] = 25   # sets B(0)=25
```
If desired, the reaction rate constants can be modified by acting on the entries of `S.C` . Note that, independently of the case study, the network initialized in `S` reflects the following indexing
```
1:    0 -> A
2:    A -> 2A
3:    A -> 0
4:  A+B -> B
5:    A -> A+B
6:  A+B -> A+2B
7:    B -> 0
```
For instance, `S.C[3] = 0.5` sets to 0.5 the rate constant of the reaction with index 3, i.e. A -> 0 . Note that setting the rate to `0.` is equivalent to removing the corresponding reaction. Alternatively, you can quickly initialize a system with zero rates and zero initial condition with
```
S = PathMI.create_network() 
``` 
and afterwards specify the initial condition and activate the reactions you want to consider by setting a rate greater than zero.

The two main functions of PathMI.jl are:
- SSA
```julia 
SSA(S::System, t::Vector{Float64})
```
performs one SSA simulation for the system `S` combined with the online computation of the associated sample of path mutual information. By default, it returns a vector with the path-MI sample computed at each entry of `t` and an object of type `::Path` which stores all the data generated by the simulation (e.g. firing times, reaction indexes, species trajectories, latent conditional moments, etc...). The keyword `output_path=false` can be passed if it is desired to return only the path mutual information sample.
Usage example:
```julia
S = PathMI.initialize1()  # for instance, initialize the first case study
t = collect(0:1.:200)     # time points where to output the path mutual information

# return path-MI sample and all simulation data
pmi_sample, Xt = PathMI.SSA(S,t)

# return only the path-MI sample
pmi_sample2 = PathMI.SSA(S,t,output-path=false)

```

- estimate
```julia
estimate(S::System, t::Vector{Float64}, N::Int64)
```
performs `N` SSA simulations of system `S` and returns the Monte Carlo estimate of the path mutual information. It returns the following four vectors: path-MI estimate, its standard deviation, the average pathMI rate, and its standard deviation. All of them are evaluated at the time points specified by `t`.
Usage example:
```julia
S = PathMI.initialize1()  # for instance, initialize the first case study
t = collect(0:1.:200)     # time points where to output the path mutual information
N = 1000                  # let's do one thousand simulations

# return path-MI estimate
PMI, PMIstd, PMIrate, PMIratestd = PathMI.estimate(S,t,N)
```

Further clarifications can be found by reading the commented code or by contacting the authors at duso@mpi-cbg.de or zechner@mpi-cbg.de . 


