# PathMI
Julia code for calculating mutual information between complete trajectories of molecular components.
For further information, contact Lorenzo Duso (duso@mpi-cbg.de) & Christoph Zechner (zechner@mpi-cbg.de)


# INTRODUCTION

All the numerical simulations in "Path mutual information for a class of biochemical reaction networks" ( http://arxiv.org/abs/1904.01988 ) have been performed with the programming language julia ( http://julialang.org/ ). Here we want to share our code for reproducibility reasons. The provided code is customized to reproduce the case studies shown in the paper and is not meant to be a general tool.

Here we give instructions on how to reproduce the plots of the paper by using the module PathMI.jl and the prepared routines in the associated launch.jl modules. Please note that the launch.jl modules need PathMI.jl to be saved in the same folder.


# INSTALLATION AND REQUIREMENTS

The provided code is compatible with all julia versions v0.6 and v1. In particular:
- if you are using julia v0.6, you will need to include launch_ver06.jl in your julia-terminal
- if you are using julia v1., you will need to include launch_ver1.jl in your julia-terminal 

The module PathMI.jl itself is compatible with all the versions and does not require the installation of any package to be used. Instead, the launch.jl modules require the package PyPlot in order to plot the results.

Adding the package PyPlot.
In case you have just installed julia or you do not have the package PyPlot installed, please follow these intstructions:

- for version v0.6, open the julia terminal and paste
```
Pkg.add("PyPlot")
using PyPlot
```

- for version v1., open the julia terminal and paste
```
import Pkg
Pkg.add("PyPlot")
using PyPlot
```
 
Note: on macOS, in case the plotting of some outputs in launch.jl fails or is too slow, it may be helpful to run the following commands and include again the launch.jl module:
```
using Conda
Conda.add("pyqt")
using PyPlot
```


# USAGE

The launch.jl modules simply interface with the module PathMI.jl, where the actual simulation code is contained, and are responsible to plot the outputs. The functions in launch.jl are named accordingly to the figure of the paper they reproduce. For instance, launch.figure1b() will run the code to output figure 1b. We give here below an example.

Example of usage. Reproducing figure 1b in julia version v1. :  
1) open the julia terminal and access the folder where PathMI.jl and launch_ver1.jl are located
2) load the module launch_ver1.jl by running
```
include("launch_ver1.jl")
```
3) run the function to reproduce figure 1b
```
launch.figure1b()
```

If it is desired to get precisely the identical figures as reported in the paper, it is possible to pass the argument ``` paper=true ```, for example ```launch.figure1b(paper=true)``` . This command sets the seed of the random number generator to the same value used when making the figures.

# ADVANCED USAGE

It is also possible to directly use the module PathMI.jl by running in the julia terminal :
```
include("PathMI.jl")
```
This allows you to manually perform simulations and change system settings such as the reaction network, the rate constants, the initial condition and the simulation time. Further clarifications can be found by reading the commented code and usage examples in PathMI.jl or by contacting the authors at duso@mpi-cbg.de or zechner@mpi-cbg.de . 


