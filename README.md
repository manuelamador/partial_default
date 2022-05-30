# Code for "Reputation and Partial Default" 

by Manuel Amador and Christopher Phelan. 


This repository contains code that replicates the numerical results of ["Reputation and Partial Default"](https://manuelamador.me/files/reputationpartial.pdf) by Manuel Amador and Christopher Phelan. 

## Software requirements


- [Julia](https://julialang.org/) version 1.7.2. 
- The file `Project.toml` lists all of the necessary Julia packages.
- The figures are generated using `PGFPlotsX`, a Julia package. Note that a working running `LaTeX` installation is necessary (with the `lualatex` command and the `standalone` package). See [PGFPlotsX](https://kristofferc.github.io/PGFPlotsX.jl/v1/) for details. 
  
## Memory and runtime requirements 

The main script takes less than 10 minutes to run in 4 core Intel machine (i7-1065G7 CPU @ 1.30GHz) with 
16GB of RAM, on Windows 11. 

## Description of code

1. The directory `src` contains the main source code. The file `partial_default.jl` contains the code that runs the simulations. 

2. The directory `scripts` contains the script that generates the figures shown in the paper. 

3. The directory `output` contains the figures generated as an output to the scripts.

## Instructions to replications

Make sure you have installed Julia (a version older than 1.7.2). 

Navigate to the `scripts` subdirectory in the command line and run:  

    > julia paper_figures.jl 

(Alternatively, you can run jupyter, and open and run the `paper_figure.ipynb` notebook.) 

The script should download and precompile the necessary packages the first time it is run (This may take some time).

## Output 

The script runs the simulations and generates all of the figures in the paper. 

## License 

The code is licensed under a MIT license. See `LICENSE.txt` for details.