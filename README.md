# Numerical analysis for "Reputation and Partial Default" 

by Manuel Amador and Christopher Phelan. 

This repository contains code that replicates the numerical results of ["Reputation and Partial Default"](https://manuelamador.me/files/reputationpartial.pdf) by Manuel Amador and Christopher Phelan. 

## Description

The code is in Julia. The code was run using [Julia](https://julialang.org/) 1.7.2. 


1. The directory `src` contains the main source code. 

2. The directory `scripts` contains the script that generates the figures shown in the paper. 

3. The directory `output` contains the figures.

## Running to the code

Make sure you have installed Julia (a version older than 1.7.2). 

Navigate to the `scripts` subdirectory in the command line and run:  

    > julia paper_figures.jl 

(Alternatively, you can run jupyter, and open and run the `paper_figure.ipynb` notebook.) 

The script will run the simulations and generate the figures in the output directory. 

The script should download and precompile the necessary packages the first time it is run (This may take some time). 