# Model reduction for Turing patterns

This repository contains the Matlab scripts for recreating the figures of the JR Soc Interface article *Model reduction enables Turing instability analysis of large reaction-diffusion models*. 

Each figure has its own script.
- Figure1_BrusselatorExample
- Figure2_BrusselatorDispersion
- Figure3_TuringExample
- Figure4_TuringDispersion
- Figure5_BorekExample
- Figure6_BorekDispersion

## Dependencies
- All Matlab code was tested in release 2016a. 
- To create PDFs of figures 2, 4 and 6, you will require the **save2pdf** function from [here](https://uk.mathworks.com/matlabcentral/fileexchange/16179-save2pdf). Additionally, you will need to uncomment the calls to **save2pdf** in each script.
- Mathematica is used to evaluate equilibria of the Borek model, which are stored in CSV files.

## Contributors
- Neil Dalchau, Microsoft Research (repository owner)
- Stephen Smith, University of Edinburgh