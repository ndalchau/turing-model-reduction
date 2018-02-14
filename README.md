# Model reduction for Turing patterns

This repository contains the Matlab scripts for recreating the figures of the JR Soc Interface article *Model reduction enables Turing instability analysis of large reaction-diffusion models*. 

To reproduce each figure, you should run the relevant script in Matlab:
- Figure1_BrusselatorExample.m
- Figure2_BrusselatorDispersion.m
- Figure3_TuringExample.m
- Figure4_TuringDispersion.m
- Figure5_BorekExample.m
- Figure6_BorekDispersion.m

## Dependencies
- All Matlab code was tested in release 2016a. 
- To create PDFs of figures 2, 4 and 6, you will require the **save2pdf** function from [here](https://uk.mathworks.com/matlabcentral/fileexchange/16179-save2pdf). Additionally, you will need to uncomment the calls to **save2pdf** in each script.
- Mathematica is used to evaluate equilibria of the Borek model, which are stored in CSV files.

## Contributors
- Neil Dalchau, Microsoft Research (repository owner)
- Stephen Smith, University of Edinburgh