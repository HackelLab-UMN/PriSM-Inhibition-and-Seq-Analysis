# PriSM-Inhibition-and-Seq-Analysis

Created by: Abbigael Harthorn, vandu054@umn.edu
Corresponding author: Benjamin Hackel, hackel@umn.edu 
Hackel Lab, University of Minnesota - Twin Cities, 2021

This project contains Python3 scripts for deep sequencing analysis of protein-small molecule conjugates, as well as MATLAB scripts for determination of apparent Ki. 

For deep sequencing analysis, sequences for each of the campaigns are found as pickle files and found in the seq_files folder. Sequence files are labeled by cysteine site and CA target. These files must be in the same directory as the analysis codes to run properly. Included python scripts were used to analyze sequence data and generate figures. In detail, scripts analyze Hamming distance, relative sequence frequency, and sequence enrichment compared to naive library.

For determination of apparent Ki, there is one MATLAB script to fit the Morrison Equation for tight binding inhibitors. Absorbance data is saved in an excel file and is imported into MATLAB.
