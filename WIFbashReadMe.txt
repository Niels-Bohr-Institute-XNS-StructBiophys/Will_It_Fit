**Introduction and functionality**
 
Welcome to WillItFit(WIF) bash edition
 
This branch of WIF offers to following functionality:
 
Fits can only be performed on one model and with only one fitting routine. The default model is that of a nanodisc and the default fitting routine is the Levenberg-Marquardt algorithm.
 
This reduction in functionality reduces the amount of arguments that is passed when calling the c part of WIF, which makes it easier to call in bash mode.
 
**Calling WIF**
 
All call to WIF takes the following arguments:
 
	-c path to the datafiles - .card
	-s path to the samples - .dat
	-p path to the parameter file -.mcp
	-n minimum value for q -default 0.0
	-x maximum value for q -default 1.0
	-t Number of steps for fitting routine (When using resolution files)
	-e path to Resolution file
	-d path to PDB file
 
To call WIF in bash mode, the first must be compiled first. The following is an example a command used to compile WillItFit.c:
 
 
The compiled file can now be called:
 
./wifcmd.com -c=cardfile.card -s=datfile.dat -p=parameterfile.mcp e=N/A -d=N/A 
 
**Model and parameter selection**
 
A different model can be chosen by altering the path given in Auxillary/ModelLocation.h. This source code needs to be recompiled for this change to take effect.
 
The free parameters used for fitting can be changed in the .mcp file passed to the compiled program. A parameter is free if the number in the fourth column is zero. Examples of .mcp files for different models can be found in the Models folder.
 
**Output**
 
The program outputs results in the following files and directories:
 
Results.wif
 
This file contains information about the files used when fitting, the final parameter values and the quality of the fit.
 
cardfile.card-results
 
This directory contains the following files:
 
ParametersAfter.mcp
 
This file essentially contains the same information as the .mcp file used when running the program.
 
dat.dat files
 
These files contain information about the individual data points in the spectra used when fitting. One is output for each spectra.
 
fits.dat files
 
These files contain information about the individual data points in the spectra obtained when fitting. One is output for each spectra.
 
 
 
 
 
 
 

