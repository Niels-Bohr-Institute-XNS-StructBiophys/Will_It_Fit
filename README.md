# Will_It_Fit
******************************
**       Will It Fit        **
**                          **
**     Copyright 2013,      **
** University of Copenhagen **
**                          **
**  Martin Cramer Pedersen  **
**       mcpe@nbi.dk        **
******************************


**************************************************************************
** This file is part of WillItFit.                                      **
**                                                                      **
** WillItFit is free software: you can redistribute it and/or modify    **
** it under the terms of the GNU General Public License as published by **
** the Free Software Foundation, either version 3 of the License, or    **
** (at your option) any later version.                                  **
**                                                                      **
** WillItFit is distributed in the hope that it will be useful,         **
** but WITHOUT ANY WARRANTY; without even the implied warranty of       **
** MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the        **
** GNU General Public License for more details.                         **
**                                                                      **
** You should have received a copy of the GNU General Public License    **
** along with WillItFit. If not, see <http://www.gnu.org/licenses/>.    **
**************************************************************************

If you use WillItFit in your work, please cite:

Pedersen et. al
Journal of Applied Crystallography 46(6), 1894-1898

*********************
* Table of Contents *
*********************

 - Setups
 - Dependencies
 - Fit-a-long - Walkthrough of a refinement
 - Data structure
 - Sample information
 - Adding models
 - Fitting routines
 - Instrumental effects


**********
* Setups *
**********

Currently, WillItFit has been reported working on the following operating systems:

Linux:
 - Ubuntu 10.04 - 13.10
 - CentOS 6.4 

Mac:
 - OSX 10.5 - 10.8

Microsoft:
 - Windows XP
 - Windows Vista
 - Windows 7

Did you make WIllItFit run on a platform not mentioned here? 
Let others know - mcpe@nbi.dk!


****************
* Dependencies *
****************

In order for the interface to function correctly, the following packages must be available:
 - Python 2.7
   http://www.python.org/download/ 

 - WxPython 2.8
   http://www.wxpython.org/download.php
   (Users of Red Hat-based systems can find an .rpm-version of wxPython at http://pkgs.org/).

 - MatPlotLib
   http://matplotlib.org/downloads.html

 - PyLab (which is part of SciPy)
   http://www.scipy.org/Download

Alternatively, Enthought provides an out-of-the-box Python installation with the essential 
packages already installed
   https://www.enthought.com/

In order for the software to compile and run correctly, the following must be in available:
 - A C-compiler - e.g. GCC for Linux/OSX:
   http://gcc.gnu.org/

   or MinGW for Windows:
   http://www.mingw.org/

 - Make sure that the environment/path variables are set up correct - and that the 
   CompilerInfo.cfg-file in the directory is correctly set up.

Some models might require additional support from various libraries - e.g.:
 - The GNU Scientific Library
   http://www.gnu.org/software/gsl/

Optionally, one needs the following library in order to utilize the implemented parallelization:
 - The OpenMP-extension for C
   http://www.openmp.org/


********************************************
* Fit-a-long - Walkthrough of a refinement *
********************************************

The following is an example of how to perform a fit to the data located in the Example-folder.

 - Install the dependencies described in the previous chapter.
   (If OpenMP is not available, the -fopenmp flag should not be raised).


 - Open the Python interface - e.g. via a shell using the command:

   python WillItFit.py   


 - Assign, which data-file should be fitted by clicking the top-most Browse-button and navigate the Examples/Data-folder and select the .card-file.


 - Assign, which sample information-file should be used by clicking the Browse-button below the previous and select the .dat-file in the same directory.


 - Assign, which model should be used by selecting the "Nanodisc with Tags"-model from the dropdown-menu. This causes the list of parameters describing the nanodisc to appear in the right side of the interface. This should enable the button labelled "Execute" in the bottom of the interface.


 - Try computing the scattering profile from the default parameters by clicking the Execute button. Once the C-part of the program finishes, it returns the chisquare. Try clicking the Plot-button to see the data and the scattering profile from the default parameters.


 - Try to fit the backgrounds of the datasets - by checking the parameters BackN100, BackN0 and BackX off in the column labelled "Fit?". Try to use the Levenberg-Marquardt-algorithm for this. Press the Execute button. When the interface regains control, press the plot button again - hopefully, the high-q region of the fit should look slightly more reasonable now.


 - Try to set more parameters free. In order to achieve the fit presented in:

   Pedersen et. al
   Journal of Applied Crystallography 46(6), 1894-1898

   the parameters AxisRatio, AreaPerHeadgroup, NumberOfLipids, RGOfTag, XRoughness, CVOfBelts, CVOfLipids and NNoughness were refined along with the backgrounds. The meaning of these parameters are all described in:

   Skar-Gislinge et. al
   Journal of American Chemical Society 132(39), 13713-22.

   Rerun one of the fitting algorithms in order to improve the refinement of these parameters. Be aware that some of the algorithms are quite time-consuming! Also, Profile Likelihood is not an optimization algorithm, so it should not be used before a reasonable minimum has been reached.


 - Upon (correct) conclusion of the program, all progress will be saved in a folder in the directory of the .card-file. Should one with to continue from a previously obtained set of parameters, one can load the LastFit.par-file in this directory via the button labelled "Load parameter file".


******************
* Data structure *
******************

The data-structure used by this software is perhaps explained the best by the example found in the Examples-folder.

In short, a .card-file should contain a list of all datasets to be included in the fitting process. Each dataset must be denoted by a concentration (which should be supplied in mM - upon read-in, this quantity is converted to 1/cm^3), a contrast (0 - 100 for SANS-data representing the D2O-content in the solvent - any other value is considered an X-ray-dataset) as well as a scaling and a background, which will be applied to the data on input.

Consider rebinning your data with e.g. the software available at Sourceforge.net/Projects/WillItRebin. The fewer points in a dataset, the faster the model computations will be executed.

Please be very aware of whitespaces, special characters and escape characters in the filenames (in e.g. your .card-file). Depending on your operating system, these may be treated in very different manners. Be particularly aware of the directory delimiter (backslash - \) in Windows - as this character is also a meta-character in C, this character should be escaped using another backslash on Windows OS'es.


**********************
* Sample information *
**********************

In order to use the implemented models, the program must be fed information on the scattering and volume of the various components. These informations are contained in a single file, examples of which are kept in the Examples-folder - one matching each model is available.
The format is rather strict and should not be altered, as this will lead to errors in the I/O-routines.

PDB-files are read in as in the associated structure named "Protein" - for details on this struct, see Structs.h in the Auxillary-folder. The protein is read as atoms as well as residues for which the center of scattering is computed as well as the center of volume.

Likewise, the declaration of the structure "UserDefinedStructure" can be found in Structs.h The intention is that user can put whichever variables, they find necessary, in the struct, and the software will transport them to the model computations.


*****************
* Adding models *
*****************

In order to add a model to the list of available models, one should create a new folder in the Models-directory. A ModelInfo.h-file in the correct format must be located in this directory - and this file must #include sufficient headers files so that the functions ComputeConstraints, Model and OutputData are defined at compile-time.

Please consult the available example located in the Examples-directory in order to get the syntax right.

Correspondingly, if one wishes to remove a model from the list of available models in the interface, one simply needs to remove the associated folder from the Models-directory (and restart the software).


********************
* Fitting routines *
********************

The software features the following fitting routines:
 - Levenberg-Marquardt
 - Broyden-Fletcher-Goldfarb-Shanno
 - Particle Swarm Optimizer
 - Genetic Algorithm

The appropriate references for all of these algorithms can be found in:

Pedersen et. al
Journal of Applied Crystallography 46(6), 1894-1898

Generally speaking, the Swarm algorithm and the Genetic algorithm will outperform the two Quasi-Newtonian algorithms in complicated optimization landscapes, whereas the LM- and the BFGS-algorithm executes considerably faster and is better suited for searches in the vicinity of a minimum.


************************
* Instrumental effects *
************************

In order to implement the effects of instrumental resolution in a fit, one must supply a .res-file (similar to the one found in the Examples-folder). 

This information is used as outlined in:

Pedersen et al.
Journal of Applied Crystallography, 23(4), 321-333.

to smear the fit.

This is a very time-comsuming process, as the model is computed N times per value of q, where N is the number of requested smearing folds.


