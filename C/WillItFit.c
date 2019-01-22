///       Will It Fit        ///
///                          ///
///     Copyright 2013,      ///
/// University of Copenhagen ///
///                          ///
///  Martin Cramer Pedersen  ///
///       mcpe@nbi.dk        ///

/// This file is part of WillItFit.                                      ///
///                                                                      ///
/// WillItFit is free software: you can redistribute it and/or modify    ///
/// it under the terms of the GNU General Public License as published by ///
/// the Free Software Foundation, either version 3 of the License, or    ///
/// (at your option) any later version.                                  ///
///                                                                      ///
/// WillItFit is distributed in the hope that it will be useful,         ///
/// but WITHOUT ANY WARRANTY; without even the implied warranty of       ///
/// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the        ///
/// GNU General Public License for more details.                         ///
///                                                                      ///
/// You should have received a copy of the GNU General Public License    ///
/// along with WillItFit. If not, see <http://www.gnu.org/licenses/>.    ///

/// If you use this software in your work, please cite: ///
///                                                     ///
/// Pedersen, M. C., Arleth, L. & Mortensen, K.         ///
/// J. Appl. Cryst. 46(6), 1894-1898                    ///

/// Defines
#define pi 3.14159265
#define EstimatedNumberOfStepsInLikelihoodProfile 10
#define MaxNumberOfConstraints 50
#define MaxSizeOfCardfile 1024

/// Inclusions
// Include build-in libraries
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include <complex.h>
#include <time.h>
#include <string.h>


// create directory under Windows
#if defined(_WIN32)
#include <direct.h>
#endif


// Include parallelisation
#include "Auxillary/Parallelisation.h"

// Include supporting headers
#include "Auxillary/Structs.h"
#include "Auxillary/Allocation.h"
#include "Auxillary/AuxillaryFunctions.h"
#include "Auxillary/IncludeUserDefined.h"

// Include mathematical headers
/*
#include "Formfactors/FormfactorDebye.h"
#include "Formfactors/FormfactorHalfBelt.h"
#include "Formfactors/FormfactorSphere.h"
#include "Formfactors/FormfactorCylinder.h"
#include "Formfactors/FormfactorCylinderDisplaced.h"
#include "Formfactors/FormfactorHammouda.h"
#include "Formfactors/FormfactorEllipticCylinderWithEllipsoidalEndcaps.h"
#include "Formfactors/FormfactorRodWithoutEndcaps.h"
*/

// Include model
#include "Auxillary/ModelLocation.h"

// Include I/O
#include "InputOutput/ImportSpectra.h"
#include "InputOutput/ImportParameters.h"
#include "InputOutput/ImportSampleInfo.h"
#include "InputOutput/ImportPDBFile.h"
#include "InputOutput/OutputSpectra.h"
#include "InputOutput/OutputParameters.h"

// Include functions accounting for instrumental smearing
#include "Resolution/SigmaOfQ.h"
#include "Resolution/Resolution.h"

// Import fitting routines
#include "FittingRoutines/LevenbergMarquardtSupportingFunctions.h"
#include "FittingRoutines/BFGSSupportingFunctions.h"
#include "FittingRoutines/SwarmSupportingFunctions.h"
#include "FittingRoutines/GeneticSupportingFunctions.h"
#include "FittingRoutines/LevenbergMarquardt.h"
#include "FittingRoutines/BFGS.h"
#include "FittingRoutines/ComputeModel.h"
#include "FittingRoutines/Swarm.h"
#include "FittingRoutines/Genetic.h"
#include "FittingRoutines/GridsearchBFGS.h"
#include "FittingRoutines/GridsearchLM.h"
#include "FittingRoutines/ProfileLikelihood.h"
#include "FittingRoutines/ProfileLikelihoodSingleParameter.h"

/// Main program
int main(int argc, char *argv[])
{
	/// Declarations

	// Variables describing the spectras
	struct Dataset * Data;
	char CardFileLocation[256] ;

	int NumberOfSpectra;
	int HighestNumberOfDatapoints;
	int TotalNumberOfDatapoints = 0;

	// Variables describing the sample info
	char SamplesFileLocation[256] ;
	int NumberOfSampleInformations;

	char ResultsDirectory[256] ;
	char LogFileLocation[256] ;

	double *VolumesOfMolecules;
	char ModificationName[4] ;

	char PDBFileLocation[256] ;
	struct Protein ProteinStructure;
	struct UserDefined UserDefinedStructure;

	// Variables describing the parameters
	struct Parameter * Parameters;
	char ParameterFileLocation[256] ;
	int NumberOfParameters = 0;
	int NumberOfFreeParameters = 0;

	// Variables describing the properties of the fit
	double QMin = 0.0;
	double QMax = 1.0;
	double DeltaForDifferentiations = 0.001;
	double ChiSquare;
	double ChiSquareFractile;

	int ChooseFittingRoutine = 0;
	int FittingRoutineArgument1 = 50;
	int FittingRoutineArgument2 = 10;
	int FittingRoutineArgument3 = 32;
	int FittingRoutineError = -1;

	bool PrintCovarianceMatrix = false;
	char Message[256] ;

	// Variables describing the resolution
	char ResolutionFileLocation[256] ;
	int NumberOfSmearingFolds = 0;
	bool IncludeResolutionEffects = false;

	// CMD mode or py
	bool CMD = true ;

	// how to write a log(file), by default write most important output to stdout (terminal)
	int WriteLog = -1 ;
	FILE *logfile ;
	int LOWER_MAX_INDEX = 0 ; 
	int UPPER_MIN_INDEX = 0 ;

	/// Obtain arguments from program or request them in console
	ErrorCheck( AssignArguments(argc, argv, CardFileLocation, SamplesFileLocation, ParameterFileLocation, PDBFileLocation, &QMin, &QMax, &ChooseFittingRoutine, &FittingRoutineArgument1, &FittingRoutineArgument2, &FittingRoutineArgument3, &IncludeResolutionEffects, &NumberOfSmearingFolds, ResolutionFileLocation, &PrintCovarianceMatrix, &ChiSquareFractile, &CMD, &WriteLog), "when assigning the arguments from the command line in AssignArguments()", NULL, -1, stdout) ;


	// Define & create directory for output
	// This is platform specific, on Mac it works without, on Linux and maybe Windows it has to be done explicitely
	//
	// ResultsDirectory includes final "/"
	sprintf( ResultsDirectory, "%s-results/", CardFileLocation) ; // sprintf appends automatically 0

	#if defined(_WIN32)
	_mkdir( ResultsDirectory ) ;
	#elif defined(__linux__)
	mkdir( ResultsDirectory, 0700) ;
	#else
	// for Mac do nothing
	#endif


	/*
		logging

		WriteLog  < 0 logfile points to stdout (terminal)
		WriteLog == 0 logfile is NULL
		WriteLog  > 0 logfile points to a log-file in the generated results under the .card-file folder

		abs(WriteLog) == 0 -> no output at all

		abs(WriteLog) == 1 -> most important output written to file pointer logfile, includes:
					-parameter file
					-excerpt of spectra, ProteinStructure
					-ChiSquare and parameters at the fitting iteration steps

		abs(WriteLog) == 2 -> full output written to file pointer logfile, includes in addition to 1:
					full output of Spectra, ProteinStructure
					output from fitting algorithm (currently Levenberg-Marquardt and Compute Model supported)

		default is -1
	*/
	if ( WriteLog > 0 )
	{
		sprintf( LogFileLocation, "%s%s", ResultsDirectory, "logfile.log") ; // sprintf appends automatically 0
		printf( "\n") ;
		printf( "Write logfile at %s\n", LogFileLocation) ;
		printf( "\n") ;
		fflush( stdout ) ;

		logfile = fopen( LogFileLocation, "w") ;
	}
	else if ( WriteLog < 0 )
	{
		logfile = stdout ;
	}
	else
	{
		printf( "Warning: No log(file) will be printed or written!\n") ;
		fflush( stdout ) ;

		logfile = NULL ;
	}


	if ( abs(WriteLog) > 0 )
	{
		ClearScreen( logfile ) ;
		fprintf( logfile, "                                 Will It Fit                                    \n\n") ;
		fprintf( logfile, "                   Fitting Macromolecules with Modifications                    \n\n") ;
		fprintf( logfile, "                            Version 0.02 2018-04-12                             \n\n") ;
		fprintf( logfile, "\n") ;
		fprintf( logfile, "Based on the Will It Fit software (M.C. Pedersen, L. Arleth and K. Mortensen, \"WillItFit: A framework for fitting of constrained models to small-angle scattering data\", J. Appl. Cryst. 46, 1894-1898 (2013) doi:10.1107/S0021889813026022)\n") ;
		fprintf( logfile, "\n") ;
		fprintf( logfile, "Code adaptations implemented by Nicholas Skar-Gislinge, Asger Neesgaard Sand and Martin Schmiele\n") ;
		fprintf( logfile, "\n") ;
		fprintf( logfile, "LINX contacts: Erik Brok and Martin Schmiele\n") ;
		fprintf( logfile, "\n") ;
		ClearScreen( logfile ) ;
		fprintf( logfile, "\n\n") ;
		ClearScreen( logfile ) ;
		fprintf( logfile, "\n") ;

		fflush( logfile ) ;
	}

	if ( abs(WriteLog) > 0 )
	{
		fprintf ( logfile, "Program called with: ") ;
		for ( int i=0; i<argc; ++i) { fprintf ( logfile, "%s ", argv[i]) ; }
		fprintf( logfile, "\n") ;
		fprintf( logfile, "\n\n") ;
		fprintf( logfile, "Summarizing assigned arguments:\n") ;
		fprintf( logfile, "\n") ;
		fprintf( logfile, "\tCMD = %d\n", (int)CMD) ;
		fprintf( logfile, "\n") ;
		fprintf( logfile, "\tCardFileLocation = %s\n", CardFileLocation) ;
		fprintf( logfile, "\tSamplesFileLocation = %s\n", SamplesFileLocation) ;
		fprintf( logfile, "\tParameterFileLocation = %s\n", ParameterFileLocation) ;
		fprintf( logfile, "\tPDBFileLocation = %s\n", PDBFileLocation) ;
		fprintf( logfile, "\n") ;
		fprintf( logfile, "\tResolutionFileLocation = %s\n", ResolutionFileLocation) ;
		fprintf( logfile, "\tIncludeResolutionEffects = %d\n", (int)IncludeResolutionEffects) ;
		fprintf( logfile, "\tNumberOfSmearingFolds = %d\n", NumberOfSmearingFolds) ;
		fprintf( logfile, "\n") ;
		fprintf( logfile, "\tQ = %lf - %lf [1/A]\n", QMin, QMax) ;
		fprintf( logfile, "\n") ;
		fprintf( logfile, "\tFit routine number = %d\n", ChooseFittingRoutine) ;
		fprintf( logfile, "\tFittingRoutineArgument1 = %d\n", FittingRoutineArgument1) ;
		fprintf( logfile, "\tFittingRoutineArgument2 = %d\n", FittingRoutineArgument2) ;
		fprintf( logfile, "\tFittingRoutineArgument3 = %d\n", FittingRoutineArgument3) ;
		fprintf( logfile, "\n") ;
		fprintf( logfile, "\tPrintCovarianceMatrix = %d\n", (int)PrintCovarianceMatrix) ;
		fprintf( logfile, "\tChiSquareFractile = %lf\n", ChiSquareFractile) ;
		fprintf( logfile, "\n\n") ;

		fflush( logfile ) ;
	}



	/// Retrieve parameters from .par-file
	if ( abs(WriteLog) > 0 )
	{
		ClearScreen( logfile ) ;
		fprintf( logfile, "\n") ;
		fprintf( logfile, "Read parameters:\n") ;
		fprintf( logfile, "\n") ;
		fprintf( logfile, "\tRead parameters from %s\n", ParameterFileLocation) ;
		fprintf( logfile, "\n") ;

		fflush( logfile ) ;
	}

	NumberOfParameters = CheckNumberOfParameters( ParameterFileLocation ) ;
	ErrorCheck(NumberOfParameters, "when reading the parameter file in CheckNumberOfParameters()", ResultsDirectory, WriteLog, logfile) ;

	AllocateParameters(&Parameters, NumberOfParameters) ;
	ImportParameters(Parameters, ParameterFileLocation) ;

	for ( int i = 0; i < NumberOfParameters; ++i)
	{
		if (Parameters[i].iParameter == true) { NumberOfFreeParameters += 1 ; }
	}

	if ( abs(WriteLog) > 0 )
	{
		fprintf( logfile, "\tFound %d initial fit parameters, thereof %d are free parameters\n", NumberOfParameters, NumberOfFreeParameters) ;
		fprintf( logfile, "\n") ;
		fprintf( logfile, "\tParameters[i]\n") ;
		fprintf( logfile, "\n") ;
		fprintf( logfile, "\t\ti\n") ;
		fprintf( logfile, "\t\t.Name\n") ;
		fprintf( logfile, "\t\t.MinValue\n") ;
		fprintf( logfile, "\t\t.Value\n") ;
		fprintf( logfile, "\t\t.MaxValue\n") ;
		fprintf( logfile, "\t\t.(int)iParameter\n") ;
		fprintf( logfile, "\t\t.Error\n") ;
		fprintf( logfile, "\n") ;
		for ( int i = 0; i < NumberOfParameters; ++i)
		{
			fprintf( logfile, "\t\t%2d %10s %-15g %-15g %-15g %d %-15g\n", i, Parameters[i].Name, Parameters[i].MinValue, Parameters[i].Value, Parameters[i].MaxValue, (int)Parameters[i].iParameter, Parameters[i].Error) ;
		}
		fprintf( logfile, "\n\n") ;

		fflush( logfile ) ;
	}

	if ( NumberOfFreeParameters == 0 && ChooseFittingRoutine > 0 )
	{
		ErrorCheck( -1, "since there are no free parameters but a fitting algorithm was chosen", ResultsDirectory, WriteLog, logfile) ;
	}



	/// Import data from .card-file and read spectra
	if ( abs(WriteLog) > 0 )
	{
		ClearScreen( logfile ) ;
		fprintf( logfile, "\n") ;
		fprintf( logfile, "Read .card-file and import the spectra:\n") ;
		fprintf( logfile, "\n") ;
		fprintf( logfile, "\tRead HighestNumberOfDatapoints and NumberOfSpectra via CheckSizeOfData() from %s\n", CardFileLocation) ;
		fprintf( logfile, "\n") ;

		fflush( logfile ) ;
	}

	HighestNumberOfDatapoints = CheckSizeOfData( CardFileLocation, &NumberOfSpectra, CMD, WriteLog, logfile) ;
	ErrorCheck( HighestNumberOfDatapoints, "when reading the .card-file or the data files in CheckSizeOfData()", ResultsDirectory, WriteLog, logfile) ;

	AllocateData(&Data, NumberOfSpectra) ;

	for ( int i = 0; i < NumberOfSpectra; ++i)
	{
		Initialize1DArray(&Data[i].QValues,           HighestNumberOfDatapoints) ;
		Initialize1DArray(&Data[i].IValues,           HighestNumberOfDatapoints) ;
		Initialize1DArray(&Data[i].FitValues,         HighestNumberOfDatapoints) ;
		Initialize1DArray(&Data[i].SigmaValues,       HighestNumberOfDatapoints) ;
		Initialize1DArray(&Data[i].SigmaQValues,      HighestNumberOfDatapoints) ;
		Initialize2DArray(&Data[i].ResolutionWeights, HighestNumberOfDatapoints, NumberOfSmearingFolds) ;
		Initialize1DArray(&Data[i].Constraints,       MaxNumberOfConstraints) ;

		Data[i].IncludeResolutionEffects = false;
	}

	if ( abs(WriteLog) > 0 )
	{
		fprintf( logfile, "\n") ;
		fprintf( logfile, "\tImporting spectra via ImportSpectra() from %s\n", CardFileLocation) ;
		fprintf( logfile, "\n") ;

		fflush( logfile ) ;
	}
	ErrorCheck( ImportSpectra( Data, CardFileLocation, NumberOfSpectra, CMD, WriteLog, logfile) , "when reading the data files in ImportSpectra()", ResultsDirectory, WriteLog, logfile) ;

	for ( int i = 0; i < NumberOfSpectra; ++i) { TotalNumberOfDatapoints += Data[i].NumberOfDatapoints ; }

	if ( abs(WriteLog) > 0 ) { fprintf( logfile, "\n") ; }



	/// Decide fitting range
	if ( abs(WriteLog) > 0 )
	{
		ClearScreen( logfile ) ;
		fprintf( logfile, "\n") ;
		fprintf( logfile, "Assign fitting ranges in spectra:\n") ;
		fprintf( logfile, "\n") ;
		fprintf( logfile, "\tAssignFittingRanges()\n") ;
		fprintf( logfile, "\n") ;

		fflush( logfile ) ;
	}
	AssignFittingRanges( Data, QMin, QMax, NumberOfSpectra, WriteLog, logfile) ;




	/// Include or exclude resolution effects
	if ( abs(WriteLog) > 0 )
	{
		ClearScreen( logfile ) ;
		fprintf( logfile, "\n") ;
		fprintf( logfile, "Resolution effects in Q:\n") ;
		fprintf( logfile, "\n") ;

		fflush( logfile ) ;
	}

	if ( IncludeResolutionEffects == true )
	{
		if ( abs(WriteLog) > 0 )
		{
			fprintf( logfile, "\tIncluding resolution effects via Resolution() from %s\n", ResolutionFileLocation) ;
			fprintf( logfile, "\n") ;

			fflush( logfile ) ;
		}
		ErrorCheck( Resolution(Data, NumberOfSpectra, ResolutionFileLocation, NumberOfSmearingFolds, WriteLog, logfile), "when reading the resolution information-file in Resolution()", ResultsDirectory, WriteLog, logfile) ;

		if ( abs(WriteLog) > 0 ) { fprintf( logfile, "\n") ; }
	}
	else
	{
		if ( abs(WriteLog) > 0 )
		{
			fprintf( logfile, "\tExcluding resolution effects\n") ;
			fprintf( logfile, "\n\n") ;

			fflush( logfile ) ;
		}
	}



	/// Print Data
	// skip ResolutionWeights, Constraints and FitValues
	// ScatteringLengths is printed later after assigning them from sample information file
	if ( abs(WriteLog) > 0 )
	{
		/*
		struct Dataset {
		    double * QValues;
		    double * IValues;
		    double * SigmaValues;
		    double * SigmaQValues;
		    double ** ResolutionWeights;
		    double * FitValues;
		    double Concentration;
		    double Contrast;
		    double * Constraints;
		    double * ScatteringLengths;
		    bool IncludeResolutionEffects;
		    int NMin;
		    int NMax;
		    int NumberOfDatapoints;
		    char Filename[256];
		};
		*/
		ClearScreen( logfile ) ;
		fprintf( logfile, "\n") ;
		fprintf( logfile, "Print Data information:\n") ;
		fprintf( logfile, "\n") ;


		LOWER_MAX_INDEX = 0 ; 
		UPPER_MIN_INDEX = 0 ;

		fprintf( logfile, "\tNumberOfSpectra = %d\n", NumberOfSpectra) ;
		fprintf( logfile, "\tHighestNumberOfDatapoints = %d\n", HighestNumberOfDatapoints) ;

		for ( int i = 0; i < NumberOfSpectra; ++i)
		{
			fprintf( logfile, "\n") ;
			fprintf( logfile, "\tSpectra no. %d\n", i+1) ;
			fprintf( logfile, "\n") ;
			fprintf( logfile, "\t\tData[%d]\n", i) ;
			fprintf( logfile, "\t\t\t.Filename                 = %s\n", Data[i].Filename) ;
			fprintf( logfile, "\t\t\t.NumberOfDatapoints       = %d\n", Data[i].NumberOfDatapoints) ;
			fprintf( logfile, "\t\t\t.NMin                     = %d\n", Data[i].NMin) ;
			fprintf( logfile, "\t\t\t.NMax                     = %d\n", Data[i].NMax) ;
			fprintf( logfile, "\t\t\t.IncludeResolutionEffects = %d\n", (int)Data[i].IncludeResolutionEffects) ;
			fprintf( logfile, "\t\t\t.Concentration            = %lf\n", Data[i].Concentration) ;
			fprintf( logfile, "\t\t\t.Contrast                 = %lf\n", Data[i].Contrast) ;
			fprintf( logfile, "\n") ;
			fprintf( logfile, "\t\t\tj\n") ;
			fprintf( logfile, "\t\t\t.QValues[j]\n") ;
			fprintf( logfile, "\t\t\t.IValues[j]\n") ;
			fprintf( logfile, "\t\t\t.SigmaValues[j]\n") ;
			fprintf( logfile, "\t\t\t.SigmaQValues[j]\n") ;

			if ( abs(WriteLog) < 2 ) 
			{
				LOWER_MAX_INDEX = min( 3, Data[i].NumberOfDatapoints) ; 
				UPPER_MIN_INDEX = max( Data[i].NumberOfDatapoints - LOWER_MAX_INDEX, LOWER_MAX_INDEX) ;
			}
			for ( int j = 0; j < max( 0, LOWER_MAX_INDEX); ++j)
			{
				fprintf( logfile, "\t\t\t%3d %-15f %-15f %-15f %-15f\n", j, Data[i].QValues[j], Data[i].IValues[j], Data[i].SigmaValues[j], Data[i].SigmaQValues[j]) ;
			}
			if ( abs(WriteLog) < 2 ) { fprintf( logfile, "\n\t\t\t... (for full output use -l=2 or -l=-2 option) ...\n\n") ; }
			for ( int j = max( 0, UPPER_MIN_INDEX); j < Data[i].NumberOfDatapoints; ++j)
			{
				fprintf( logfile, "\t\t\t%3d %-15f %-15f %-15f %-15f\n", j, Data[i].QValues[j], Data[i].IValues[j], Data[i].SigmaValues[j], Data[i].SigmaQValues[j]) ;
			}

		}
		fprintf( logfile, "\n\n") ;

		fflush( logfile ) ;
	}



	/// Retrieve sample informations 
	if ( abs(WriteLog) > 0 )
	{
		ClearScreen( logfile ) ;
		fprintf( logfile, "\n") ;
		fprintf( logfile, "Read sample informations:\n") ;
		fprintf( logfile, "\n") ;
		fprintf( logfile, "\tRead NumberOfSampleInformations via CheckSizeOfSampleInformation() from %s\n", SamplesFileLocation) ;

		fflush( logfile ) ;
	}

	// get NumberOfSampleInformations
	NumberOfSampleInformations = CheckSizeOfSampleInformation( SamplesFileLocation ) ;
	ErrorCheck( NumberOfSampleInformations, "when reading the sample information file in CheckSizeOfSampleInformation()", ResultsDirectory, WriteLog, logfile) ;

	if ( abs(WriteLog) > 0 )
	{
		fprintf( logfile, "\tFound %d sample informations\n", NumberOfSampleInformations) ;
		fprintf( logfile, "\n") ;

		fflush( logfile ) ;
	}

	// initializations and import sample informations
	// e.g. volume of solvent molecules, SLD of solvent molecules (in case of neutron scattering weighted with the provided sample contrast)
	Initialize1DArray( &VolumesOfMolecules, NumberOfSampleInformations) ;
	for ( int i = 0; i < NumberOfSpectra; ++i) { Initialize1DArray( &Data[i].ScatteringLengths, NumberOfSampleInformations) ; }

	if ( abs(WriteLog) > 0 )
	{
		fprintf( logfile, "\tRead sample informations via ImportSampleInformation() from %s\n", SamplesFileLocation) ;
		fprintf( logfile, "\t\tVolumes of molecules will be stored in array VolumesOfMolecules[]\n") ;
		fprintf( logfile, "\t\tScattering lengths will be assigned to corresponding dataset in Data[].ScatteringLengths (in case of neutron scattering weighted with the provided sample contrast)\n") ;
		fprintf( logfile, "\n") ;

		fflush( logfile ) ;
	}
	ErrorCheck( ImportSampleInformation( Data, VolumesOfMolecules, SamplesFileLocation, NumberOfSampleInformations, NumberOfSpectra, ModificationName), "when reading the sample informationfile in ImportSampleInformation()", ResultsDirectory, WriteLog, logfile) ;

	if ( abs(WriteLog) > 0 )
	{
		fprintf( logfile, "\t\ti\n") ;
		fprintf( logfile, "\t\tVolumesOfMolecules[i]\n") ;
		fprintf( logfile, "\n") ;
		for ( int i = 0; i < NumberOfSampleInformations; ++i)
		{
			fprintf( logfile, "\t\t%2d %-15f\n", i, VolumesOfMolecules[i]) ;
		}
		fprintf( logfile, "\n\n") ;
		fprintf( logfile, "\t\ti\n") ;
		fprintf( logfile, "\t\tj\n") ;
		fprintf( logfile, "\t\tData[i].ScatteringLengths[j]\n") ;
		fprintf( logfile, "\n") ;
		for ( int i = 0; i < NumberOfSpectra; ++i)
		{
			for ( int j = 0; j < NumberOfSampleInformations; ++j)
			{
				fprintf( logfile, "\t\t%2d %2d %-15g\n", i, j, Data[i].ScatteringLengths[j]) ;
			}
		}
		fprintf( logfile, "\n") ;
		fprintf( logfile, "\n") ;
		fprintf( logfile, "\t\tModification label is '%s'\n", ModificationName) ;
		fprintf( logfile, "\n\n") ;

		fflush( logfile ) ;
	}



	/// Import the PDB-file

	// allocate ProteinStructure by reading ProteinStructure.NumberOfResidues, ProteinStructure.NumberOfAtoms
	ProteinStructure.NumberOfAtoms = 0 ;

	if (strcmp(PDBFileLocation, "N/A") != 0)
	{
		if ( abs(WriteLog) > 0 )
		{
			ClearScreen( logfile ) ;
			fprintf( logfile, "\n") ;
			fprintf( logfile, "Importing PDB-file:\n") ;
			fprintf( logfile, "\n") ;

			fflush( logfile ) ;
		}

		sprintf(ProteinStructure.PDBFileLocation, "%s", PDBFileLocation) ;

		if ( abs(WriteLog) > 0 )
		{
			fprintf( logfile, "\tRead ProteinStructure.NumberOfResidues via CheckNumberOfResiduesInPDBFile() from %s\n", PDBFileLocation) ;
			fprintf( logfile, "\n") ;

			fflush( logfile ) ;
		}
		ProteinStructure.NumberOfResidues = CheckNumberOfResiduesInPDBFile( PDBFileLocation ) ;
		ErrorCheck(ProteinStructure.NumberOfResidues, "when reading the PDB-file in CheckNumberOfResiduesInPDBFile()", ResultsDirectory, WriteLog, logfile) ;

		if ( abs(WriteLog) > 0 )
		{
			fprintf( logfile, "\tRead ProteinStructure.NumberOfAtoms via CheckNumberOfAtomsInPDBFile() from %s\n", PDBFileLocation) ;
			fprintf( logfile, "\n") ;

			fflush( logfile ) ;
		}
		ProteinStructure.NumberOfAtoms = CheckNumberOfAtomsInPDBFile( PDBFileLocation ) ;
		ErrorCheck(ProteinStructure.NumberOfAtoms, "when reading the PDB-file in CheckNumberOfAtomsInPDBFile()", ResultsDirectory, WriteLog, logfile) ;

		if ( abs(WriteLog) > 0 )
		{
			fprintf( logfile, "\tFound %d atoms distributed amongst %d residues\n", ProteinStructure.NumberOfAtoms, ProteinStructure.NumberOfResidues) ;
			fprintf( logfile, "\n") ;
			fprintf( logfile, "\tAllocateProteinStructure()\n") ;
			fprintf( logfile, "\n") ;

			fflush( logfile ) ;
		}

		AllocateProteinStructure( &ProteinStructure, ProteinStructure.NumberOfResidues, ProteinStructure.NumberOfAtoms) ;
	}

	// fill with ModificationName
	strcpy( ProteinStructure.ModificationName, ModificationName) ;

	// import PDB structure
	if (strcmp(PDBFileLocation, "N/A") != 0)
	{
		if ( abs(WriteLog) > 0 )
		{
			fprintf( logfile, "\tImport ProteinStructure.Residues via ImportResiduesFromPDBFile() from %s\n", PDBFileLocation) ;
			fprintf( logfile, "\n") ;

			fflush( logfile ) ;
		}
		ImportResiduesFromPDBFile( PDBFileLocation, &ProteinStructure, ResultsDirectory, WriteLog, logfile) ;

		// really necessary to import Atoms, too, since it has been done in residues ???
		if ( abs(WriteLog) > 0 )
		{
			fprintf( logfile, "\tImportAtomsFromPDBFile()\n") ;
			fprintf( logfile, "\n\n") ;

			fflush( logfile ) ;
		}
		ImportAtomsFromPDBFile( PDBFileLocation, &ProteinStructure, WriteLog, logfile) ;

		// write ProteinStruct to logfile
		if ( abs(WriteLog) > 0 )
		{
			fprintf( logfile, "\tProteinStructure\n") ;
			fprintf( logfile, "\t\t.PDBFileLocation           = %s\n", ProteinStructure.PDBFileLocation) ;
			fprintf( logfile, "\t\t.NumberOfResidues          = %d\n", ProteinStructure.NumberOfResidues) ;
			fprintf( logfile, "\t\t.NumberOfAtoms             = %d\n", ProteinStructure.NumberOfAtoms) ;
			fprintf( logfile, "\t\t.NumberOfHAtoms            = %d\n", ProteinStructure.NumberOfHAtoms) ;
			fprintf( logfile, "\t\t.NumberOfCAtoms            = %d\n", ProteinStructure.NumberOfCAtoms) ;
			fprintf( logfile, "\t\t.NumberOfNAtoms            = %d\n", ProteinStructure.NumberOfNAtoms) ;
			fprintf( logfile, "\t\t.NumberOfOAtoms            = %d\n", ProteinStructure.NumberOfOAtoms) ;
			fprintf( logfile, "\t\t.NumberOfPAtoms            = %d\n", ProteinStructure.NumberOfPAtoms) ;
			fprintf( logfile, "\t\t.NumberOfSAtoms            = %d\n", ProteinStructure.NumberOfSAtoms) ;
			// fprintf( logfile, "\t\t.NumberOfIAtoms         = %d\n", ProteinStructure.NumberOfIAtoms) ;
			fprintf( logfile, "\t\t.NumberOfZNAtoms           = %d\n", ProteinStructure.NumberOfZNAtoms) ;
			fprintf( logfile, "\t\t.NumberOfCLAtoms           = %d\n", ProteinStructure.NumberOfCLAtoms) ;
			fprintf( logfile, "\t\t.NumberOfNAAtoms           = %d\n", ProteinStructure.NumberOfNAAtoms) ;
			fprintf( logfile, "\t\t.NumberOfCAAtoms           = %d\n", ProteinStructure.NumberOfCAAtoms) ;
			fprintf( logfile, "\t\t.NumberOfFEAtoms           = %d\n", ProteinStructure.NumberOfFEAtoms) ;
			fprintf( logfile, "\t\t.NumberOfQAtoms            = %d\n", ProteinStructure.NumberOfQAtoms) ;
			fprintf( logfile, "\t\t.ModificationName          = %s\n", ProteinStructure.ModificationName) ;
			fprintf( logfile, "\t\t.NumberOfModificationAtoms = %d\n", ProteinStructure.NumberOfModificationAtoms) ;
			fprintf( logfile, "\t\t.Weight                    = %lf\n", ProteinStructure.Weight) ;
			fprintf( logfile, "\t\t.Residues[i]\n") ;
			fprintf( logfile, "\t\t\ti\n") ;
			fprintf( logfile, "\t\t\t.Name\n") ;
			fprintf( logfile, "\t\t\t.AtomName\n") ;
			fprintf( logfile, "\t\t\t.ResidueID\n") ;
			fprintf( logfile, "\t\t\t.Volume\n") ;
			fprintf( logfile, "\t\t\t.Weight\n") ;
			fprintf( logfile, "\t\t\t.XRayScatteringLength\n") ;
			fprintf( logfile, "\t\t\t.NeutronScatteringLength\n") ;
			fprintf( logfile, "\t\t\t.xVolume\n") ;
			fprintf( logfile, "\t\t\t.yVolume\n") ;
			fprintf( logfile, "\t\t\t.zVolume\n") ;
			fprintf( logfile, "\t\t\t.xXRayScattering\n") ;
			fprintf( logfile, "\t\t\t.yXRayScattering\n") ;
			fprintf( logfile, "\t\t\t.zXRayScattering\n") ;
			fprintf( logfile, "\t\t\t.xNeutronScattering\n") ;
			fprintf( logfile, "\t\t\t.yNeutronScattering\n") ;
			fprintf( logfile, "\t\t\t.zNeutronScattering\n") ;
			fprintf( logfile, "\n") ;

			LOWER_MAX_INDEX = 0 ; 
			UPPER_MIN_INDEX = 0 ;

			if ( abs(WriteLog) < 2 ) 
			{
				LOWER_MAX_INDEX = min( 7,  ProteinStructure.NumberOfResidues) ; 
				UPPER_MIN_INDEX = max( ProteinStructure.NumberOfResidues - LOWER_MAX_INDEX, LOWER_MAX_INDEX) ;
			}

			for ( int i = 0; i < max( 0, LOWER_MAX_INDEX); ++i)
			{
				fprintf( logfile, "\t\t\t%5d %s %s %d %-8.3f %-8.3f %-.4G %-.4G %-8.3f %-8.3f %-8.3f %-8.3f %-8.3f %-8.3f %-8.3f %-8.3f %-8.3f\n", i, ProteinStructure.Residues[i].Name, ProteinStructure.Residues[i].AtomName, ProteinStructure.Residues[i].ResidueID, ProteinStructure.Residues[i].Volume, ProteinStructure.Residues[i].Weight, ProteinStructure.Residues[i].XRayScatteringLength, ProteinStructure.Residues[i].NeutronScatteringLength, ProteinStructure.Residues[i].xVolume, ProteinStructure.Residues[i].yVolume, ProteinStructure.Residues[i].zVolume, ProteinStructure.Residues[i].xXRayScattering, ProteinStructure.Residues[i].yXRayScattering, ProteinStructure.Residues[i].zXRayScattering, ProteinStructure.Residues[i].xNeutronScattering, ProteinStructure.Residues[i].yNeutronScattering, ProteinStructure.Residues[i].zNeutronScattering) ;
			}
			if ( abs(WriteLog) < 2 ) { fprintf( logfile, "\n\t\t\t... (for full output use -l=2 or -l=-2 option) ...\n\n") ; }
			for ( int i = max( 0, UPPER_MIN_INDEX); i < ProteinStructure.NumberOfResidues; ++i)
			{
				fprintf( logfile, "\t\t\t%5d %s %s %d %-8.3f %-8.3f %-.4G %-.4G %-8.3f %-8.3f %-8.3f %-8.3f %-8.3f %-8.3f %-8.3f %-8.3f %-8.3f\n", i, ProteinStructure.Residues[i].Name, ProteinStructure.Residues[i].AtomName, ProteinStructure.Residues[i].ResidueID, ProteinStructure.Residues[i].Volume, ProteinStructure.Residues[i].Weight, ProteinStructure.Residues[i].XRayScatteringLength, ProteinStructure.Residues[i].NeutronScatteringLength, ProteinStructure.Residues[i].xVolume, ProteinStructure.Residues[i].yVolume, ProteinStructure.Residues[i].zVolume, ProteinStructure.Residues[i].xXRayScattering, ProteinStructure.Residues[i].yXRayScattering, ProteinStructure.Residues[i].zXRayScattering, ProteinStructure.Residues[i].xNeutronScattering, ProteinStructure.Residues[i].yNeutronScattering, ProteinStructure.Residues[i].zNeutronScattering) ;
			}
			fprintf( logfile, "\n\n") ;

			fflush( logfile ) ;
		}
	}





	/// Initialize the user-defined structure
	if ( abs(WriteLog) > 0 )
	{
		ClearScreen( logfile ) ;
		fprintf( logfile, "\n") ;
		fprintf( logfile, "Initialize the user-defined structure:\n") ;
		fprintf( logfile, "\n") ;
		fprintf( logfile, "\tInitializeUserDefinedStructure()\n") ;
		fprintf( logfile, "\n\n") ;

		fflush( logfile ) ;
	}
	InitializeUserDefinedStructure( &UserDefinedStructure ) ;



	/// do a test-computation Constraints for all datasets and write them to logfile
	// not essential to run WIF
	// it is crucial that ComputeConstraints() function does only assign values to Data[i].Constraints, but does not change all other input, e.g. Parameters
	if ( abs(WriteLog) > 0 )
	{
		double* DummyParameters ;
		Initialize1DArray(&DummyParameters, NumberOfParameters) ;
		for ( int i = 0; i < NumberOfParameters; ++i) { DummyParameters[i] = Parameters[i].Value ; }


		ClearScreen( logfile ) ;
		fprintf( logfile, "\n") ;
		fprintf( logfile, "Test computation of constraints for all datasets:\n") ;
		fprintf( logfile, "\n") ;


		for ( int i = 0; i < NumberOfSpectra; ++i)
		{
			fprintf( logfile, "\tComputeConstraints() for Data[%d].Constraints\n", i) ;

			ComputeConstraints( DummyParameters, VolumesOfMolecules, Data[i].ScatteringLengths, Data[i].Contrast, Data[i].Concentration, Data[i].Constraints, ProteinStructure, &UserDefinedStructure) ;
		}
		fprintf( logfile, "\n") ;

		fprintf( logfile, "\tTable\n") ;
		fprintf( logfile, "\t\t.j\n") ;
		for ( int i = 0; i < NumberOfSpectra; ++i)
		{
			fprintf( logfile, "\t\t.Data[%d].Constraints[j]\n", i) ;
		}
		fprintf( logfile, "\n") ;
		for ( int j = 0; j < MaxNumberOfConstraints; ++j)
		{
			fprintf( logfile, "\t\t%2d", j) ;
			for ( int i = 0; i < NumberOfSpectra; ++i)
			{
				fprintf( logfile, " %-.4G", Data[i].Constraints[j]) ;
			}
			fprintf( logfile, "\n") ;
		}
		fprintf( logfile, "\n\n") ;

		fflush( logfile ) ;

		free(DummyParameters) ;
	}


	/// Run fitting routine
	if ( abs(WriteLog) > 0 )
	{
		ClearScreen( logfile ) ;
		fprintf( logfile, "\n") ;
		fprintf( logfile, "Run fitting routine ") ;
	}

	// currently (full) logfile support only for ComputeValue and LevenbergMarquardt routines
	switch ( ChooseFittingRoutine )
	{
		case 0:
		if ( abs(WriteLog) > 0 ) { fprintf( logfile, "no. %d ComputeModel():\n\n", ChooseFittingRoutine) ; fflush( logfile ) ; }
		FittingRoutineError = ComputeModel(Data, NumberOfSpectra, Parameters, NumberOfParameters, &ChiSquare, NumberOfSmearingFolds, VolumesOfMolecules, ProteinStructure, &UserDefinedStructure, TotalNumberOfDatapoints, NumberOfFreeParameters, WriteLog, logfile) ;
		break ;

		case 1:
		// FittingRoutineError = LevenbergMarquardt(Data, NumberOfSpectra, Parameters, NumberOfParameters, FittingRoutineArgument1, &ChiSquare, NumberOfSmearingFolds, VolumesOfMolecules, false, PrintCovarianceMatrix, ProteinStructure, &UserDefinedStructure, DeltaForDifferentiations, NumberOfSampleInformations, TotalNumberOfDatapoints, NumberOfFreeParameters, HighestNumberOfDatapoints) ;
		if ( abs(WriteLog) > 0 ) { fprintf( logfile, "no. %d LevenbergMarquardt():\n\n", ChooseFittingRoutine) ; fflush( logfile ) ; }
		FittingRoutineError = LevenbergMarquardt(Data, NumberOfSpectra, TotalNumberOfDatapoints, HighestNumberOfDatapoints, Parameters, NumberOfParameters, NumberOfFreeParameters, VolumesOfMolecules, NumberOfSampleInformations, ProteinStructure, &UserDefinedStructure, FittingRoutineArgument1, NumberOfSmearingFolds, false, DeltaForDifferentiations, &ChiSquare, PrintCovarianceMatrix, ResultsDirectory, WriteLog, logfile) ;
		break ;

		case 2:
//		FittingRoutineError = GridsearchLM(Data, NumberOfSpectra, Parameters, NumberOfParameters, FittingRoutineArgument1, &ChiSquare, NumberOfSmearingFolds, VolumesOfMolecules, FittingRoutineArgument2, ProteinStructure, &UserDefinedStructure, DeltaForDifferentiations, true, NumberOfSampleInformations, TotalNumberOfDatapoints, NumberOfFreeParameters, HighestNumberOfDatapoints) ;
		if ( abs(WriteLog) > 0 ) { fprintf( logfile, "no. %d GridsearchLM():\n\n", ChooseFittingRoutine) ; fflush( logfile ) ; }
		FittingRoutineError = GridsearchLM(Data, NumberOfSpectra, TotalNumberOfDatapoints, HighestNumberOfDatapoints, Parameters, NumberOfParameters, NumberOfFreeParameters, VolumesOfMolecules, NumberOfSampleInformations, ProteinStructure, &UserDefinedStructure, FittingRoutineArgument1, FittingRoutineArgument2, NumberOfSmearingFolds, DeltaForDifferentiations, &ChiSquare, true, ResultsDirectory, WriteLog, logfile) ;
		break ;

		case 3:
		if ( abs(WriteLog) > 0 ) { fprintf( logfile, "no. %d BFGS():\n\n", ChooseFittingRoutine) ; fflush( logfile ) ; }
		FittingRoutineError = BFGS(Data, NumberOfSpectra, Parameters, NumberOfParameters, FittingRoutineArgument1, &ChiSquare, NumberOfSmearingFolds, VolumesOfMolecules, ProteinStructure, &UserDefinedStructure, DeltaForDifferentiations, false, TotalNumberOfDatapoints, NumberOfFreeParameters) ;
		break ;

		case 4:
		if ( abs(WriteLog) > 0 ) { fprintf( logfile, "no. %d GridsearchBFGS():\n\n", ChooseFittingRoutine) ; fflush( logfile ) ; }
		FittingRoutineError = GridsearchBFGS(Data, NumberOfSpectra, Parameters, NumberOfParameters, FittingRoutineArgument1, &ChiSquare, NumberOfSmearingFolds, VolumesOfMolecules, FittingRoutineArgument2, ProteinStructure, &UserDefinedStructure, DeltaForDifferentiations, true, TotalNumberOfDatapoints, NumberOfFreeParameters) ;
		break ;

		case 5:
		if ( abs(WriteLog) > 0 ) { fprintf( logfile, "no. %d Swarm():\n\n", ChooseFittingRoutine) ; fflush( logfile ) ; }
		FittingRoutineError = Swarm(Data, NumberOfSpectra, Parameters, NumberOfParameters, &ChiSquare, NumberOfSmearingFolds, VolumesOfMolecules, ProteinStructure, &UserDefinedStructure, FittingRoutineArgument3, FittingRoutineArgument2, FittingRoutineArgument1, false, HighestNumberOfDatapoints, NumberOfSampleInformations, TotalNumberOfDatapoints, NumberOfFreeParameters) ;
		break ;

		case 6:
		if ( abs(WriteLog) > 0 ) { fprintf( logfile, "no. %d Genetic():\n\n", ChooseFittingRoutine) ; fflush( logfile ) ; }
		FittingRoutineError = Genetic(Data, NumberOfSpectra, Parameters, NumberOfParameters, &ChiSquare, NumberOfSmearingFolds, VolumesOfMolecules, ProteinStructure,
		&UserDefinedStructure, FittingRoutineArgument3, FittingRoutineArgument2, HighestNumberOfDatapoints, NumberOfSampleInformations, TotalNumberOfDatapoints, NumberOfFreeParameters) ;
		break ;

		case 7:
//		FittingRoutineError = ProfileLikelihood(Data, NumberOfSpectra, Parameters, NumberOfParameters, FittingRoutineArgument1, &ChiSquare, NumberOfSmearingFolds, VolumesOfMolecules, ProteinStructure, &UserDefinedStructure, DeltaForDifferentiations, ChiSquareFractile, FittingRoutineArgument3, CardFileLocation, NumberOfSampleInformations, TotalNumberOfDatapoints, NumberOfFreeParameters, HighestNumberOfDatapoints) ;
		if ( abs(WriteLog) > 0 ) { fprintf( logfile, "no. %d ProfileLikelihood():\n\n", ChooseFittingRoutine) ; fflush( logfile ) ; }
		FittingRoutineError = ProfileLikelihood(CardFileLocation, Data, NumberOfSpectra, TotalNumberOfDatapoints, HighestNumberOfDatapoints, Parameters, NumberOfParameters, NumberOfFreeParameters, VolumesOfMolecules, NumberOfSampleInformations, ProteinStructure, &UserDefinedStructure, FittingRoutineArgument1, NumberOfSmearingFolds, DeltaForDifferentiations, &ChiSquare, ChiSquareFractile, FittingRoutineArgument3, ResultsDirectory, WriteLog, logfile) ;
		break;

		case 8:
//		FittingRoutineError = ProfileLikelihoodSingleParameter(Data, NumberOfSpectra, Parameters, NumberOfParameters, FittingRoutineArgument1, &ChiSquare, NumberOfSmearingFolds, VolumesOfMolecules, ProteinStructure, &UserDefinedStructure, DeltaForDifferentiations, ChiSquareFractile, FittingRoutineArgument3, CardFileLocation, FittingRoutineArgument2, NumberOfSampleInformations, TotalNumberOfDatapoints, NumberOfFreeParameters, HighestNumberOfDatapoints) ;
		if ( abs(WriteLog) > 0 ) { fprintf( logfile, "no. %d ProfileLikelihoodSingleParameter():\n\n", ChooseFittingRoutine) ; fflush( logfile ) ; }
		FittingRoutineError = ProfileLikelihoodSingleParameter(CardFileLocation, Data, NumberOfSpectra, TotalNumberOfDatapoints, HighestNumberOfDatapoints, Parameters, NumberOfParameters, NumberOfFreeParameters, VolumesOfMolecules, NumberOfSampleInformations, ProteinStructure, &UserDefinedStructure, FittingRoutineArgument1, NumberOfSmearingFolds, DeltaForDifferentiations, &ChiSquare, ChiSquareFractile, FittingRoutineArgument3, FittingRoutineArgument2, ResultsDirectory, WriteLog, logfile) ;
		break ;
	}

	ErrorCheck( FittingRoutineError, "when running the selected fitting routine", ResultsDirectory, WriteLog, logfile) ;


	// Conclusion
	if ( abs(WriteLog) > 0 )
	{
		fprintf( logfile, "\tFit complete !\n") ;
		fprintf( logfile, "\n\n") ;

		fflush( logfile ) ;
	}



	/// Output data and parameters
	if ( abs(WriteLog) > 0 )
	{
		ClearScreen( logfile ) ;
		fprintf( logfile, "\n") ;
		fprintf( logfile, "Output data and parameters:\n") ;
		fprintf( logfile, "\n") ;
		fprintf( logfile, "\tSummarize fit results and parameters in %sResults.wif via OutputData()\n", ResultsDirectory) ;
		fprintf( logfile, "\n") ;

		fflush( logfile ) ;
	}
	// OutputData is part of the Model !
	OutputData( ChiSquare, QMin, QMax, Parameters, NumberOfParameters, Data, NumberOfSpectra, CardFileLocation, ProteinStructure, UserDefinedStructure, SamplesFileLocation, ResultsDirectory) ;

	if ( abs(WriteLog) > 0 )
	{
		fprintf( logfile, "\tWrite experimental and fitted spectra via OutputSpectra()\n") ;
		fprintf( logfile, "\n") ;

		fflush( logfile ) ;
	}
	ErrorCheck( OutputSpectra(Data, NumberOfSpectra, ResultsDirectory, WriteLog, logfile), "when writing the spectras in OutputSpectra()", ResultsDirectory, WriteLog, logfile) ;

	if ( abs(WriteLog) > 0 )
	{
		fprintf( logfile, "\tWrite fitted parameters in parameter file and to logfile via OutputParameters()\n") ;
		fprintf( logfile, "\n") ;

		fflush( logfile ) ;
	}
	ErrorCheck( OutputParameters(Parameters, NumberOfParameters, ChooseFittingRoutine, ChiSquare, ResultsDirectory, WriteLog, logfile), "when writing the fitted parameters in OutputParameters()", ResultsDirectory, WriteLog, logfile) ;

	if ( abs(WriteLog) > 0 ) { fprintf( logfile, "\n\n") ; }



	/// Free variables
	if ( abs(WriteLog) > 0 )
	{
		ClearScreen( logfile ) ;
		fprintf( logfile, "\n") ;
		fprintf( logfile, "Free variables\n") ;
		fprintf( logfile, "\n\n") ;

		fflush( logfile ) ;
	}

	// Free dynamically allocated arrays in Data struct
	for ( int i = 0; i < NumberOfSpectra; ++i)
	{
		free(Data[i].QValues) ;
		free(Data[i].IValues) ;
		free(Data[i].FitValues) ;
		free(Data[i].SigmaValues) ;
		free(Data[i].SigmaQValues) ;
		free(Data[i].Constraints) ;
		free(Data[i].ScatteringLengths) ;

		for ( int j = 0; j < HighestNumberOfDatapoints; ++j) { free(Data[i].ResolutionWeights[j]) ; }

		free(Data[i].ResolutionWeights) ;
	}
	free(Data) ;

	// Free other dynamically allocated arrays
	free(VolumesOfMolecules) ;
	free(Parameters) ;
	FreeUserDefined( &UserDefinedStructure ) ;

	if ( strcmp(PDBFileLocation, "N/A") != 0 )
	{
		free(ProteinStructure.Residues) ;
		free(ProteinStructure.Atoms) ;
	}



	/// Return the chisquare if the algorithm executes correctly
	sprintf( Message, "%g", ChiSquare) ;
	if ( !CMD ) { ReturnMessage( Message, ResultsDirectory) ; }

	if ( abs(WriteLog) > 0 ) { fclose(logfile) ; }

	return 0;
}
