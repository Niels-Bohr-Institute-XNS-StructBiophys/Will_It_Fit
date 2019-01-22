/* void OutputData(double ChiSquare, double QMin, double QMax, struct Parameter * Parameters, int NumberOfParameters,
                struct Dataset * Data, int NumberOfSpectra, char cardfilename[128], struct Protein ProteinStructure,
                struct UserDefined UserDefinedStructure, char SampleFilename[256], char *ResultsDirectory) */
/* void OutputData(double ChiSquare, double QMin, double QMax, struct Parameter * Parameters, int NumberOfParameters,
                struct Dataset * Data, int NumberOfSpectra, char cardfilename[256], struct Protein ProteinStructure,
                struct UserDefined UserDefinedStructure, char SampleFilename[256], char *ResultsDirectory)*/

// changed Nicholas version such that char arrays are not static in size

// basically only Parameter, Data, ProteinStructure, UserDefinedStructure are provided 
void OutputData(double ChiSquare, double QMin, double QMax, struct Parameter * Parameters, int NumberOfParameters, struct Dataset * Data, int NumberOfSpectra, char *cardfilename, struct Protein ProteinStructure, struct UserDefined UserDefinedStructure, char *SampleFilename, char *ResultsDirectory)
{
	/// Declarations

	// Variables describing the output file
	FILE *fp ;
	//const char *filename ;
	char Filename[256] ;
	sprintf(Filename, "%sResults.wif", ResultsDirectory) ;

	/// I/O
	// Create the output file
	fp = fopen(Filename, "w+") ;

	// Print filenames to file
	fprintf( fp, ".card-file:\n") ;
	fprintf( fp, "\t%s \n", cardfilename) ;

	fprintf( fp, "\n") ;
	fprintf( fp, "Datafiles:\n") ;

	for( int i = 0; i < NumberOfSpectra; ++i) { fprintf(fp, "\t%s \n", Data[i].Filename) ; }
	fprintf( fp, "\n") ;

	// Print location of sample file
	fprintf( fp, "Sample info-file:\n") ;
	fprintf( fp, "\t%s \n", SampleFilename) ;
	fprintf( fp, "\n") ;

	// Print fit quality to file
	fprintf( fp, "Final chisq = %g \n", ChiSquare) ;
	fprintf( fp, "\n") ;

	// Print range of q to file
	fprintf( fp, "Lower limit on q = %g \n", QMin) ;
	fprintf( fp, "Upper limit on q = %g \n", QMax) ;
	fprintf( fp, "\n") ;

	// Print parameters and properties of parameters to file
	fprintf( fp, "Parameters:\n\n") ;
	fprintf( fp, "\t     Value             Error             Name\n") ;
	for( int i = 0; i < NumberOfParameters; ++i)
	{
		if ( Parameters[i].iParameter == true ) { fprintf( fp, "\t%2d   %-15g   %-15g   %s\n", i, Parameters[i].Value, Parameters[i].Error, Parameters[i].Name) ; }
		else { fprintf( fp, "\t%2d   %-15g   Fixed             %s\n", i, Parameters[i].Value, Parameters[i].Name) ; }
	}
	fprintf( fp, "\n\n") ;





	// Print info on protein file

	// scaling with HYDR (for WAT) and GLYCV (for MOD) must applied here, since the displaced volumes in the ProteinStructure are not changed (in Model() each time a copy CurrentAtom is created for the scaled version of the current atom)

	fprintf( fp, "Atom and residue counts:\n") ;
	fprintf( fp, "\n") ;
	fprintf( fp, "\tNumber of Atoms                    : %d\n", ProteinStructure.NumberOfAtoms) ;
	fprintf( fp, "\tNumber of Atoms (Modification %s) : %d\n", ProteinStructure.ModificationName, ProteinStructure.NumberOfModificationAtoms) ;
	fprintf( fp, "\tNumber of Residues                 : %d\n", ProteinStructure.NumberOfResidues) ;
	fprintf( fp, "\n") ;
	fprintf( fp, "\tNumber of H  : %d\n", ProteinStructure.NumberOfHAtoms) ;
	fprintf( fp, "\tNumber of D  : %d\n", ProteinStructure.NumberOfDAtoms) ;
	fprintf( fp, "\tNumber of C  : %d\n", ProteinStructure.NumberOfCAtoms) ;
	fprintf( fp, "\tNumber of N  : %d\n", ProteinStructure.NumberOfNAtoms) ;
	fprintf( fp, "\tNumber of O  : %d\n", ProteinStructure.NumberOfOAtoms) ;
	fprintf( fp, "\tNumber of P  : %d\n", ProteinStructure.NumberOfPAtoms) ;
	fprintf( fp, "\tNumber of S  : %d\n", ProteinStructure.NumberOfSAtoms) ;
	// fprintf( fp, "Number of I: %d\n", ProteinStructure.NumberOfIAtoms) ;
	fprintf( fp, "\n") ;
	fprintf( fp, "\tNumber of Zn : %d\n", ProteinStructure.NumberOfZNAtoms) ;
	fprintf( fp, "\tNumber of Cl : %d\n", ProteinStructure.NumberOfCLAtoms) ;
	fprintf( fp, "\tNumber of Na : %d\n", ProteinStructure.NumberOfNAAtoms) ;
	fprintf( fp, "\tNumber of Ca : %d\n", ProteinStructure.NumberOfCAAtoms) ;
	fprintf( fp, "\tNumber of Fe : %d\n", ProteinStructure.NumberOfFEAtoms) ;
	fprintf( fp, "\n") ;
	fprintf( fp, "\tNumber of Q  : %d\n", ProteinStructure.NumberOfQAtoms) ;
	fprintf( fp, "\n\n") ;



	// Calculate molecular weigths, displaced volumes, scattering length densities
	// these quantities do not depend on the spectra, i.e. the contrasts
	double Weight = 0.0 ;
	double WeightHOH = 0.0 ;
	double WeightWAT = 0.0 ;
	double WeightWater = 0.0 ;
	double WeightModification = 0.0 ;
	double WeightWithoutWater = 0.0 ;
	double WeightWithoutWaterAndModification = 0.0 ;

	double Volume = 0.0 ;
	double VolumeHOH = 0.0 ;
	double VolumeWAT = 0.0 ;
	double VolumeWater = 0.0 ;
	double VolumeModification = 0.0 ;
	double VolumeWithoutWater = 0.0 ;
	double VolumeWithoutWaterAndModification = 0.0 ;

	for ( int i = 0; i < ProteinStructure.NumberOfAtoms; ++i)
	{
		Volume += ProteinStructure.Atoms[i].Volume ;
		Weight += ProteinStructure.Atoms[i].Weight ;

		if( strcmp( ProteinStructure.Atoms[i].ResidueName, "WAT") == 0 )
		{
			VolumeWAT += ProteinStructure.Atoms[i].Volume * ( 1./(1.+Parameters[HYDR].Value) ) ;
			WeightWAT += ProteinStructure.Atoms[i].Weight ;
		}
		if( strcmp( ProteinStructure.Atoms[i].ResidueName, "HOH") == 0 )
		{
			VolumeHOH += ProteinStructure.Atoms[i].Volume * ( 1./(1.+Parameters[HYDR].Value) ) ;
			WeightHOH += ProteinStructure.Atoms[i].Weight ;
		}
	 	if( strcmp(ProteinStructure.Atoms[i].ResidueName, ProteinStructure.ModificationName) == 0 )
		{
			VolumeModification += ProteinStructure.Atoms[i].Volume * ( 1./(1.+Parameters[GLYCV].Value) ) ;
			WeightModification += ProteinStructure.Atoms[i].Weight ;
		}
	}

	WeightWater = WeightWAT + WeightHOH ;

	WeightWithoutWater = Weight - WeightWater ;
	WeightWithoutWaterAndModification = Weight - WeightWater - WeightModification ;


	VolumeWater = VolumeWAT + VolumeHOH ;

	VolumeWithoutWater = Volume - VolumeWater ;
	VolumeWithoutWaterAndModification = Volume - VolumeWater - VolumeModification ;

	fprintf(fp, "Molecular weights [u] and displaced volumes [A^3]:\n") ;
	fprintf( fp, "\n") ;
	fprintf(fp, "\tMolecular weight                             : %lf\n", Weight) ;
	fprintf(fp, "\tMolecular weight (only WAT)                  : %lf\n", WeightWAT) ;
	fprintf(fp, "\tMolecular weight (only HOH)                  : %lf\n", WeightHOH) ;
	fprintf(fp, "\tMolecular weight (only WAT + HOH)            : %lf\n", WeightWater) ;
	fprintf(fp, "\tMolecular weight (only %s)                  : %lf\n", ProteinStructure.ModificationName, WeightModification) ;
	fprintf(fp, "\tMolecular weight (without WAT + HOH)         : %lf\n", WeightWithoutWater) ;
	fprintf(fp, "\tMolecular weight (without WAT + HOH and %s) : %lf\n", ProteinStructure.ModificationName, WeightWithoutWaterAndModification) ;
	fprintf( fp, "\n") ;
	fprintf(fp, "\tVolume                                       : %lf\n", Volume) ;
	fprintf(fp, "\tVolume (only WAT)                            : %lf\n", VolumeWAT) ;
	fprintf(fp, "\tVolume (only HOH)                            : %lf\n", VolumeHOH) ;
	fprintf(fp, "\tVolume (only WAT + HOH)                      : %lf\n", VolumeWater) ;
	fprintf(fp, "\tVolume (only %s)                            : %lf\n", ProteinStructure.ModificationName, VolumeModification) ;
	fprintf(fp, "\tVolume (without WAT + HOH)                   : %lf\n", VolumeWithoutWater) ;
	fprintf(fp, "\tVolume (without WAT + HOH and %s)           : %lf\n", ProteinStructure.ModificationName, VolumeWithoutWaterAndModification) ;
	fprintf( fp, "\n") ;
	fprintf(fp, "\tSpecific volume of %s (cm^3/g)              : ", ProteinStructure.ModificationName) ;
	if ( WeightModification != 0.0 ) { fprintf(fp, "%lf\n", VolumeModification * 1.0e-24 * 6.022e23 / WeightModification) ; } else { fprintf(fp, "0.0\n") ; }

	fprintf( fp, "\n\n") ;


	// in the following put all computations involving scattering length (densities) since they depend on the spectras contrasts
	// use quantities such as molecular volumes from above
	double Rg2, dsld, sld ;
	double Rg2HOH, dsldHOH, sldHOH ;
	double Rg2WAT, dsldWAT, sldWAT ;
	double Rg2Modification, dsldModification, sldModification ;

	double r0 = 2.818e-13 ; // classical electron radius [cm]

	double Rg2Water, sldWater ;
	double Rg2WithoutWater, sldWithoutWater ;
	double Rg2WithoutWaterAndModification, sldWithoutWaterAndModification ;

	double CoM[3], CoMHOH[3], CoMWAT[3], CoMModification[3], CoMWithoutWater[3], CoMWithoutWaterAndModification[3] ;

	fprintf( fp, "Radii of gyration (Rg [A]), centers of mass (CoM [A]) and scattering length densities (SLD [1e-6 A^-2]), for X-rays also electron densities (ED [e-/nm^-3]) are given in brackets:\n\n") ;
	for ( int i = 0 ; i < NumberOfSpectra ; ++i )
	{
		fprintf( fp, "\tSpectra no. %d:\n\n", i) ;

		// Calculate radii of gyration and CoM (both weighted by excess SLDs) as well as scattering lengths
		// Note that selections in the calculations need to be done in exactly the same way as above 
		// Rg2, CoM, sld and dsld are reset to 0.0 in ComputeRg2AndCoMAndSL at first
		// units:
		// -Data[i].Constraints[SLDWater] are in [cm/Å^3]
		// -ScatteringLengths are in cm
		// -(displaced) atom volumes are in Å^3
		// https://stackoverflow.com/questions/20443699/pass-by-pointer-pass-by-reference-in-c
		ComputeRg2AndCoMAndSL( &Rg2, CoM, &dsld, &sld, ProteinStructure, Parameters, Data[i].Contrast, Data[i].Constraints[SLDWater], "") ;
		ComputeRg2AndCoMAndSL( &Rg2HOH, CoMHOH, &dsldHOH, &sldHOH, ProteinStructure, Parameters, Data[i].Contrast, Data[i].Constraints[SLDWater], "HOH") ;
		ComputeRg2AndCoMAndSL( &Rg2WAT, CoMWAT, &dsldWAT, &sldWAT, ProteinStructure, Parameters, Data[i].Contrast, Data[i].Constraints[SLDWater], "WAT") ;
		// for now Rg2Modification is unused
		ComputeRg2AndCoMAndSL( &Rg2Modification, CoMModification, &dsldModification, &sldModification, ProteinStructure, Parameters, Data[i].Contrast, Data[i].Constraints[SLDWater], ProteinStructure.ModificationName) ;


		// calculate missing Rg and CoM (all in [Å])
		Rg2Water = ( Rg2HOH * dsldHOH + Rg2WAT * dsldWAT ) / ( dsldHOH + dsldWAT ) ;
		Rg2WithoutWater = ( Rg2 * dsld - Rg2HOH * dsldHOH - Rg2WAT * dsldWAT ) / ( dsld - dsldHOH - dsldWAT ) ;
		Rg2WithoutWaterAndModification = ( Rg2 * dsld - Rg2HOH * dsldHOH - Rg2WAT * dsldWAT - Rg2Modification * dsldModification ) / ( dsld - dsldHOH - dsldWAT - dsldModification ) ;

		CoMWithoutWater[0] = ( CoM[0] * dsld - CoMHOH[0] * dsldHOH - CoMWAT[0] * dsldWAT ) / ( dsld - dsldHOH - dsldWAT ) ;
		CoMWithoutWater[1] = ( CoM[1] * dsld - CoMHOH[1] * dsldHOH - CoMWAT[1] * dsldWAT ) / ( dsld - dsldHOH - dsldWAT ) ;
		CoMWithoutWater[2] = ( CoM[2] * dsld - CoMHOH[2] * dsldHOH - CoMWAT[2] * dsldWAT ) / ( dsld - dsldHOH - dsldWAT ) ;

		CoMWithoutWaterAndModification[0] = ( CoM[0] * dsld - CoMHOH[0] * dsldHOH - CoMWAT[0] * dsldWAT - CoMModification[0] * dsldModification ) / ( dsld - dsldHOH - dsldWAT - dsldModification ) ;
		CoMWithoutWaterAndModification[1] = ( CoM[1] * dsld - CoMHOH[1] * dsldHOH - CoMWAT[1] * dsldWAT - CoMModification[1] * dsldModification ) / ( dsld - dsldHOH - dsldWAT - dsldModification ) ;
		CoMWithoutWaterAndModification[2] = ( CoM[2] * dsld - CoMHOH[2] * dsldHOH - CoMWAT[2] * dsldWAT - CoMModification[2] * dsldModification ) / ( dsld - dsldHOH - dsldWAT - dsldModification ) ;


		// calculate SLDs (all in [1e-6 Å^-2]) from scattering lengths [cm=1e8 Å] and volumes [Å^3]
		// first calculate missing scattering lengths
		sldWater = sldWAT + sldHOH ;
 		sldWithoutWater = sld - sldWater ;
		sldWithoutWaterAndModification = sld - sldWater - sldModification ;

		// now scale all scattering lengths accordingly with their volume to obtain the SLDs
		if ( Volume != 0.0 ) { sld *= ( 1.0e14 / Volume ) ; } else { sld = 0.0 ; }
		if ( VolumeHOH != 0.0 ) { sldHOH *= ( 1.0e14 / VolumeHOH ) ; } else { sldHOH = 0.0 ; }
		if ( VolumeWAT != 0.0 ) { sldWAT *= ( 1.0e14 / VolumeWAT ) ; } else { sldWAT = 0.0 ; }
		if ( VolumeWater != 0.0 ) { sldWater *= ( 1.0e14 / VolumeWater ) ; } else { sldWater = 0.0 ; }
		if ( VolumeModification != 0.0 ) { sldModification *= ( 1.0e14 / VolumeModification ) ; } else { sldModification = 0.0 ; }

		if ( VolumeWithoutWater != 0.0 ) { sldWithoutWater *= ( 1.0e14 / VolumeWithoutWater ) ; } else { sldWithoutWater = 0.0 ; }
		if ( VolumeWithoutWaterAndModification != 0.0 ) { sldWithoutWaterAndModification *= ( 1.0e14 / VolumeWithoutWaterAndModification ) ; } else { sldWithoutWaterAndModification = 0.0 ; }


		// write information
		if ( Data[i].Contrast < 0.0 ) { fprintf( fp, "\t\tX-ray                ") ; }
		else { fprintf( fp, "\t\tNeutron (%5.1lf%% D2O) ", Data[i].Contrast) ; }
		fprintf( fp, "Rg                              : %7.2lf\n", sqrt(Rg2)) ;

		if ( Data[i].Contrast < 0.0 ) { fprintf( fp, "\t\tX-ray                ") ; }
		else { fprintf( fp, "\t\tNeutron (%5.1lf%% D2O) ", Data[i].Contrast) ; }
		fprintf( fp, "Rg  (without WAT + HOH)         : %7.2lf\n", sqrt(Rg2WithoutWater)) ;

		if ( Data[i].Contrast < 0.0 ) { fprintf( fp, "\t\tX-ray                ") ; }
		else { fprintf( fp, "\t\tNeutron (%5.1lf%% D2O) ", Data[i].Contrast) ; }
		fprintf( fp, "Rg  (without WAT + HOH and %s) : %7.2lf\n", ProteinStructure.ModificationName, sqrt(Rg2WithoutWaterAndModification)) ;

		if ( Data[i].Contrast < 0.0 ) { fprintf( fp, "\t\tX-ray                ") ; }
		else { fprintf( fp, "\t\tNeutron (%5.1lf%% D2O) ", Data[i].Contrast) ; }
		fprintf( fp, "Rg  (only WAT + HOH)            : %7.2lf\n", sqrt(Rg2Water)) ;

		fprintf( fp, "\n") ;


		if ( Data[i].Contrast < 0.0 ) { fprintf( fp, "\t\tX-ray                ") ; }
		else { fprintf( fp, "\t\tNeutron (%5.1lf%% D2O) ", Data[i].Contrast) ; }
		fprintf( fp, "SLD                             : %7.3lf", sld) ;
		if ( Data[i].Contrast < 0.0 ) { fprintf( fp, " (ED = %7.2lf)\n", (sld/r0*1.0e-11)) ; } else { fprintf( fp, "\n") ; }

		if ( Data[i].Contrast < 0.0 ) { fprintf( fp, "\t\tX-ray                ") ; }
		else { fprintf( fp, "\t\tNeutron (%5.1lf%% D2O) ", Data[i].Contrast) ; }
		fprintf( fp, "SLD (only WAT)                  : %7.3lf", sldWAT) ;
		if ( Data[i].Contrast < 0.0 ) { fprintf( fp, " (ED = %7.2lf)\n", (sldWAT/r0*1.0e-11)) ; } else { fprintf( fp, "\n") ; }

		if ( Data[i].Contrast < 0.0 ) { fprintf( fp, "\t\tX-ray                ") ; }
		else { fprintf( fp, "\t\tNeutron (%5.1lf%% D2O) ", Data[i].Contrast) ; }
		fprintf( fp, "SLD (only HOH)                  : %7.3lf", sldHOH) ;
		if ( Data[i].Contrast < 0.0 ) { fprintf( fp, " (ED = %7.2lf)\n", (sldHOH/r0*1.0e-11)) ; } else { fprintf( fp, "\n") ; }

		if ( Data[i].Contrast < 0.0 ) { fprintf( fp, "\t\tX-ray                ") ; }
		else { fprintf( fp, "\t\tNeutron (%5.1lf%% D2O) ", Data[i].Contrast) ; }
		fprintf( fp, "SLD (sample solvent)            : %7.3lf", 1.0e14*Data[i].Constraints[SLDWater]) ;
		if ( Data[i].Contrast < 0.0 ) { fprintf( fp, " (ED = %7.2lf)\n", (Data[i].Constraints[SLDWater]/r0*1.0e3)) ; } else { fprintf( fp, "\n") ; }

		if ( Data[i].Contrast < 0.0 ) { fprintf( fp, "\t\tX-ray                ") ; }
		else { fprintf( fp, "\t\tNeutron (%5.1lf%% D2O) ", Data[i].Contrast) ; }
		fprintf( fp, "SLD (only %s)                  : %7.3lf", ProteinStructure.ModificationName, sldModification) ;
		if ( Data[i].Contrast < 0.0 ) { fprintf( fp, " (ED = %7.2lf)\n", (sldModification/r0*1.0e-11)) ; } else { fprintf( fp, "\n") ; }

		if ( Data[i].Contrast < 0.0 ) { fprintf( fp, "\t\tX-ray                ") ; }
		else { fprintf( fp, "\t\tNeutron (%5.1lf%% D2O) ", Data[i].Contrast) ; }
		fprintf( fp, "SLD (without WAT + HOH)         : %7.3lf", sldWithoutWater) ;
		if ( Data[i].Contrast < 0.0 ) { fprintf( fp, " (ED = %7.2lf)\n", (sldWithoutWater/r0*1.0e-11)) ; } else { fprintf( fp, "\n") ; }

		if ( Data[i].Contrast < 0.0 ) { fprintf( fp, "\t\tX-ray                ") ; }
		else { fprintf( fp, "\t\tNeutron (%5.1lf%% D2O) ", Data[i].Contrast) ; }
		fprintf( fp, "SLD (without WAT + HOH and %s) : %7.3lf", ProteinStructure.ModificationName, sldWithoutWaterAndModification) ;
		if ( Data[i].Contrast < 0.0 ) { fprintf( fp, " (ED = %7.2lf)\n", (sldWithoutWaterAndModification/r0*1.0e-11)) ; } else { fprintf( fp, "\n") ; }

		fprintf( fp, "\n") ;


		if ( Data[i].Contrast < 0.0 ) { fprintf( fp, "\t\tX-ray                ") ; }
		else { fprintf( fp, "\t\tNeutron (%5.1lf%% D2O) ", Data[i].Contrast) ; }
		fprintf( fp, "CoM                             : (%lf %lf %lf)\n", CoM[0], CoM[1], CoM[2]) ;

		if ( Data[i].Contrast < 0.0 ) { fprintf( fp, "\t\tX-ray                ") ; }
		else { fprintf( fp, "\t\tNeutron (%5.1lf%% D2O) ", Data[i].Contrast) ; }
		fprintf( fp, "CoM (without WAT + HOH)         : (%lf %lf %lf)\n", CoMWithoutWater[0], CoMWithoutWater[1], CoMWithoutWater[2]) ;

		if ( Data[i].Contrast < 0.0 ) { fprintf( fp, "\t\tX-ray                ") ; }
		else { fprintf( fp, "\t\tNeutron (%5.1lf%% D2O) ", Data[i].Contrast) ; }
		fprintf( fp, "CoM (without WAT + HOH and %s) : (%lf %lf %lf)\n", ProteinStructure.ModificationName, CoMWithoutWaterAndModification[0], CoMWithoutWaterAndModification[1], CoMWithoutWaterAndModification[2]) ;

		fprintf( fp, "\n") ;
	}


	// Close file and end program
	fclose(fp) ;
}
