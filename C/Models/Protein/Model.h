double Model(double q, double * Parameters, double * Constraints, double Contrast, struct Protein ProteinStructure, struct UserDefined * UserDefinedStructure)
{
	// Conc and SLDWater indices for CONSTRAINTS are globally defined in Constraints.h
	// GLYCV, HYDR etc indices for PARAMETERS are globally defined in Constraints.h


	// Scaling and background comes via Parameters array
	// Concentration of sample comes via Constraints which previously was assigned with Data[i].Concentration


	// Variables describing the sample
	struct Atom CurrentAtom;
	double Intensity;
	double Scaling;
	double Background;

	// get current sample concentration from Constraints array
	double ConcentrationOfSample          = Constraints[Conc] ; // [1/cm^3]
	// get current sample SLD of solvent from Constraints array
	// e.g. for water
	// X-rays		9.412E-14 [cm/Å^3]
	// neutrons  42% D2O	2.361E-14 [cm/Å^3]
	// neutrons 100% D2O	6.394E-14 [cm/Å^3] = 6.394 [10^-6 Å^-2]
	double ScatteringLengthDensityOfWater = Constraints[SLDWater] ; // [cm/Å^3]

	double complex ** Beta = ComplexArray( NumberOfHarmonics + 1, NumberOfHarmonics + 1) ;

	//Loop over residues
	for ( int j = 0; j < ProteinStructure.NumberOfAtoms; j++)
	{
		// for AA based calculation
		CopyAtom( &CurrentAtom, &ProteinStructure.Atoms[j]) ;

		// treat dummy waters (WAT residues with atoms Q), crystal waters (HOH, atoms H and O) and modifications separately

		// dummy waters (WAT residues, atoms Q)
		if ( strcmp( CurrentAtom.ResidueName, "WAT") == 0 )
		{
			CurrentAtom.Volume = CurrentAtom.Volume * ( 1./(1.+Parameters[HYDR]) ) ;

			if ( Contrast < 0.0 )
			{
				// for X-rays is has been assigned already in ImportPDBFile.h
				// e.g. 4.133 * 10 * 2.818 e-13 = 1.165 e-11 [cm]
			}
			else // ( Contrast >= 0.0 && Contrast <= 100.0)
			{
				// for neutron scattering update now here the scattering length of the waters with the contrast value
				// e.g. Contrast =   42 -> 2.9208e-12 [cm]
				// e.g. Contrast =  100 -> 7.9126e-12 [cm]
				CurrentAtom.NeutronScatteringLength = 4.133 *( 2.0 *((Contrast* 6.671e-13 / 100.) + ((100.-Contrast)* -3.741e-13 / 100.)) + 5.803e-13) ;
			}
		}
		// crystal waters (residues HOH, atoms H and O) consist of and are treated like normal atoms, i.e. 2xH + 1xO atoms
		// scale volume for both atom types, adapt NeutronScatteringLength for hydrogen according to contrast
		if ( strcmp( CurrentAtom.ResidueName, "HOH") == 0 )
		{
			CurrentAtom.Volume = CurrentAtom.Volume * ( 1./(1.+Parameters[HYDR]) ) ;

			if ( Contrast < 0.0 )
			{
				// for X-rays is has been assigned already in ImportPDBFile.h
			}
			else // ( Contrast >= 0.0 && Contrast <= 100.0)
			{
				CurrentAtom.NeutronScatteringLength = ( Contrast * 6.671e-13 + (100.-Contrast) * -3.741e-13 ) / 100. ;
			}
		}
		// modification
		if ( strcmp(CurrentAtom.ResidueName, ProteinStructure.ModificationName) == 0 )
		{
			CurrentAtom.Volume = CurrentAtom.Volume * ( 1./(1.+Parameters[GLYCV]) ) ;
		}

		AddScatteringFromAtom( Beta, q, CurrentAtom, Contrast, ScatteringLengthDensityOfWater, 1.0) ;
		//AddScatteringFromSolvent(BetaSolvent, q, CurrentAtom, Contrast, ScatteringLengthDensityOfWater, 1.0); // Parameters[PROTSCALE]);
	}

	// Calculate intensity
	Intensity = 0.0;
	for (int l = 0; l < NumberOfHarmonics + 1; l++)
	{
		for (int m = 0; m <= l; m++)
		{
			Intensity += ((m > 0) + 1) * pow(cabs(Beta[l][m]), 2) ;
			//IntensitySolvent += ((m > 0) + 1) * pow(cabs(BetaSolvent[l][m]), 2);
		}
	}

	// Rescale result and return
	// if (Contrast >= 0.0 && Contrast <= 100.0)
	if ( Contrast < 0.0 )
	{
		Scaling    = Parameters[SCALECONC] * Parameters[SCALEX] ;
		Background = Parameters[BACKX] ;
	}
	else
	{
		Scaling    = Parameters[SCALECONC] * (Contrast / 100.0 * Parameters[SCALEN100] + (100.0 - Contrast) / 100.0 * Parameters[SCALEN0]) ;
		Background = Contrast / 100.0 * Parameters[BACKN100] + (100.0 - Contrast) / 100.0 * Parameters[BACKN0] ;
	}

	Intensity = Scaling *  Intensity * ConcentrationOfSample ; //;/pow(2.818E-13,2)

	Intensity = Intensity - Background ;

	// printf( "%lf, %g\n", q, Intensity) ;

	return Intensity ;
}
