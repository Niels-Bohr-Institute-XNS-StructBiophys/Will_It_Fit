int Sign(double x) {
	int Result;

	if (x < 0.0) { Result = -1 ; }
	else if (x > 0.0) { Result = 1 ; }
	else { Result = 0; }

	return Result;
}


double complex PolarComplexNumber(double Radius, double Phi)
{
	if (Phi == 0)
		 return Radius + I * 0.0;
	else
		 return Radius * (cos(Phi) + I * sin(Phi));
}


// computes the Rg^2 and CoM from the given ProteinStructure.Atoms
// if ResidueMatch is not empty, only atoms belonging to a residue whose name matches ResidueMatch will be considered in the sums
void ComputeRg2AndCoMAndSL( double * Rg2, double * CoM, double * dsld, double * sl, struct Protein ProteinStructure, struct Parameter * Parameters, double Contrast, double ScatteringLengthDensityOfSolvent, char const * ResidueMatch)
{
	double x, y, z, r2, V ;
	double ScatteringLength ;

	// reset to 0.0
	*Rg2 = 0.0 ;

	CoM[0] = 0.0 ;
	CoM[1] = 0.0 ;
	CoM[2] = 0.0 ;

	*dsld = 0.0 ;
	*sl = 0.0 ;

	for ( int i = 0; i< ProteinStructure.NumberOfResidues; ++i)
	{
		// skip atoms of not matching residues, select for example only WAT, HOH or modifications
		if ( strcmp( ResidueMatch, "") != 0 )
		{
			if ( strcmp( ProteinStructure.Residues[i].Name, ResidueMatch) != 0 ) { continue ; }
		}

		if ( Contrast < 0.0 )
		{
			x = ProteinStructure.Residues[i].xXRayScattering ;
			y = ProteinStructure.Residues[i].yXRayScattering ;
			z = ProteinStructure.Residues[i].zXRayScattering ;

			ScatteringLength = ProteinStructure.Residues[i].XRayScatteringLength ;
		}
		else
		{
			x = ProteinStructure.Residues[i].xNeutronScattering ;
			y = ProteinStructure.Residues[i].yNeutronScattering ;
			z = ProteinStructure.Residues[i].zNeutronScattering ;

			ScatteringLength = ProteinStructure.Residues[i].NeutronScatteringLength ; 

			// note that crystal waters (residues HOH) consist of and are treated like normal atoms, i.e. O and 2*H atoms
			// so far only dummy waters are adjusted with the right contrast 
			if ( strcmp( ProteinStructure.Residues[i].Name, "WAT") == 0 )
			{
				ScatteringLength = 4.133 *( 2.0 *((Contrast* 6.671e-13 / 100.) + ((100.-Contrast)* -3.741e-13 / 100.)) + 5.803e-13) ;
			}
			// if ( strcmp( ProteinStructure.Residues[i].Name, "HOH") == 0 )
			// {
			// 	ScatteringLength = 2.0 * ((Contrast* 6.671e-13 / 100.) + ((100.-Contrast)* -3.741e-13 / 100.)) + 5.803e-13 ;
			// }
		}

		r2  = x * x + y * y + z * z ;

		V  = ProteinStructure.Residues[i].Volume ;

		// Scale certain volumes 
		// so far only dummy waters are scaled
		if ( strcmp( ProteinStructure.Residues[i].Name, "WAT") == 0 )
		{
			V *= ( 1./(1.+Parameters[HYDR].Value) ) ;
		}
		// if ( strcmp( ProteinStructure.Residues[i].Name, "WAT") == 0 || strcmp( ProteinStructure.Residues[i].Name, "HOH") == 0 )
		// {
		// 	V *= ( 1./(1.+Parameters[HYDR].Value) ) ;
		// }
		if ( strcmp( ProteinStructure.Residues[i].Name, "  X") == 0 )
		{
			V *= ( 1./(1.+Parameters[GLYCV].Value) ) ;
		}

		*Rg2 += ( ScatteringLength / V - ScatteringLengthDensityOfSolvent ) * r2 ;

		CoM[0] += ( ScatteringLength / V - ScatteringLengthDensityOfSolvent ) * x ;
		CoM[1] += ( ScatteringLength / V - ScatteringLengthDensityOfSolvent ) * y ;
		CoM[2] += ( ScatteringLength / V - ScatteringLengthDensityOfSolvent ) * z ;

		*dsld += ( ScatteringLength / V - ScatteringLengthDensityOfSolvent ) ;
		*sl += ScatteringLength ;
	}

	// set those quantities affected / weighted by dsld to 0.0 
	if ( *dsld != 0.0 ) { *Rg2 /= *dsld ; CoM[0] /= *dsld ; CoM[1] /= *dsld ; CoM[2] /= *dsld ; } else { *Rg2 = 0.0 ; CoM[0] = 0.0 ; CoM[1] = 0.0 ; CoM[2] = 0.0 ; }
}
//NEW if ImportPDBFile has been fixed:
/*
void ComputeRg2AndCoM( double Rg2, double * CoM, double dsld, struct Protein ProteinStructure, struct Parameter * Parameters, double Contrast, double ScatteringLengthDensityOfSolvent, char* ResidueMatch)
{
	double x, y, z, r2, V ;
	double ScatteringLength ;

	// reset to 0.0
	*Rg2 = 0.0 ;
	CoM[0] = 0.0 ;
	CoM[1] = 0.0 ;
	CoM[2] = 0.0 ;
	*dsld = 0.0 ;

	for ( int i = 0; i< ProteinStructure.NumberOfAtoms)
	{
		// skip atoms of matching residues such as WAT, HOH or modifications
		if ( strcmp( ResidueMatch, "") != 0 )
		{
			if ( strcmp( ProteinStructure.Atoms[i].ResidueName, ResidueMatch) != 0 ) { continue ; }
		}

		if ( Contrast < 0.0 ) 
		{
			ScatteringLength = ProteinStructure.Atoms[i].XRayScatteringLength ;
		}
		else
		{
			ScatteringLength = ProteinStructure.Atoms[i].NeutronScatteringLength ; 

			// note that crystal waters (residues HOH) consist of and are treated like normal atoms, i.e. O and 2*H atoms
			// so far only dummy waters are adjusted with the right contrast 
			if ( strcmp( ProteinStructure.Atoms[i].ResidueName, "WAT") == 0 )
			{
				ScatteringLength = 4.133 *( 2.0 *((Contrast* 6.671e-13 / 100.) + ((100.-Contrast)* -3.741e-13 / 100.)) + 5.803e-13) ;
			}
			// if ( strcmp( ProteinStructure.Atoms[i].ResidueName, "HOH") == 0 )
			// {
			// 	ScatteringLength = 2.0 * ((Contrast* 6.671e-13 / 100.) + ((100.-Contrast)* -3.741e-13 / 100.)) + 5.803e-13 ;
			// }
		}

		x = ProteinStructure.Atoms[i].x ;
		y = ProteinStructure.Atoms[i].y ;
		z = ProteinStructure.Atoms[i].z ;

		r2  = x * x + y * y + z * z ;

		V  = ProteinStructure.Atoms[i].Volume ;
		// so far only dummy waters are scaled 
		if ( strcmp( ProteinStructure.Atoms[i].ResidueName, "WAT") == 0 )
		{
			V *= ( 1./(1.+Parameters[HYDR].Value) ) ;
		}
		// if ( strcmp( ProteinStructure.Atoms[i].ResidueName, "WAT") == 0 || strcmp( ProteinStructure.Atoms[i].ResidueName, "HOH") == 0 )
		// {
		// 	V *= ( 1./(1.+Parameters[HYDR].Value) ) ;
		// }
		if ( strcmp( ProteinStructure.Atoms[i].ResidueName, "  X") == 0 )
		{
			V *= ( 1./(1.+Parameters[GLYCV].Value) ) ;
		}

		*Rg2 += ( ScatteringLength / V - ScatteringLengthDensityOfSolvent ) * r2 ;
		CoM[0] += ( ScatteringLength / V - ScatteringLengthDensityOfSolvent ) * x ;
		CoM[1]+= ( ScatteringLength / V - ScatteringLengthDensityOfSolvent ) * y ;
		CoM[2] += ( ScatteringLength / V - ScatteringLengthDensityOfSolvent ) * z ;
		*dsld += ( ScatteringLength / V - ScatteringLengthDensityOfSolvent ) ;
	}

	if ( *dsld != 0.0 ) { *Rg2 /= *dsld ; CoM[0] /= *dsld ; CoM[1] /= *dsld ; CoM[2] /= *dsld ; } else { 	*Rg2 = 0.0 ; CoM[0] = 0.0 ; CoM[1] = 0.0 ; CoM[2] = 0.0 ; }
*/



void AddScatteringFromResidue(double complex **Beta, double q, struct Residue CurrentResidue, double Contrast, double ScatteringLengthDensityOfSolvent, double DeltaB) {
	// Søren Kynde, 2012 (rewritten by Martin Cramer Pedersen, 2015)
	// This function adds the scattering of a residue to the scattering amplitude expandended by the spherical harmonic coefficients Beta_lm.
	int l ;
	int m ;

	double x ;
	double y ;
	double z ;

	double Radius ;
	double Theta ;
	double Phi ;

	double ScatteringLengthOfResidue ;
	double ScatteringLengthOfDisplacedSolvent = CurrentResidue.Volume * ScatteringLengthDensityOfSolvent; // in case of WAT and modifications X scaling is applied to volume
	double ExcessScatteringLength ;

	if ( Contrast < 0.0 )
	{
		ScatteringLengthOfResidue = CurrentResidue.XRayScatteringLength ;
		ExcessScatteringLength    = ScatteringLengthOfResidue - ScatteringLengthOfDisplacedSolvent ;

		//x = CurrentResidue.xVolume ;
		x = (CurrentResidue.xXRayScattering * ScatteringLengthOfResidue - CurrentResidue.xVolume * ScatteringLengthOfDisplacedSolvent) / ExcessScatteringLength;
		//y = CurrentResidue.yVolume ;
		y = (CurrentResidue.yXRayScattering * ScatteringLengthOfResidue - CurrentResidue.yVolume * ScatteringLengthOfDisplacedSolvent) / ExcessScatteringLength;
		// z = CurrentResidue.zVolume ;
		z = (CurrentResidue.zXRayScattering * ScatteringLengthOfResidue - CurrentResidue.zVolume * ScatteringLengthOfDisplacedSolvent) / ExcessScatteringLength;
	}
	else
	{
		ScatteringLengthOfResidue = CurrentResidue.NeutronScatteringLength;
		ExcessScatteringLength    = ScatteringLengthOfResidue - ScatteringLengthOfDisplacedSolvent;

		//x = CurrentResidue.xVolume ; 
		x = (CurrentResidue.xNeutronScattering * ScatteringLengthOfResidue - CurrentResidue.xVolume * ScatteringLengthOfDisplacedSolvent) / ExcessScatteringLength;
		//y = CurrentResidue.yVolume ; 
		y = (CurrentResidue.yNeutronScattering * ScatteringLengthOfResidue - CurrentResidue.yVolume * ScatteringLengthOfDisplacedSolvent) / ExcessScatteringLength;
		// z = CurrentResidue.zVolume ;
		z = (CurrentResidue.zNeutronScattering * ScatteringLengthOfResidue - CurrentResidue.zVolume * ScatteringLengthOfDisplacedSolvent) / ExcessScatteringLength;
	}
	//printf("%s, %0.3e, %f ,%f, %f\n",ScatteringLengthOfResidue,CurrentResidue.Name,x,y,z);
	Radius = sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2)) ;
	Theta  = acos(z / Radius) ;
	Phi    = acos(x / (Radius * sin(Theta))) * Sign(y) ;
	//printf("Rad=%lf The=%lf Phi=%lf ", Radius, Theta, Phi) ;




	/* GSL < 2.0 */
	/*
	double Legendre[NumberOfHarmonics + 1] ;
	double Bessel[NumberOfHarmonics + 1] ;

	// Calculate spherical Bessel functions for l = 0 to NumberOfHarmonics
	gsl_sf_bessel_jl_array( NumberOfHarmonics, q * Radius, Bessel) ;

	// Calculate all associated Legendre polynomials $P_l^m(cos(Theta))$ of degree l = m ... NumberOfHarmonics including prefactor $\sqrt{(2l+1)/(4\pi)} \sqrt{(l-m)!/(l+m)!}$ - store the values in Legendre[m], Legendre[m + 1], ..., Legendre[NumberOfHarmonics]
	for ( m = 0; m < NumberOfHarmonics + 1; m++)
	{
		// http://home.thep.lu.se/~jari/documents/gsl-ref.html/Associated-Legendre-Polynomials-and-Spherical-Harmonics.html
		// int gsl_sf_legendre_sphPlm_array (int lmax, int m, double x, double result_array[]) 	Function
		// These functions compute an array of spherical harmonic associated Legendre polynomials $\sqrt{(2l+1)/(4\pi)} \sqrt{(l-m)!/(l+m)!} P_l^m(x)$, and optionally their derivatives, for m >= 0, l = |m|, ..., lmax, |x| <= 1.0
		// i.e. they provide the unnormalized associated Legendre polynomials P_l^m(x) times the prefactor !!!
		// \sqrt{(2l+1)/(4\pi)} \sqrt{(l-m)!/(l+m)!} P_l^m( cos(theta) ) for m >= 0 and l = m ... NumberOfHarmonics
		// fills Legendre[m,...,NumberOfHarmonics] from &Legendre[m] on
		gsl_sf_legendre_sphPlm_array( NumberOfHarmonics, m, cos(Theta), &Legendre[m]) ;

		for ( l = m; l < NumberOfHarmonics + 1; l++)
		{
			printf("l=%d m=%d Bessel[l]=%g Legendre[l]=%lf ", l, m, Bessel[l], Legendre[l]) ;
			printf("dBeta=%g\n", sqrt(DeltaB)*sqrt(4.0 * M_PI) * cpow(I, l) * ExcessScatteringLength * Bessel[l] * Legendre[l] * PolarComplexNumber(1.0, -m * Phi)) ;

			Beta[l][m] += sqrt(DeltaB)*sqrt(4.0 * M_PI) * cpow(I, l) * ExcessScatteringLength * Bessel[l] * Legendre[l] * PolarComplexNumber(1.0, -m * Phi) ;
		}
	}
	*/




	/* GSL >= 2.0 */

	// https://www.gnu.org/software/gsl/doc/html/specfunc.html?highlight=gsl_sf_bessel_jl_array#associated-legendre-polynomials-and-spherical-harmonics

	/*
	NEW for GSL >= 2.0
	size_t gsl_sf_legendre_array_n(const size_t lmax)
	This function returns the minimum array size for maximum degree lmax needed for the array versions of the associated Legendre functions. 
	Size is calculated as the total number of P_l^m(x) functions, plus extra space for precomputing multiplicative factors used in the recurrence relations.
	*/
	size_t LegendreSize = gsl_sf_legendre_array_n( NumberOfHarmonics ) ;
	double Legendre[LegendreSize] ;
	double Bessel[NumberOfHarmonics + 1] ;
	//printf("LegendreSize=%d\n", LegendreSize) ;


	// Calculate spherical Bessel functions for l = 0 to NumberOfHarmonics
	/*
	Same as for GSL < 2.0
	int gsl_sf_bessel_jl_array(int lmax, double x, double result_array[])
	This routine computes the values of the regular spherical Bessel functions j_l(x) for l from 0 to lmax inclusive for lmax >= 0 and x >= 0, storing the results in the array result_array. The values are computed using recurrence relations for efficiency, and therefore may differ slightly from the exact values.
	*/
	gsl_sf_bessel_jl_array( NumberOfHarmonics, q * Radius, Bessel) ;

	// Calculate all associated Legendre polynomials $P_l^m(cos(Theta))$ of degree l = m ... NumberOfHarmonics including prefactor $\sqrt{(2l+1)/(4\pi)} \sqrt{(l-m)!/(l+m)!}$
	/*
	NEW for GSL >= 2.0
	int gsl_sf_legendre_array(const gsl_sf_legendre_t norm, const size_t lmax, const double x, double result_array[])
	These functions calculate all normalized associated Legendre polynomials for 0 <= l <= lmax and 0 <= m <= l for |x| <= 1.
	The norm parameter specifies which normalization is used.
	The normalized P_l^m(x) values are stored in result_array, whose minimum size can be obtained from calling gsl_sf_legendre_array_n().
	The array index of P_l^m(x) is obtained from calling gsl_sf_legendre_array_index(l, m).
	To include or exclude the Condon-Shortley phase factor of (-1)^m, set the parameter csphase to either -1 or 1 respectively in the _e function.
	This factor is included by default.

	gsl_sf_legendre_t:
		Value				Description
		GSL_SF_LEGENDRE_NONE		The unnormalized associated Legendre polynomials P_l^m(x)
		GSL_SF_LEGENDRE_SCHMIDT		The Schmidt semi-normalized associated Legendre polynomials S_l^m(x)
	--->	GSL_SF_LEGENDRE_SPHARM		The spherical harmonic associated Legendre polynomials Y_l^m(x) 	<---
		GSL_SF_LEGENDRE_FULL		The fully normalized associated Legendre polynomials N_l^m(x)
	*/
	const gsl_sf_legendre_t LegendreNorm = GSL_SF_LEGENDRE_SPHARM ; // provides $\sqrt{(2l+1)/(4\pi)} \sqrt{(l-m)!/(l+m)!} P_l^m(x)$ if the Condon-Shortley phase factor (-1)^m is omitted via csphase = 1.0
	const double csphase = 1.0 ;

	gsl_sf_legendre_array_e( LegendreNorm, NumberOfHarmonics, cos(Theta), csphase, Legendre) ;

	size_t LegendreIndex ;
	for ( m = 0; m < NumberOfHarmonics + 1; m++)
	{
		for ( l = m; l < NumberOfHarmonics + 1; l++)
		{
			/*
			NEW for GSL >= 2.0
			size_t gsl_sf_legendre_array_index(const size_t l, const size_t m)
			This function returns the index into result_array, result_deriv_array, or result_deriv2_array corresponding to P_l^m(x), P_l^{'m}(x), or P_l^{''m}(x).
			The index is given by l(l+1)/2 + m.
			*/
			LegendreIndex = gsl_sf_legendre_array_index( l, m) ;
			//printf("l=%d m=%d LegendreIndex=%d Bessel[l]=%g Legendre[LegendreIndex]=%lf ", l, m, LegendreIndex, Bessel[l], Legendre[LegendreIndex]) ;
			//printf("dBeta=%g\n", sqrt(DeltaB)*sqrt(4.0 * M_PI) * cpow(I, l) * ExcessScatteringLength * Bessel[l] * Legendre[LegendreIndex] * PolarComplexNumber(1.0, -m * Phi)) ;
			Beta[l][m] += sqrt(DeltaB)*sqrt(4.0 * M_PI) * cpow(I, l) * ExcessScatteringLength * Bessel[l] * Legendre[LegendreIndex] * PolarComplexNumber(1.0, -m * Phi) ;
		}
	}


}


void CopyResidue(struct Residue * Original, struct Residue * Copy)
{
	Copy->xVolume = Original->xVolume;
	Copy->yVolume = Original->yVolume;
	Copy->zVolume = Original->zVolume;

	Copy->xXRayScattering = Original->xXRayScattering;
	Copy->yXRayScattering = Original->yXRayScattering;
	Copy->zXRayScattering = Original->zXRayScattering;

	Copy->xNeutronScattering = Original->xNeutronScattering;
	Copy->yNeutronScattering = Original->yNeutronScattering;
	Copy->zNeutronScattering = Original->zNeutronScattering;

	Copy->XRayScatteringLength    = Original->XRayScatteringLength;
	Copy->NeutronScatteringLength = Original->NeutronScatteringLength;

	Copy->Volume    = Original->Volume;
	Copy->Weight    = Original->Weight ; // added by MS

	Copy->Name[0]   = Original->Name[0];
	Copy->Name[1]   = Original->Name[1];
	Copy->Name[2]   = Original->Name[2];
	Copy->Name[3]   = Original->Name[3];

	Copy->AtomName[0]   = Original->AtomName[0]; // added by MS
	Copy->AtomName[1]   = Original->AtomName[1]; // added by MS
	Copy->AtomName[2]   = Original->AtomName[2]; // added by MS

	Copy->ResidueID = Original->ResidueID;
}


double complex pol(double r, double phi)
{
	if ( phi == 0 )
		return r+I*0 ;
	else
		return r*(cos(phi)+I*sin(phi)) ;
}

double complex **ComplexArray( int dim1, int dim2)
{
	double complex ** arr;
	arr = (double complex**)malloc(dim1*sizeof(double complex*)) ;
	for (int i = 0; i< dim1; i++)
	{
		arr[i] = (double complex*) malloc((dim2)*sizeof(double complex)) ;
		for (int ii = 0; ii < dim2; ii++) { arr[i][ii] = 0 ; }
	}
	return arr ;
}

void FreeComplexArray(double complex **alpha, int dim1, int dim2)
{
	for(int i=0;i<dim1;i++) { free(alpha[i]) ; }
	free(alpha) ;
}


