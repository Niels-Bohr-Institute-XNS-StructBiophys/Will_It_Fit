double Model(double q, double * Parameters, double * Constraints, double Contrast, struct Protein ProteinStructure, struct UserDefined * UserDefinedStructure) 
{
	// Dummy integers    
    int j;
    int l;
    int m;

	// Get scattering length densities
    double ScatteringLengthDensityOfWater = Constraints[0];
    double ConcentrationOfSample          = Constraints[1];

    // Parameters describing the solution
    double Dummy;
    double Roughness;
    double Scaling;
    double Background;
    double Intensity;

	// Parameters describing the protein structure
    double complex ** Beta = ComplexArray(NumberOfHarmonics + 1, NumberOfHarmonics + 1);
    struct Residue CurrentResidue; 

    // Compute the roughness from the parameters
    if (Contrast >= 0.0 && Contrast <= 100.0) {
        Dummy      = q * Parameters[1];
        Roughness  = exp(-(Dummy * Dummy));
        Scaling    = Parameters[8] * (Contrast / 100.0 * Parameters[3] + (100.0 - Contrast) / 100.0 * Parameters[5]);
        Background = Contrast / 100.0 * Parameters[2] + (100.0 - Contrast) / 100.0 * Parameters[4];
    } else {
        Dummy      = q * Parameters[0];
        Roughness  = exp(-(Dummy * Dummy));
        Scaling    = Parameters[7] * Parameters[8];
        Background = Parameters[6];
    }
    //printf("Contrast computed \n");
    // Construct complex matrix
    for(j = 0; j < ProteinStructure.NumberOfResidues; j++){
        //printf("Residue nr: %d \n", j);
        CopyResidue(&ProteinStructure.Residues[j], &CurrentResidue);
        //printf("Residues copied\n");
		AddScatteringFromResidue(Beta, q, CurrentResidue, Contrast, ScatteringLengthDensityOfWater);
        //printf("Scattering added\n");
    }
   
    // Calculate intensity
    Intensity = 0.0;

    for(l = 0; l < NumberOfHarmonics + 1; l++) {

        for (m = 0; m <= l; m++) {
            Intensity += ((m > 0) + 1) * pow(cabs(sqrt(Roughness) * Beta[l][m]), 2);
        }
    }

    // Scale the sum
    Intensity = Intensity * ConcentrationOfSample * Scaling - Background;

    // Free the arrays
    FreeComplexArray(Beta, NumberOfHarmonics + 1, NumberOfHarmonics + 1);

    return Intensity;
}
