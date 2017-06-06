double Model(double q, double * Parameters, double * Constraints, double Contrast, struct Protein ProteinStructure, struct UserDefined * UserDefinedStructure)
{
    /// Declarations
    // Dummy integers
    int j;
    int l;
    int m;
    // Variables describing the sample
    double Intensity;
    double I0Prot;
    double I0;
    double Scaling;
    double Background;
    double ConcentrationOfSample          =  Constraints[Conc];
    struct Residue CurrentResidue;
    double ScatteringLengthDensityOfWater =  Constraints[SLDWater];
    double complex ** Beta = ComplexArray(NumberOfHarmonics + 1, NumberOfHarmonics + 1);


	//Loop over residues
      for(j = 0; j < ProteinStructure.NumberOfResidues; j++){ // for AA based calculation
	
	     
    	      CopyResidue(&ProteinStructure.Residues[j], &CurrentResidue); // For AA based calculation
// 	      if (q ==0.5 ){ //Print statement for housekeeping
//                            printf("current residue: %c Volume %0.2e \n", *CurrentResidue.Name, CurrentResidue.Volume);
//			   }
             
	      if (*CurrentResidue.Name == 'W') {
                 				CurrentResidue.Volume = CurrentResidue.Volume/Parameters[HYDR];
           				      	}

              if (*CurrentResidue.Name == 'X') {
    	    					CurrentResidue.Volume = CurrentResidue.Volume/Parameters[GLYCV];
          					}
	     AddScatteringFromResidue(Beta, q, CurrentResidue, Contrast, ScatteringLengthDensityOfWater, 1.0); // Parameters[PROTSCALE]);
	}


    // Calculate intensity
    Intensity = 0.0;
    I0Prot = 0.0;
    I0 = 0.0;

    for(l = 0; l < NumberOfHarmonics + 1; l++) {

      for (m = 0; m <= l; m++) {
        Intensity += ((m > 0) + 1) * pow(cabs(Beta[l][m]), 2); //Calculate intensity

      }
    }



   /// Rescale result and return
    if (Contrast >= 0.0 && Contrast <= 100.0) {
        Scaling    = Parameters[SCALECONC] * (Contrast / 100.0 * Parameters[SCALEN100] + (100.0 - Contrast) / 100.0 * Parameters[SCALEN0]);
        Background = Contrast / 100.0 * Parameters[BACKN100] + (100.0 - Contrast) / 100.0 * Parameters[BACKN0];
    } else {
        Scaling    = Parameters[SCALECONC] * Parameters[SCALEX];
        Background = Parameters[BACKX];
    }
    I0 = ConcentrationOfSample*(I0Prot) * Scaling;
    Intensity = Scaling*Intensity/pow(2.82e-13,2);
    Intensity = Intensity  - Background;     return Intensity; //;
}
