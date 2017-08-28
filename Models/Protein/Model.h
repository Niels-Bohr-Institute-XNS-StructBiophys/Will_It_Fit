double Model(double q, double * Parameters, double * Constraints, double Contrast, struct Protein ProteinStructure, struct UserDefined * UserDefinedStructure)
{
    /// Declarations
    // Dummy integers
    int j;
    int l;
    int m;
    // Variables describing the sample
    double Intensity;
    double IntensitySolvent;
    double I0Prot;
    double I0;
    double Scaling;
    double Background;
    double ConcentrationOfSample          =  Constraints[Conc];
    struct Residue CurrentResidue;
    double ScatteringLengthDensityOfWater =  Constraints[SLDWater];
    double complex ** Beta = ComplexArray(NumberOfHarmonics + 1, NumberOfHarmonics + 1);
    double complex ** BetaSolvent = ComplexArray(NumberOfHarmonics + 1, NumberOfHarmonics + 1);


	//Loop over residues
      for(j = 0; j < ProteinStructure.NumberOfResidues; j++){ // for AA based calculation
	
	     
    	      CopyResidue(&ProteinStructure.Residues[j], &CurrentResidue); // For AA based calculation
	      
// 	      if (q ==0.5 ){ //Print statement for housekeeping
//                            printf("current residue: %s Volume %0.2e %f,%f,%f\n",CurrentResidue.Name, CurrentResidue.Volume, CurrentResidue.xVolume, CurrentResidue.yVolume,CurrentResidue.zVolume);
//			   }
             
	      if (strcmp(CurrentResidue.Name, "WAT") == 0 ) {
                 				CurrentResidue.Volume = CurrentResidue.Volume*(1./(1.+Parameters[HYDR]));
						if (Contrast >= 0.0 && Contrast <= 100.0) {
						
						CurrentResidue.NeutronScatteringLength = 4.133 *(2*((Contrast* 6.671e-13 / 100.) +  ((100.-Contrast)* -3.741e-13 / 100.)) + 5.803); 
						}
	      }

              if (strcmp(CurrentResidue.Name,  "  X") == 0 ) {
    	    					CurrentResidue.Volume = CurrentResidue.Volume*(1./(1.+Parameters[GLYCV]));
          					}
	     AddScatteringFromResidue(Beta, q, CurrentResidue, Contrast, ScatteringLengthDensityOfWater, 1.0); // Parameters[PROTSCALE]);
	     //AddScatteringFromSolvent(BetaSolvent, q, CurrentResidue, Contrast, ScatteringLengthDensityOfWater, 1.0); // Parameters[PROTSCALE]);
}


    // Calculate intensity
    Intensity = 0.0;
    IntensitySolvent = 0.0;
    I0Prot = 0.0;
    I0 = 0.0;

    for(l = 0; l < NumberOfHarmonics + 1; l++) {

      for (m = 0; m <= l; m++) {
        Intensity += ((m > 0) + 1) * pow(cabs(Beta[l][m]), 2); //Calculate intensity
        //IntensitySolvent += ((m > 0) + 1) * pow(cabs(BetaSolvent[l][m]), 2); //Calculate intensity

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
    Intensity = Scaling*(Intensity)*ConcentrationOfSample; //;/pow(2.818E-13,2)
    Intensity = Intensity- Background;     
    return Intensity;
}
