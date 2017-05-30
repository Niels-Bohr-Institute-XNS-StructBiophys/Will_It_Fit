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
    double complex SFBeta;
    double Scaling;
    double Background;
    double ConcentrationOfSample          =  Constraints[Conc];
    struct Residue CurrentResidue;
    double ScatteringLengthDensityOfWater =  Constraints[SLDWater];
    double complex ** Beta = ComplexArray(NumberOfHarmonics + 1, NumberOfHarmonics + 1);
    double complex ** BetaProt = ComplexArray(NumberOfHarmonics + 1, NumberOfHarmonics + 1);
    double Q0 = 0.00001;
    double I0Agg;
    double Sagg;
    double Dm = Parameters[FracDim];
    double RgAgg = Parameters[AggRg];
//    double r0 = 2.*Constraints[RZero];
    double r0 = Parameters[RZeroPar];
    double XiAgg = sqrt(2 * pow(RgAgg, 2) / (Dm * (Dm + 1)));



      for(j = 0; j < ProteinStructure.NumberOfResidues; j++){
          CopyResidue(&ProteinStructure.Residues[j], &CurrentResidue);
//        printf("current residue: %c Volume %0.2e \n", *CurrentResidue.Name, CurrentResidue.Volume);
          if (*CurrentResidue.Name == 'W') {
          CurrentResidue.Volume = CurrentResidue.Volume/Parameters[HYDR];

    }
    if (*CurrentResidue.Name == 'X') {
    CurrentResidue.Volume = CurrentResidue.Volume/Parameters[GLYCV];
}
		      AddScatteringFromResidue(Beta, q, CurrentResidue, Contrast, ScatteringLengthDensityOfWater, Parameters[PROTSCALE]);
          AddScatteringFromResidue(BetaProt, Q0, CurrentResidue, Contrast, ScatteringLengthDensityOfWater, Parameters[PROTSCALE]);

      }


    int NL=1000;
    double L = Parameters[ContourLength];
    double b = Parameters[KuhnLength];


    // Calculate intensity
    Intensity = 0.0;
    I0Prot = 0.0;
    I0 = 0.0;
    SFBeta = 0.0 +I*0.0;

    for(l = 0; l < NumberOfHarmonics + 1; l++) {

      for (m = 0; m <= l; m++) {
        Intensity += ((m > 0) + 1) * pow(cabs(Beta[l][m]), 2); //Calculate intensity
        I0Prot += ((m > 0) + 1) * pow(cabs(BetaProt[l][m]), 2); //Calculate intensity at Q0 for normalization
        SFBeta += ((m > 0) + 1) * Beta[l][m]; //Calculate complex amplitude

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
    Sagg = 0.0;
    Sagg = Dm*tgamma(Dm-1) * sin((Dm-1) * atan(q*XiAgg));
    Sagg /= (pow(q*r0, Dm)*pow(1 + 1/pow(q*XiAgg,2), (Dm-1)/2));
    Sagg += 1.;
    I0 = ConcentrationOfSample*(I0Prot) * Scaling;
    Intensity = Scaling*Intensity/I0Prot;
    Intensity = Intensity  - Background;     return Intensity; //;
}
