double Model(double q, double * Parameters, double * Constraints, double Contrast, struct Protein ProteinStructure, struct UserDefined * UserDefinedStructure)
{

    // Get scattering length densities
    double ScatteringLengthDensityOfCaps = Constraints[0];
    double ScatteringLengthDensityOfCore = Constraints[1];
    double ScatteringLengthDensityOfMethyl = Constraints[2];
    double ScatteringLengthDensityOfBelt = Constraints[3];
    double ScatteringLengthDensityOfWater = Constraints[4];
    
    // Get the parameters describing the geometry of the disc
    double HeightOfBelt = Parameters[2];
    double HeightOfLipids = Constraints[6];
    double HeightOfCore = Constraints[7];
    double HeightOfMethyl = Constraints[8];
    double VerticalAxisOfEllipsoid = Constraints[9];
    double ScaleFactorOfCaps = Constraints[10];
    double MajorSemiAxisOfCore = Constraints[23];
    double MinorSemiAxisOfCore = Constraints[24];
    double ThicknessOfBelt = Constraints[27];

    // Get the parameters describng the membrane protein
    double CorrectionToMPVolume = Parameters[11];
    double MPTranslation[3];
    double MPRotation[3][3];
    SetTranslationVector(MPTranslation, Parameters[24], Parameters[25], Parameters[26]);
    SetRotationMatrix(MPRotation, Parameters[27], Parameters[28], Parameters[29]);

    // Parameters describing the solution
    double Dummy;
    double Roughness;
    double Scaling;
    double Background;
    double Intensity;
    double ConcentrationOfSample = Constraints[21];

    // Compute the roughness from the parameters
    if (Contrast >= 0.0 && Contrast <= 100.0) {
        Dummy      = q * Parameters[18];
        Roughness  = exp(-(Dummy * Dummy));
        Scaling    = Parameters[20] * (Contrast / 100.0 * Parameters[13] + (100.0 - Contrast) / 100.0 * Parameters[15]);
        Background = Contrast / 100.0 * Parameters[12] + (100.0 - Contrast) / 100.0 * Parameters[14];
    } else {
        Dummy      = q * Parameters[8];
        Roughness  = exp(-(Dummy * Dummy));
        Scaling    = Parameters[17] * Parameters[20];
        Background = Parameters[16];
    }
    
    int j,l,m;
    double complex ** beta=ComplexArray(Nh+1,Nh+1); //Array for holding harmonic coefficients for bead based model of membrane protein
    double complex ** alpha=ComplexArray(Nh+1,Nh+1); // Array for holding harmonic coefficients for analytical model of nanodisc and tags
    double bg;
    struct Residue CurrentResidue;

    /******************Calculate alpha*************************/
    NanodiscPDBModel(alpha,q, ScatteringLengthDensityOfCaps, ScatteringLengthDensityOfCore, ScatteringLengthDensityOfMethyl,
                     ScatteringLengthDensityOfBelt, ScatteringLengthDensityOfWater,
                     MajorSemiAxisOfCore, MinorSemiAxisOfCore, ThicknessOfBelt, 
                     HeightOfLipids, HeightOfCore, HeightOfMethyl, HeightOfBelt, 
                     VerticalAxisOfEllipsoid, ScaleFactorOfCaps);

    /******************Calculate beta*************************/
    for(j=0; j<ProteinStructure.NumberOfResidues; j++){
        CopyResidue(&ProteinStructure.Residues[j],&CurrentResidue);

        Orient(&CurrentResidue.xVolume, &CurrentResidue.yVolume, &CurrentResidue.zVolume, MPRotation, MPTranslation);
        if(Contrast<0){
            Orient(&CurrentResidue.xXRayScattering, &CurrentResidue.yXRayScattering, &CurrentResidue.zXRayScattering, MPRotation, MPTranslation);
        }
        else{
            Orient(&CurrentResidue.xNeutronScattering, &CurrentResidue.yNeutronScattering, &CurrentResidue.zNeutronScattering, MPRotation, MPTranslation);
        }

       AddScatteringFromResidue(beta, q, CurrentResidue, Contrast,
           ScatteringLengthDensityOfCaps, ScatteringLengthDensityOfCore, ScatteringLengthDensityOfMethyl,
           ScatteringLengthDensityOfBelt, ScatteringLengthDensityOfWater,
           HeightOfLipids, HeightOfCore, HeightOfMethyl, HeightOfBelt, CorrectionToMPVolume,
           MajorSemiAxisOfCore, MinorSemiAxisOfCore, 
           VerticalAxisOfEllipsoid, ScaleFactorOfCaps);
    }
   
    /***********Calculate combined intensity from alpha and beta**********/

    Intensity=0;
    for(l=0;l<Nh+1;l++){
        for(m=0;m<=l;m++){
            Intensity += ((m>0)+1)*pow( cabs(sqrt(Roughness)*alpha[l][m]+beta[l][m]),2 );
        }
    }
    
    //****************Debugging**********************************
    //if(isnan(Intensity)){
    //    printf("Warning Intensity is Not A Number:%g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g\n", 
    //            q, Roughness,
    //            ScatteringLengthDensityOfCaps, ScatteringLengthDensityOfCore, ScatteringLengthDensityOfMethyl,
    //            ScatteringLengthDensityOfBelt, ScatteringLengthDensityOfWater, 
    //            HeightOfLipids, HeightOfCore, HeightOfMethyl, HeightOfBelt, 
    //            CorrectionToMPVolume,
    //            MajorSemiAxisOfCore, MinorSemiAxisOfCore,
    //            VerticalAxisOfEllipsoid, ScaleFactorOfCaps);
    //}
    //**********************************************************

    // Scale the sum
    Intensity = Intensity * ConcentrationOfSample;
    Intensity = Intensity * Scaling - Background;

    // Free the arrays
    FreeComplexArray(alpha, Nh+1, Nh+1);
    FreeComplexArray(beta, Nh+1, Nh+1);

    return Intensity;
}
