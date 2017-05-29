void ComputeConstraints(double * Parameters, double * VolumesOfMolecules, double * ScatteringLengths, double Contrast,
                        double Concentration, double * Constraints, struct Protein ProteinStructure, struct UserDefined * UserDefinedStructure)
{
    // Declare dummy variables needed in function
    double Dummy1;
    double Dummy2;
    double AreaOfDisc;
    int j,ss=0;;

    // Variables describing the protein orientation
    double VA,VM,VH,VT;
    struct Residue CurrentResidue;
    double translationvector[3];
    double rotationmatrix[3][3];
    SetTranslationVector(translationvector, Parameters[24], Parameters[25], Parameters[26]);
    SetRotationMatrix(rotationmatrix, Parameters[27], Parameters[28], Parameters[29]);

    // Variables describing the lipid core
    double HeightOfCore;
    double HeightOfBelt;
    double HeightOfLipids;
    double HeightOfMethyl;

    // Variables describing scattering lengths
    double ScatteringLengthDensityOfCaps;
    double ScatteringLengthDensityOfCore;
    double ScatteringLengthDensityOfBelt;
    double ScatteringLengthDensityOfWater;
    double ScatteringLengthDensityOfMethyl;

    double ScatteringLengthOfBelt;
    double ScatteringLengthOfCore;
    double ScatteringLengthOfWater;
    double ScatteringLengthOfCaps;
    double ScatteringLengthOfMethyl;

    // Variables describing the volumes
    double VolumeOfWater;
    double VolumeOfHead;
    double VolumeOfHeads;
    double VolumeOfTail;
    double VolumeOfTails;
    double VolumeOfMethyl;
    double VolumeOfMethyls;
    double VolumeOfBelt;

    double CorrectionToVolumeOfLipid;
    double CorrectionToVolumeOfWater;
    double CorrectionToVolumeOfBelt;
    double CorrectionToVolumeOfMP;

    // Variables describing water
    double NumberOfWaterAtBelt;
    double NumberOfWaterAtHead;

    // Variables describing the properties of the disc
    double RatioBetweenAxis;
    double ThicknessOfBelt;

    double SemiMajorAxisOfCore;
    double SemiMinorAxisOfCore;

    // Variables describing the lipids
    double NumberOfLipids;
    double DisplacedAlkyl0;
    double DisplacedAlkyl=14;
    double DisplacedMethyl0;
    double DisplacedMethyl=14;
    double DisplacedHeads=14;

    // Parameters describing the endcaps
    double VerticalAxisOfEllipsoid;
    double VerticalShiftOfEllipsoidCenter;
    double ScaleFactorOfCaps;
    double SemiMajorAxisOfEndcaps;
    double SemiMinorAxisOfEndcaps;
    double MaxThicknessOfCaps;

    /// Get parameters from Parameters
    // Scale factors
    CorrectionToVolumeOfBelt = Parameters[9];
    CorrectionToVolumeOfLipid = Parameters[10];
    CorrectionToVolumeOfWater = Parameters[19];
    CorrectionToVolumeOfMP   = Parameters[11];

    // Hydration numbers
    NumberOfWaterAtHead = fabs(Parameters[5]);
    NumberOfWaterAtBelt = Parameters[6];

    /// Get parameters from variable VolumesOfMolecules
    // Volume of water
    VolumeOfWater = VolumesOfMolecules[0];
    // Volume of the hydrophilic head groups - including hydration water
    VolumeOfHead = VolumesOfMolecules[1] * CorrectionToVolumeOfLipid + NumberOfWaterAtHead * VolumesOfMolecules[0] * CorrectionToVolumeOfWater;
    // Volume of lipid tail group without the methyl endgroups
    VolumeOfTail = VolumesOfMolecules[2] * CorrectionToVolumeOfLipid;
    // Volue of methyl end groups
    VolumeOfMethyl = VolumesOfMolecules[3] * CorrectionToVolumeOfLipid;
    // Volume of two belts, i.e. the total belt volume
    VolumeOfBelt = 2.0 * VolumesOfMolecules[4] * CorrectionToVolumeOfBelt + 2.0 * NumberOfWaterAtBelt * VolumesOfMolecules[0];

    // Obtain scattering length densities
    // Get the scattering length of hydration water!
    ScatteringLengthOfWater = ScatteringLengths[0];
    // Calculate the scattering length of the hydrophilic head groups, INCLUDING hydration water!
    ScatteringLengthOfCaps = ScatteringLengths[1] + NumberOfWaterAtHead * ScatteringLengths[0];
    // Get the scattering length of the hydrophobic tails minus methyl groups
    ScatteringLengthOfCore = ScatteringLengths[2];
    // Get the scattering length of the hydrophobic tails minus methyl groups
    ScatteringLengthOfMethyl = ScatteringLengths[3];
    // Scattering length of two belts, i.e. the total belt volume
    ScatteringLengthOfBelt = 2.0f * (ScatteringLengths[4] + NumberOfWaterAtBelt * ScatteringLengths[0]);

    /// Derive scattering lengths
    // Calculate the excess scattering length density of the hydration water
    ScatteringLengthDensityOfWater = ScatteringLengthOfWater / VolumeOfWater;
    // Calculate the excess scattering length density of the head
    ScatteringLengthDensityOfCaps = ScatteringLengthOfCaps / VolumeOfHead;
    // Calculate the excess scattering length density of the hydrophobic core
    ScatteringLengthDensityOfCore = ScatteringLengthOfCore / VolumeOfTail;
    // Calculate the excess scattering length density of the methyl
    ScatteringLengthDensityOfMethyl = ScatteringLengthOfMethyl / VolumeOfMethyl;
    // Calculate the excess scattering length density of the belt
    ScatteringLengthDensityOfBelt = ScatteringLengthOfBelt / VolumeOfBelt;
        
    // Assign values of the constraints
    // Values to use in computing model
    Constraints[0] = ScatteringLengthDensityOfCaps;
    Constraints[1] = ScatteringLengthDensityOfCore;
    Constraints[2] = ScatteringLengthDensityOfMethyl;
    Constraints[3] = ScatteringLengthDensityOfBelt;
    Constraints[4] = ScatteringLengthDensityOfWater;
    Constraints[5] = 0.0;

    // Use fabs to avoid problems with negative squareroots later in the program
    RatioBetweenAxis = fabs(Parameters[0]);
    // Height of belt
    HeightOfBelt = Parameters[2];
    HeightOfCore = Parameters[2] + Parameters[1];

    // Number of lipids
    NumberOfLipids = Parameters[3];
    //HeightOfLipidCylinder = Parameters[1];
    //AreaPrHeadgroup = Parameters[1];
    MaxThicknessOfCaps = fabs(Parameters[21]);
    ScaleFactorOfCaps = Parameters[22];

    Dummy1 = sqrt(pow(ScaleFactorOfCaps,2)-1.)/ScaleFactorOfCaps;
    VerticalAxisOfEllipsoid = MaxThicknessOfCaps/(1.-Dummy1);
    VerticalShiftOfEllipsoidCenter = Dummy1*MaxThicknessOfCaps/(1.-Dummy1);

    Dummy2=pow(ScaleFactorOfCaps,2)*VerticalAxisOfEllipsoid * (2./3.-Dummy1+1./3.*pow(Dummy1,3));

    do{
        DisplacedAlkyl0=DisplacedAlkyl;
        DisplacedMethyl0=DisplacedMethyl;

        // Volume of the lipid layers including the parts displaced by the membrane protein
        VolumeOfHeads = VolumeOfHead * (NumberOfLipids+DisplacedHeads);
        VolumeOfTails = VolumeOfTail * (NumberOfLipids+DisplacedAlkyl); 
        VolumeOfMethyls = VolumeOfMethyl * (NumberOfLipids+DisplacedMethyl);

        HeightOfMethyl = (HeightOfCore+2*Dummy2)/(VolumeOfTails/VolumeOfMethyls + 1);
        if (HeightOfMethyl < 0.0) {
            printf("Warning: Height of methyl layer is below zero!!! \n");
        }
        SemiMinorAxisOfCore = sqrt(VolumeOfMethyls/(HeightOfMethyl*pi*RatioBetweenAxis));
        SemiMajorAxisOfCore = SemiMinorAxisOfCore * RatioBetweenAxis;
        AreaOfDisc=pi*SemiMinorAxisOfCore*SemiMajorAxisOfCore;
        SemiMajorAxisOfEndcaps = ScaleFactorOfCaps * SemiMajorAxisOfCore;
        SemiMinorAxisOfEndcaps = ScaleFactorOfCaps * SemiMinorAxisOfCore;
        HeightOfLipids = VolumeOfHeads / AreaOfDisc + HeightOfCore;

        VA=0;VM=0;VH=0;VT=0;
        for(j=0;j<ProteinStructure.NumberOfResidues;j++){
            CopyResidue(&ProteinStructure.Residues[j], &CurrentResidue);
            Orient(&CurrentResidue.xVolume, &CurrentResidue.yVolume, &CurrentResidue.zVolume, rotationmatrix, translationvector);
            if(ResidueIsInMethylLayer(CurrentResidue.zVolume,HeightOfMethyl))
                VM+=CurrentResidue.Volume; // Residue is in methyl layer
            else if(ResidueIsInAlkylLayer(CurrentResidue.xVolume, CurrentResidue.yVolume, CurrentResidue.zVolume, SemiMajorAxisOfEndcaps, SemiMinorAxisOfEndcaps, VerticalAxisOfEllipsoid, -VerticalShiftOfEllipsoidCenter, HeightOfCore))
                VA+=CurrentResidue.Volume; // Residue is in alkyl layer
            else if(ResidueIsInLipidLayer(CurrentResidue.xVolume, CurrentResidue.yVolume, CurrentResidue.zVolume, SemiMajorAxisOfEndcaps, SemiMinorAxisOfEndcaps, VerticalAxisOfEllipsoid, -VerticalShiftOfEllipsoidCenter, HeightOfLipids))
                VH+=CurrentResidue.Volume;  // Residue is in head group layer
            VT+=CurrentResidue.Volume;
        }
        DisplacedHeads=VH/VolumeOfHead;
        DisplacedAlkyl=VA/VolumeOfTail;
        DisplacedMethyl=VM/VolumeOfMethyl;
        //printf("DisplacedHeads=%.4lf,DisplacedAlkyl=%.4lf, DispalcedMethyl=%.4lf, Total Vol=%.4lf\n",DisplacedHeads,DisplacedAlkyl,DisplacedMethyl,VT);
        //printf("DisplacedHeads=%.4lf,DisplacedAlkyl=%.4lf, DispalcedMethyl=%.4lf, Total Vol=%.4lf\n",VH/VT,VA/VT,VM/VT,VT);
        ss++;

    }while(fabs(DisplacedMethyl-DisplacedMethyl0)/DisplacedMethyl + fabs(DisplacedAlkyl-DisplacedAlkyl0)/DisplacedAlkyl > 0.1  && ss<10);

    // Volume of the lipid layers including the parts displaced by the membrane protein
    VolumeOfHeads = VolumeOfHead * (NumberOfLipids+DisplacedHeads);
    VolumeOfTails = VolumeOfTail * (NumberOfLipids+DisplacedAlkyl); 
    VolumeOfMethyls = VolumeOfMethyl * (NumberOfLipids+DisplacedMethyl);

    HeightOfMethyl = (HeightOfCore+2*Dummy2)/(VolumeOfTails/VolumeOfMethyls + 1);
    if (HeightOfMethyl < 0.0) {
        printf("Warning: Height of methyl layer is below zero!!! \n");
    }
    SemiMinorAxisOfCore = sqrt(VolumeOfMethyls/(HeightOfMethyl*pi*RatioBetweenAxis));
    SemiMajorAxisOfCore = SemiMinorAxisOfCore * RatioBetweenAxis;
    
    AreaOfDisc=pi*SemiMinorAxisOfCore*SemiMajorAxisOfCore;
    HeightOfLipids = VolumeOfHeads / AreaOfDisc + HeightOfCore;
    ThicknessOfBelt = - (SemiMajorAxisOfCore + SemiMinorAxisOfCore) / 2.0 +
        sqrt(pow(SemiMajorAxisOfCore + SemiMinorAxisOfCore, 2) / 4.0 +
                VolumeOfBelt / (pi * HeightOfBelt));
    // Assign values of the constraints
    Constraints[6] = HeightOfLipids;
    Constraints[7] = HeightOfCore;
    Constraints[8] = HeightOfMethyl;
    Constraints[9] = VerticalAxisOfEllipsoid;
    Constraints[10] = ScaleFactorOfCaps;
    Constraints[11] = 0.0;
    Constraints[12] = 0.0;
    Constraints[13] = 0.0;
    Constraints[14] = 0.0;
    Constraints[15] = 0.0;
    Constraints[16] = 0.0;
    Constraints[17] = 0.0;
    Constraints[18] = 0.0;
    Constraints[19] = 0.0;
    Constraints[20] = 0.0;
    Constraints[21] = Concentration;
    Constraints[22] = 0.0;
    Constraints[23] = SemiMajorAxisOfCore;
    Constraints[24] = SemiMinorAxisOfCore;
    Constraints[25] = 0.0;
    Constraints[26] = 0.0;
    Constraints[27] = ThicknessOfBelt;
    Constraints[28] = 0.0;
    Constraints[29] = 0;
    Constraints[30] = 0;
    //printf("VolTags %.4lf, DensTags %.4lf\n",Constraints[29],Constraints[30]);
    //getchar();
}
