void ComputeConstraints(double * Parameters, double * VolumesOfMolecules, double * ScatteringLengths, double Contrast,
                        double Concentration, double * Constraints, struct Protein ProteinStructure, struct UserDefined * UserDefinedStructure)
{
    /// Declarations
    // Variables describing the lipid core
    double HeightOfCore;
    double HeightOfBelt;
    double HeightOfLipids;
    double HeightOfMethyl;
    double HeightOfMembraneProtein;

    // Variables describing scattering lengths
    double ScatteringLengthDensityOfHeads;
    double ScatteringLengthDensityOfTail;
    double ScatteringLengthDensityOfBelt;
    double ScatteringLengthDensityOfWater;
    double ScatteringLengthDensityOfMethyl;
    double ScatteringLengthDensityOfMembraneProtein;

    double ScatteringLengthOfBelt;
    double ScatteringLengthOfCore;
    double ScatteringLengthOfWater;
    double ScatteringLengthOfCaps;
    double ScatteringLengthOfMethyl;
    double ScatteringLengthOfMembraneProtein;

    // Variables describing the valumes
    double VolumeOfHead;
    double VolumeOfTail;
    double VolumeOfBelt;
    double VolumeOfWater;
    double VolumeOfMethyl;
    double VolumeOfMembraneProtein;
    double AngularOffset;
    double DisplacementOfMP;

    double CorrectionToVolumeOfHead;
    double CorrectionToVolumeOfTail;
    double CorrectionToVolumeOfBelt;
    double CorrectionToVolumeOfWater;
    double CorrectionToVolumeOfMembraneProtein;

    // Variables describing water
    double NumberOfWaterAtBelt;
    double NumberOfWaterAtHead;

    // Variables describing the properties of the disc
    double RatioBetweenAxis;
    double ThicknessOfBelt;
    double RatioBetweenAxisOfMembraneProtein;
    double MinorAxisOfMembraneProtein;
    double MajorAxisOfMembraneProtein;

    double MajorAxisOfCore;
    double MinorAxisOfCore;

    // Variables describing the lipids
    double AverageAreaOfCore;
    double NumberOfLipids;
    double AreaPerHeadgroupOfLipids;
    double AreaOfMembraneProtein;


    /// Get parameters from Parameters
    // Scale factors
    CorrectionToVolumeOfBelt = Parameters[9];
    CorrectionToVolumeOfTail = Parameters[10];
    CorrectionToVolumeOfHead = Parameters[10];
    CorrectionToVolumeOfMembraneProtein = Parameters[11];
    CorrectionToVolumeOfWater = Parameters[19];

    // Other parameters
    NumberOfWaterAtHead               = fabs(Parameters[5]);
    NumberOfWaterAtBelt               = fabs(Parameters[6]);
    AreaPerHeadgroupOfLipids          = Parameters[1];
    HeightOfBelt                      = Parameters[2];
    NumberOfLipids                    = Parameters[3];
    RatioBetweenAxis                  = fabs(Parameters[0]);
    RatioBetweenAxisOfMembraneProtein = fabs(Parameters[21]);

    /// Get parameters from variable VolumesOfMolecules
    VolumeOfWater           = VolumesOfMolecules[0];
    VolumeOfHead            = VolumesOfMolecules[1] * CorrectionToVolumeOfHead + NumberOfWaterAtHead * VolumesOfMolecules[0] * CorrectionToVolumeOfWater;
    VolumeOfTail            = VolumesOfMolecules[2] * CorrectionToVolumeOfTail;
    VolumeOfMethyl          = VolumesOfMolecules[3] * CorrectionToVolumeOfTail;
    VolumeOfBelt            = 2.0 * VolumesOfMolecules[4] * CorrectionToVolumeOfBelt + 2.0 * NumberOfWaterAtBelt * VolumesOfMolecules[0];
    VolumeOfMembraneProtein = VolumesOfMolecules[5] * CorrectionToVolumeOfMembraneProtein;

    /// Compute geometry of disc
    HeightOfMembraneProtein = Parameters[22];
    AreaOfMembraneProtein   = VolumeOfMembraneProtein / HeightOfMembraneProtein;

    MinorAxisOfMembraneProtein = sqrt(fabs(AreaOfMembraneProtein / (pi * RatioBetweenAxisOfMembraneProtein)));
    MajorAxisOfMembraneProtein = RatioBetweenAxisOfMembraneProtein * MinorAxisOfMembraneProtein;

	AverageAreaOfCore = NumberOfLipids * AreaPerHeadgroupOfLipids / 2.0 + AreaOfMembraneProtein;

    MinorAxisOfCore = sqrt(fabs(AverageAreaOfCore / (pi * RatioBetweenAxis)));
    MajorAxisOfCore = MinorAxisOfCore * RatioBetweenAxis;

    HeightOfLipids = 2.0 * (VolumeOfHead + VolumeOfTail + VolumeOfMethyl) / AreaPerHeadgroupOfLipids;
    HeightOfCore   = 2.0 * (VolumeOfTail + VolumeOfMethyl) / AreaPerHeadgroupOfLipids;
    HeightOfMethyl = 2.0 * VolumeOfMethyl / AreaPerHeadgroupOfLipids;

    /// Obtain scattering length densities
    ScatteringLengthDensityOfWater           = ScatteringLengths[0];
    ScatteringLengthDensityOfHeads           = ScatteringLengths[1] + NumberOfWaterAtHead * ScatteringLengths[0];
    ScatteringLengthDensityOfTail            = ScatteringLengths[2];
    ScatteringLengthDensityOfMethyl          = ScatteringLengths[3];
    ScatteringLengthDensityOfBelt            = 2.0 * (ScatteringLengths[4] + NumberOfWaterAtBelt * ScatteringLengths[0]);
    ScatteringLengthDensityOfMembraneProtein = ScatteringLengths[6];

    /// Derive excess scattering length densities
    ScatteringLengthOfWater           = ScatteringLengthDensityOfWater / VolumeOfWater;
    ScatteringLengthOfCaps            = ScatteringLengthDensityOfHeads / VolumeOfHead - ScatteringLengthOfWater;
    ScatteringLengthOfCore            = ScatteringLengthDensityOfTail / VolumeOfTail - ScatteringLengthOfWater;
    ScatteringLengthOfMethyl          = ScatteringLengthDensityOfMethyl / VolumeOfMethyl - ScatteringLengthOfWater;
    ScatteringLengthOfBelt            = ScatteringLengthDensityOfBelt / VolumeOfBelt - ScatteringLengthOfWater;
    ScatteringLengthOfMembraneProtein = ScatteringLengthDensityOfMembraneProtein / VolumeOfMembraneProtein - ScatteringLengthOfWater;

    /// Angle between cylinders and displacement
    AngularOffset = Parameters[24];
    DisplacementOfMP = Parameters[25];

    /// Assign values of the constraints
    Constraints[0]  = ScatteringLengthOfCaps * 1e24;
    Constraints[1]  = ScatteringLengthOfCore * 1e24;
    Constraints[2]  = ScatteringLengthOfMethyl * 1e24;
    Constraints[3]  = ScatteringLengthOfBelt * 1e24;
    Constraints[4]  = 0.0;
    Constraints[5]  = ScatteringLengthOfMembraneProtein * 1e24;
    Constraints[6]  = HeightOfLipids;
    Constraints[7]  = HeightOfCore;
    Constraints[8]  = HeightOfMethyl;
    Constraints[9]  = HeightOfMembraneProtein;
    Constraints[10] = MinorAxisOfMembraneProtein;
    Constraints[11] = MajorAxisOfMembraneProtein;
    Constraints[12] = Concentration;
    Constraints[13] = VolumesOfMolecules[2] * CorrectionToVolumeOfHead;
    Constraints[14] = VolumeOfTail;
    Constraints[15] = VolumeOfMethyl;
    Constraints[16] = VolumeOfBelt;
    Constraints[17] = 0.0f;
    Constraints[18] = VolumeOfMembraneProtein;
    Constraints[19] = AngularOffset;
    Constraints[20] = DisplacementOfMP;

    /// Geometric properties of the disc
    ThicknessOfBelt = -(MajorAxisOfCore + MinorAxisOfCore) / 2.0 +
                         sqrt(pow(MajorAxisOfCore + MinorAxisOfCore, 2) / 4.0f +
                              VolumeOfBelt / (pi * HeightOfBelt));

    Constraints[21] = Concentration;
    Constraints[22] = Parameters[3];
    Constraints[23] = MajorAxisOfCore;
    Constraints[24] = MinorAxisOfCore;
    Constraints[25] = Parameters[3] * VolumeOfTail;
    Constraints[26] = Parameters[3] * VolumeOfHead;
    Constraints[27] = ThicknessOfBelt;
    Constraints[28] = Parameters[3] * VolumeOfMethyl;
}
