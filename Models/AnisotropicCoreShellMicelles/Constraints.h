void ComputeConstraints(double * Parameters, double * VolumesOfMolecules, double * ScatteringLengths, double Contrast,
                        double Concentration, double * Constraints, struct Protein ProteinStructure, struct UserDefined * UserDefinedStructure)
{
    /// Declarations
    // Variables describing scattering lengths
    double ExcessScatteringLengthDensityOfHeads;
    double ExcessScatteringLengthDensityOfTails;
    double ScatteringLengthOfHeads;
    double ScatteringLengthOfTails;

    // Variables describing the aggregates
    double VolumeOfHead;
    double VolumeOfTail;
    double VolumeOfCore;
    double VolumeOfShell;
    double CorrectionToVolumeOfHead;
    double CorrectionToVolumeOfTail;
    double MinorRadiusOfCore;
    double MajorRadiusOfCore;
    double AxisRatioOfCore;
    double NumberOfMoleculesPerAggregate;
    double ThicknessOfShell;

    // Variables describing water
    double VolumeOfWater;
    double CorrectionToVolumeOfWater;
    double NumberOfWaterAtHead;
    double ScatteringLengthOfWater;
    double ScatteringLengthDensityOfWater;

    // Dummy variablas
    double a;
    double b;
    double c;
    double d;

    double p;
    double q;
    double r;

    /// Get parameters from Parameters
    NumberOfWaterAtHead       = Parameters[3];
    CorrectionToVolumeOfHead  = Parameters[4];
    CorrectionToVolumeOfTail  = Parameters[5];
    CorrectionToVolumeOfWater = Parameters[6];
    MinorRadiusOfCore         = Parameters[7];
    AxisRatioOfCore           = Parameters[11];

    /// Get parameters from variable VolumesOfMolecules
    VolumeOfWater = VolumesOfMolecules[0] * CorrectionToVolumeOfWater;
    VolumeOfHead  = VolumesOfMolecules[1] * CorrectionToVolumeOfHead + NumberOfWaterAtHead * VolumeOfWater;
    VolumeOfTail  = VolumesOfMolecules[2] * CorrectionToVolumeOfTail;

    /// Obtain scattering length
    ScatteringLengthOfWater = ScatteringLengths[0];
    ScatteringLengthOfHeads = ScatteringLengths[1] + NumberOfWaterAtHead * ScatteringLengthOfWater;
    ScatteringLengthOfTails = ScatteringLengths[2];

    /// Derive scattering length densities
    ScatteringLengthDensityOfWater       = ScatteringLengthOfWater / VolumeOfWater;
    ExcessScatteringLengthDensityOfHeads = ScatteringLengthOfHeads / VolumeOfHead - ScatteringLengthDensityOfWater;
    ExcessScatteringLengthDensityOfTails = ScatteringLengthOfTails / VolumeOfTail - ScatteringLengthDensityOfWater;

    /// Geometry of the aggregate
    MajorRadiusOfCore             = MinorRadiusOfCore * AxisRatioOfCore;
    VolumeOfCore                  = 4.0 / 3.0 * pi * AxisRatioOfCore * pow(MinorRadiusOfCore, 3);
    NumberOfMoleculesPerAggregate = VolumeOfCore / VolumeOfTail;
    VolumeOfShell                 = NumberOfMoleculesPerAggregate * VolumeOfHead;

    // The thickness of the shell is determined by solving the following third degree equation
    // V = 4/3 pi (eR + D) (R + D) (R + D) - 4/3 pi e R^3
    // using Cardanos formula

    a = 1.0;
    b = 2.0 * MinorRadiusOfCore + AxisRatioOfCore * MinorRadiusOfCore;
    c = 2.0 * AxisRatioOfCore * pow(MinorRadiusOfCore, 2) + pow(MinorRadiusOfCore, 2);
    d = - 3.0 * VolumeOfShell / (4.0 * pi);

    p = - b / (3.0 * a);
    q = pow(p, 3) + (b * c - 3.0 * a * d) / (6.0 * pow(a, 2));
    r = c / (3.0 * a);

    ThicknessOfShell = cbrt(q + sqrt(pow(q, 2) + pow(r - pow(p, 2), 3))) +
                       cbrt(q - sqrt(pow(q, 2) + pow(r - pow(p, 2), 3))) +
                       p;

    /// Assign values of the constraints
    Constraints[0] = ExcessScatteringLengthDensityOfHeads;
    Constraints[1] = ExcessScatteringLengthDensityOfTails;
    Constraints[2] = Concentration;
    Constraints[3] = VolumeOfHead;
    Constraints[4] = VolumeOfTail;
    Constraints[5] = NumberOfMoleculesPerAggregate;
    Constraints[6] = MinorRadiusOfCore;
    Constraints[7] = MajorRadiusOfCore;
    Constraints[8] = ThicknessOfShell;
}
