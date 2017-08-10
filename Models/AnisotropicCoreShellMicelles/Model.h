double Model(double q, double * Parameters, double * Constraints, double Contrast, struct Protein ProteinStructure, struct UserDefined * UserDefinedStructure)
{
    /// Declarations
    double ScatteringLengthDensityOfTails;
    double ScatteringLengthDensityOfHeads;
    double MinorRadiusOfCore;
    double MajorRadiusOfCore;
    double ThicknessOfShell;
    double ConcentrationValue;
    double NumberOfMoleculesPerAggregate;
    double Intensity;
    double ShellRoughness = 0.0;
    double CoreRoughness  = 0.0;
    double AxisRatioOfCore;
    double AxisRatioOfShell;
    double Scaling;
    double Background;

    /// Get values from constraints
    ScatteringLengthDensityOfHeads = Constraints[0];
    ScatteringLengthDensityOfTails = Constraints[1];
    ConcentrationValue             = Constraints[2];
    NumberOfMoleculesPerAggregate  = Constraints[5];

    MinorRadiusOfCore = Constraints[6];
    MajorRadiusOfCore = Constraints[7];
    ThicknessOfShell  = Constraints[8];

    AxisRatioOfCore  = MajorRadiusOfCore / MinorRadiusOfCore;
    AxisRatioOfShell = (MajorRadiusOfCore + ThicknessOfShell) / (MinorRadiusOfCore + ThicknessOfShell);

    if (Contrast >= 0.0 && Contrast <= 100.0) {
        ShellRoughness = Parameters[14];
        CoreRoughness  = Parameters[15];
    } else {
        ShellRoughness = Parameters[12];
        CoreRoughness  = Parameters[13];
    }

    /// Compute model
    Intensity = AnisotropicMicelles(q, MinorRadiusOfCore, MajorRadiusOfCore, ScatteringLengthDensityOfHeads, ScatteringLengthDensityOfTails,
                                    ConcentrationValue / NumberOfMoleculesPerAggregate, ThicknessOfShell, ShellRoughness,
                                    CoreRoughness, AxisRatioOfCore, AxisRatioOfShell);

    /// Rescale result and return
    if (Contrast >= 0.0 && Contrast <= 100.0) {
        Scaling    = Parameters[16] * (Contrast / 100.0 * Parameters[0] + (100.0 - Contrast) / 100.0 * Parameters[1]);
        Background = Contrast / 100.0 * Parameters[8] + (100.0 - Contrast) / 100.0 * Parameters[9];
    } else {
        Scaling    = Parameters[16] * Parameters[2];
        Background = Parameters[10];
    }

    Intensity = Intensity * Scaling - Background;

    return Intensity;
}
