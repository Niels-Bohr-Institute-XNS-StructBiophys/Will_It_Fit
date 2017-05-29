void ComputeConstraints(double * Parameters, double * VolumesOfMolecules, double * ScatteringLengths, double Contrast,
                        double Concentration, double * Constraints, struct Protein ProteinStructure, struct UserDefined * UserDefinedStructure)
{
    double MinorRadius = Parameters[3];
    double MajorRadius = Parameters[4];
    double Height = Parameters[5];
    double ConcentrationValue = Concentration;

    double ScatteringLength = ScatteringLengths[1];
    double ScatteringLengthOfSolvent = ScatteringLengths[0];

    double Volume = VolumesOfMolecules[1];
    double VolumeOfSolvent = VolumesOfMolecules[0];

    double ScatteringLengthDensity = ScatteringLength / Volume;
    double ScatteringLengthDensityOfSolvent = ScatteringLengthOfSolvent / VolumeOfSolvent;

    double ExcessScatteringLengthDensity = ScatteringLengthDensity - ScatteringLengthDensityOfSolvent;

    Constraints[0] = ExcessScatteringLengthDensity;
    Constraints[1] = MinorRadius;
    Constraints[2] = MajorRadius;
    Constraints[3] = Height;
    Constraints[4] = ConcentrationValue;
}
