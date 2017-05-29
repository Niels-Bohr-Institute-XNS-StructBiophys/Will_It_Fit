void ComputeConstraints(double * Parameters, double * VolumesOfMolecules, double * ScatteringLengths, double Contrast,
                        double Concentration, double * Constraints, struct Protein ProteinStructure, struct UserDefined * UserDefinedStructure)
{
    /// Declarations
    // Dummy integer
    int i;

    // Variables describing the lipid core
    double HeightOfCore;
    double HeightOfBelt;
    double HeightOfLipids;
    double HeightOfMethyl;

    // Variables describing scattering lengths
    double ScatteringLengthDensityOfHeads;
    double ScatteringLengthDensityOfTail;
    double ScatteringLengthDensityOfBelt;
    double ScatteringLengthDensityOfWater;
    double ScatteringLengthDensityOfMethyl;

    double ScatteringLengthOfBelt;
    double ScatteringLengthOfCore;
    double ScatteringLengthOfWater;
    double ScatteringLengthOfCaps;
    double ScatteringLengthOfMethyl;

    // Variables describing the valumes
    double VolumeOfHead;
    double VolumeOfTail;
    double VolumeOfBelt;
    double VolumeOfWater;
    double VolumeOfMethyl;

    double CorrectionToVolumeOfHead;
    double CorrectionToVolumeOfTail;
    double CorrectionToVolumeOfBelt;
    double CorrectionToVolumeOfWater;

    // Variables describing water
    double NumberOfWaterAtHead;

    // Variables describing the lipids
    double SigmaNLipids;
    double AreaPerHeadgroupOfLipids;
    double NumberOfPeptidesMean;

    // Variables used in numerical integral
    double NLipids;
    double NLipidsMean;
    double NLipidsMin;
    double NLipidsMax;
    double NLipidsStep;
    int NumberOfStepsInNLipids = 50;
    double Prefactor;
    double Weight;
    double Contribution;
    double Sigma;

    /// Get parameters from Parameters
    // Scale factors
    CorrectionToVolumeOfBelt  = Parameters[9];
    CorrectionToVolumeOfTail  = Parameters[10];
    CorrectionToVolumeOfHead  = Parameters[10];
    CorrectionToVolumeOfWater = Parameters[19];

    // Other parameters
    NumberOfWaterAtHead      = fabs(Parameters[5]);
    AreaPerHeadgroupOfLipids = Parameters[1];
    HeightOfBelt             = Parameters[2];
    NLipidsMean              = Parameters[3];
    SigmaNLipids             = Parameters[4];

    /// Get parameters from variable VolumesOfMolecules
    VolumeOfWater  = VolumesOfMolecules[0];
    VolumeOfHead   = VolumesOfMolecules[1] * CorrectionToVolumeOfHead + NumberOfWaterAtHead * VolumesOfMolecules[0] * CorrectionToVolumeOfWater;
    VolumeOfTail   = VolumesOfMolecules[2] * CorrectionToVolumeOfTail;
    VolumeOfMethyl = VolumesOfMolecules[3] * CorrectionToVolumeOfTail;
    VolumeOfBelt   = VolumesOfMolecules[4] * CorrectionToVolumeOfBelt;

    /// Compute the height of the different parts of the lipids
    HeightOfLipids = 2.0 * (VolumeOfHead + VolumeOfTail + VolumeOfMethyl) / AreaPerHeadgroupOfLipids;
    HeightOfCore   = 2.0 * (VolumeOfTail + VolumeOfMethyl) / AreaPerHeadgroupOfLipids;
    HeightOfMethyl = 2.0 * VolumeOfMethyl / AreaPerHeadgroupOfLipids;

    /// Obtain scattering length densities
    ScatteringLengthDensityOfWater  = ScatteringLengths[0];
    ScatteringLengthDensityOfHeads  = ScatteringLengths[1] + NumberOfWaterAtHead * ScatteringLengths[0];
    ScatteringLengthDensityOfTail   = ScatteringLengths[2];
    ScatteringLengthDensityOfMethyl = ScatteringLengths[3];
    ScatteringLengthDensityOfBelt   = ScatteringLengths[4];

    /// Derive scattering lengths
    ScatteringLengthOfWater  = ScatteringLengthDensityOfWater / VolumeOfWater;
    ScatteringLengthOfCaps   = ScatteringLengthDensityOfHeads / VolumeOfHead - ScatteringLengthOfWater;
    ScatteringLengthOfCore   = ScatteringLengthDensityOfTail / VolumeOfTail - ScatteringLengthOfWater;
    ScatteringLengthOfMethyl = ScatteringLengthDensityOfMethyl / VolumeOfMethyl - ScatteringLengthOfWater;
    ScatteringLengthOfBelt   = ScatteringLengthDensityOfBelt / VolumeOfBelt - ScatteringLengthOfWater;

    /// Compute mean number of peptides in a disc
    NLipidsMin = NLipidsMean - 3.0 * NLipidsMean * SigmaNLipids;

    if (NLipidsMin < 0.0) {
        NLipidsMin = 0.0;
    }

    NLipidsMax = NLipidsMean + 3.0 * NLipidsMean * SigmaNLipids;
    NLipidsStep = (NLipidsMax - NLipidsMin) / (1.0 * NumberOfStepsInNLipids);

    Sigma =  NLipidsMean * SigmaNLipids;

    Prefactor = 1.0 / (Sigma * sqrt(2.0 * pi));

    NumberOfPeptidesMean = 0.0;

    double CorrectionDueToUnusedDiscs = 0.0;

    for (i = 0; i < NumberOfStepsInNLipids; ++i) {
        NLipids = (i + 0.5) * NLipidsStep + NLipidsMin;

        Weight = exp(- pow(NLipids - NLipidsMean, 2) / (2.0 * pow(Sigma, 2)));

        Contribution = (pi * HeightOfBelt * pow(Parameters[0], 2) +
                        pi * HeightOfBelt * Parameters[0] * sqrt(2.0 * NLipids * AreaPerHeadgroupOfLipids / pi)) / VolumeOfBelt;

        NumberOfPeptidesMean += Prefactor * Weight * NLipidsStep * Contribution;

        CorrectionDueToUnusedDiscs += Prefactor * Weight * NLipidsStep;
    }

    /// Compute values used to rescale intensity
    double LipidsPerPeptideInDiscs = NLipidsMean / NumberOfPeptidesMean;
    double LipidsPerPeptideInSolution = Parameters[21];
    double RatioOfPeptidesInDisc = LipidsPerPeptideInSolution / LipidsPerPeptideInDiscs;
    double VolumeOfPeptideInTrimer = VolumeOfBelt * Parameters[23];
    double HeightOfPeptideInTrimer = VolumeOfPeptideInTrimer / (pi * pow(Parameters[22], 2));

    /// Assign values of the constraints
    Constraints[0] = ScatteringLengthOfCaps;
    Constraints[1] = ScatteringLengthOfCore;
    Constraints[2] = ScatteringLengthOfMethyl;
    Constraints[3] = ScatteringLengthOfBelt;
    Constraints[4] = 0.0;
    Constraints[5] = 0.0;
    Constraints[6] = HeightOfLipids;
    Constraints[7] = HeightOfCore;
    Constraints[8] = HeightOfMethyl;
    Constraints[9] = 0.0;
    Constraints[10] = 0.0;
    Constraints[11] = 0.0;
    Constraints[12] = 0.0;
    Constraints[13] = VolumesOfMolecules[2] * CorrectionToVolumeOfHead;
    Constraints[14] = VolumeOfTail;
    Constraints[15] = VolumeOfMethyl;
    Constraints[16] = VolumeOfBelt;
    Constraints[17] = NumberOfPeptidesMean;
    Constraints[18] = RatioOfPeptidesInDisc;

    /// Geometric properties of the disc
    Constraints[20] = Parameters[3];
    Constraints[21] = RatioOfPeptidesInDisc * Concentration / NumberOfPeptidesMean / CorrectionDueToUnusedDiscs;
    Constraints[22] = (1.0 - RatioOfPeptidesInDisc) * Concentration / 3.0;
    Constraints[23] = HeightOfPeptideInTrimer;
    Constraints[24] = 0.0;
    Constraints[25] = Parameters[3] * VolumeOfTail;
    Constraints[26] = Parameters[3] * VolumeOfHead;
    Constraints[27] = 0.0;
    Constraints[28] = Parameters[3] * VolumeOfMethyl;
}
