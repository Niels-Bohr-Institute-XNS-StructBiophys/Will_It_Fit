double PeptideDiscs(double q, double HeightOfBelt, double HeightOfLipids, double HeightOfCore, double HeightOfMethyl, double VolumeOfBelt,
                double ScatteringLengthOfCaps, double ScatteringLengthOfCore, double ScatteringLengthOfBelt, double ScatteringLengthOfMethyl,
                double NumberOfLipids, double AggregationNumber, double Radlipidcore,
                double VolumeOfCore, double VolumeOfCaps, double VolumeOfMethyl, double ThicknessOfBelt)
{
    /// Declarations
    // Declarations of dummy variables used in computations
    int i;

    // Declarations used to compute return value
    double SumOverAlpha;
    double PartialSumOfNanodisc;
    double ReturnValue;

    // Declarations used to describe the formfactors
    double PsiCore;
    double PsiLipids;
    double PsiOuterBelt;
    double PsiInnerBelt;
    double PsiMethyl;

    // Declare volume weights
    double ComputedVolumeOfOuterBelt;
    double ComputedVolumeOfInnerBelt;
    double ComputedVolumeOfLipids;
    double ComputedVolumeOfCore;
    double ComputedVolumeOfMethyl;

    // Declare radii
    double InnerRadiusOfBelt;
    double OuterRadiusOfBelt;

    // Declarations used to compute amplitudes
    double FormfactorOfBelt;
    double FormfactorOfCaps;
    double FormfactorOfCore;
    double FormfactorOfMethyl;

    // Declarations of amplitudes
    double AmplitudeOfBelt;
    double AmplitudeOfCaps;
    double AmplitudeOfCore;
    double AmplitudeOfMethyl;
    double AmplitudeOfDisc;

    // Introduce the parameters used in the integrations
    const int NumberOfStepsInAlpha = 50;
    const double AlphaStepSize = pi / (2.0 * NumberOfStepsInAlpha);
    double Alpha;

    /// Begin integrations
    // Integrate over Alpha
    SumOverAlpha = 0.0;

    for (i = 1; i <= NumberOfStepsInAlpha; ++i) {
        Alpha = AlphaStepSize * i;

        /// Compute characteristics of belt
        // Compute the inner radius of the belt
        InnerRadiusOfBelt = Radlipidcore;

        // Compute the outer radius of the belt
        OuterRadiusOfBelt = Radlipidcore + ThicknessOfBelt;

        // Compute characteristics of the belt
        PsiOuterBelt = FormfactorCylinder(q, OuterRadiusOfBelt, HeightOfBelt, Alpha);
        PsiInnerBelt = FormfactorCylinder(q, InnerRadiusOfBelt, HeightOfBelt, Alpha);

        // Compute volumes
        ComputedVolumeOfOuterBelt = pi * HeightOfBelt * pow(OuterRadiusOfBelt, 2);
        ComputedVolumeOfInnerBelt = pi * HeightOfBelt * pow(InnerRadiusOfBelt, 2);

        // Compute form factor of belt
        FormfactorOfBelt = (ComputedVolumeOfOuterBelt * PsiOuterBelt - ComputedVolumeOfInnerBelt * PsiInnerBelt) /
                           (ComputedVolumeOfOuterBelt - ComputedVolumeOfInnerBelt);

        /// Compute characteristics of core
        // Compute characteristics for the methylcore
        PsiLipids = FormfactorCylinder(q, InnerRadiusOfBelt, HeightOfLipids, Alpha);
        PsiCore   = FormfactorCylinder(q, InnerRadiusOfBelt, HeightOfCore, Alpha);
        PsiMethyl = FormfactorCylinder(q, InnerRadiusOfBelt, HeightOfMethyl, Alpha);

        // Compute volumes
        ComputedVolumeOfLipids = pi * pow(InnerRadiusOfBelt, 2) * HeightOfLipids;
        ComputedVolumeOfCore   = pi * pow(InnerRadiusOfBelt, 2) * HeightOfCore;
        ComputedVolumeOfMethyl = pi * pow(InnerRadiusOfBelt, 2) * HeightOfMethyl;

        // Introduce functions to calculate the normalized amplitudes for the different parts of the disc
        FormfactorOfCaps   = (ComputedVolumeOfLipids * PsiLipids - ComputedVolumeOfCore * PsiCore) / (ComputedVolumeOfLipids - ComputedVolumeOfCore);
        FormfactorOfCore   = (ComputedVolumeOfCore * PsiCore - ComputedVolumeOfMethyl * PsiMethyl) / (ComputedVolumeOfCore - ComputedVolumeOfMethyl);
        FormfactorOfMethyl = PsiMethyl;

        /// Compute form factors of the different parts
        // Compute the different parts of the form factor
        AmplitudeOfBelt   = ScatteringLengthOfBelt * VolumeOfBelt * FormfactorOfBelt;
        AmplitudeOfCaps   = ScatteringLengthOfCaps * VolumeOfCaps * FormfactorOfCaps;
        AmplitudeOfCore   = ScatteringLengthOfCore * VolumeOfCore * FormfactorOfCore;
        AmplitudeOfMethyl = ScatteringLengthOfMethyl * VolumeOfMethyl * FormfactorOfMethyl;

        AmplitudeOfDisc = AmplitudeOfBelt + AmplitudeOfCore + AmplitudeOfCaps + AmplitudeOfMethyl;
        PartialSumOfNanodisc = pow(AmplitudeOfDisc, 2);

        SumOverAlpha += PartialSumOfNanodisc * sin(Alpha);
    }

    /// Rescale results
    ReturnValue = SumOverAlpha * AggregationNumber * AlphaStepSize;

    return ReturnValue;
}
