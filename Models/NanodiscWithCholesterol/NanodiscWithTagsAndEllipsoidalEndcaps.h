double NanodiscWithTagsAndEllipsoidalEndcaps(double q, double HeightOfBelt, double HeightOfLipids, double HeightOfCore, double HeightOfMethyl,
                                             double VolumeOfBelt, double VolumeOfMethyl, double ScatteringLengthOfCaps, double ScatteringLengthOfCore,
                                             double ScatteringLengthOfMethyl, double ScatteringLengthOfBelt, double NumberOfLipids,
                                             double AggregationNumber, double MajorAxisOfCore, double MinorAxisOfCore, double VolumeOfCore,
                                             double VolumeOfCaps, double ThicknessOfBelt, double VolumeOfTag, double RadiusOfGyrationForTag,
                                             double ScatteringLengthOfTag, double HeightOfCurve, double ScaleFactorOfEndcaps,
                                             double VerticalAxisOfEndcaps, double VerticalShiftOfEndcaps)
{
    /// Declarations
    // Declarations of dummy variables used in computations
    int i;
    int j;

    // Declarations used to compute return value
    double SumOverAlpha;
    double PartialSumOfNanodisc;
    double PartialSumWithTags;
    double RescaledSum;
    double ReturnValue;

    // Declarations of partial sums
    double SumOverBeta;
    double DiscChainCorrelation;
    double AutoCorrelation;
    double CrossCorrelation;

    // Declarations used to describe the formfactors
    double PsiCore;
    double PsiLipids;
    double PsiOuterBelt;
    double PsiInnerBelt;
    double PsiTags;
    double PsiMethyl;

    // Declare volume weights
    double ComputedVolumeOfOuterBelt;
    double ComputedVolumeOfInnerBelt;
    double ComputedVolumeOfLipids;
    double ComputedVolumeOfCore;
    double ComputedVolumeOfEndcaps;
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
    double AmplitudeOfDisc;
    double AmplitudeOfMethyl;

    // Introduce the parameters used in the integrations
    const int NumberOfStepsInAlpha = 50;
    const int NumberOfStepsInBeta = 50;

    const double AlphaStepSize = pi / (2.0 * NumberOfStepsInAlpha);
    const double BetaStepSize = pi / (2.0 * NumberOfStepsInBeta);

    double Alpha;
    double Beta;

    /// Begin integrations
    // Integrate over Alpha
    SumOverAlpha = 0.0;

    for (i = 1; i <= NumberOfStepsInAlpha; ++i) {
        Alpha = AlphaStepSize * i;
        SumOverBeta = 0.0;

        // Integrate over Beta
        for (j = 1; j <= NumberOfStepsInBeta; ++j) {
            Beta = BetaStepSize * j;

            /// Compute characteristics of belt
            // Compute the inner radius of the belt
            InnerRadiusOfBelt = sqrt(pow(MinorAxisOfCore * cos(Beta), 2) + pow(MajorAxisOfCore * sin(Beta), 2));

            // Compute the outer radius of the belt
            OuterRadiusOfBelt = sqrt(pow((MinorAxisOfCore + ThicknessOfBelt) * cos(Beta), 2) +
                                     pow((MajorAxisOfCore + ThicknessOfBelt) * sin(Beta), 2));

            // Compute characteristics of the belt
            PsiOuterBelt = FormfactorCylinder(q, OuterRadiusOfBelt, HeightOfBelt, Alpha);
            PsiInnerBelt = FormfactorCylinder(q, InnerRadiusOfBelt, HeightOfBelt, Alpha);

            // Compute volumes
            ComputedVolumeOfOuterBelt = pi * HeightOfBelt * (MinorAxisOfCore + ThicknessOfBelt) * (MajorAxisOfCore + ThicknessOfBelt);
            ComputedVolumeOfInnerBelt = pi * HeightOfBelt * MinorAxisOfCore * MajorAxisOfCore;

            // Compute form factor of belt
            FormfactorOfBelt = (ComputedVolumeOfOuterBelt * PsiOuterBelt - ComputedVolumeOfInnerBelt * PsiInnerBelt) /
                               (ComputedVolumeOfOuterBelt - ComputedVolumeOfInnerBelt);

            /// Compute characteristics of core
            // Formfactors
            PsiLipids = FormfactorCylinderWithEndcaps(q, Alpha, Beta, MinorAxisOfCore, MajorAxisOfCore, HeightOfLipids,
                                                      ScaleFactorOfEndcaps, VerticalAxisOfEndcaps);

            PsiCore = FormfactorCylinderWithEndcaps(q, Alpha, Beta, MinorAxisOfCore, MajorAxisOfCore, HeightOfCore,
                                                    ScaleFactorOfEndcaps, VerticalAxisOfEndcaps);

            PsiMethyl = FormfactorCylinder(q, InnerRadiusOfBelt, HeightOfMethyl, Alpha);

            // Compute volumes
            ComputedVolumeOfEndcaps = 2.0 * pi * MinorAxisOfCore * MajorAxisOfCore * pow(ScaleFactorOfEndcaps, 2) *
                                      (2.0 * VerticalAxisOfEndcaps / 3.0 - fabs(VerticalShiftOfEndcaps) +
                                       pow(fabs(VerticalShiftOfEndcaps), 3) / (3.0 * pow(VerticalAxisOfEndcaps, 2)));

            ComputedVolumeOfLipids  = pi * MinorAxisOfCore * MajorAxisOfCore * HeightOfLipids + ComputedVolumeOfEndcaps;
            ComputedVolumeOfCore    = pi * MinorAxisOfCore * MajorAxisOfCore * HeightOfCore   + ComputedVolumeOfEndcaps;
            ComputedVolumeOfMethyl  = pi * MinorAxisOfCore * MajorAxisOfCore * HeightOfMethyl;

            // Introduce functions to calculate the normalized amplitudes for the different parts of the disc
            FormfactorOfCaps   = (ComputedVolumeOfLipids * PsiLipids - ComputedVolumeOfCore * PsiCore) /
                                 (ComputedVolumeOfLipids - ComputedVolumeOfCore);

            FormfactorOfCore   = (ComputedVolumeOfCore * PsiCore - ComputedVolumeOfMethyl * PsiMethyl) /
                                 (ComputedVolumeOfCore - ComputedVolumeOfMethyl);

            FormfactorOfMethyl = PsiMethyl;

            /// Compute form factors of the different parts
            AmplitudeOfBelt   = ScatteringLengthOfBelt   * VolumeOfBelt   * FormfactorOfBelt;
            AmplitudeOfCaps   = ScatteringLengthOfCaps   * VolumeOfCaps   * FormfactorOfCaps;
            AmplitudeOfCore   = ScatteringLengthOfCore   * VolumeOfCore   * FormfactorOfCore;
            AmplitudeOfMethyl = ScatteringLengthOfMethyl * VolumeOfMethyl * FormfactorOfMethyl;

            AmplitudeOfDisc = AmplitudeOfBelt + AmplitudeOfCore + AmplitudeOfCaps + AmplitudeOfMethyl;
            PartialSumOfNanodisc = pow(AmplitudeOfDisc, 2);

            /// Include the his-tags at the outer rim of the nanodisc belt
            // Replaced by expression for thin rod-like shell
            PsiTags = FormfactorRodWithoutEndcaps(q, OuterRadiusOfBelt + RadiusOfGyrationForTag, HeightOfBelt, Alpha);

            // Compute chain-chain crosscorrelation function
            CrossCorrelation = 2.0 * pow(ScatteringLengthOfTag * VolumeOfTag, 2) * pow(PsiTags, 2) *
                               pow(FormfactorHammouda(q * RadiusOfGyrationForTag), 2);

            // Compute chain-chain autocorrelation function
            AutoCorrelation = 2.0 * pow(ScatteringLengthOfTag * VolumeOfTag, 2) * FormfactorDebye(q, RadiusOfGyrationForTag);

            // Compute chain-disc crosscorrelation function
            DiscChainCorrelation = 4.0 * AmplitudeOfDisc * (ScatteringLengthOfTag * VolumeOfTag) *
                                   PsiTags * FormfactorHammouda(q * RadiusOfGyrationForTag);

            // Summing over the contributions
            PartialSumWithTags = PartialSumOfNanodisc + CrossCorrelation + AutoCorrelation + DiscChainCorrelation;

            // End of the inclusion of his-tags
            SumOverBeta += PartialSumWithTags;

            /// To exclude his-tags, use this line following line instead
            //SumOverBeta += PartialSumOfNanodisc;
        }

        SumOverAlpha += SumOverBeta * sin(Alpha);
    }

    /// Rescale results
    RescaledSum = SumOverAlpha * AggregationNumber;
    ReturnValue = RescaledSum * (2.0 / pi) * BetaStepSize * AlphaStepSize;

    return ReturnValue;
}
