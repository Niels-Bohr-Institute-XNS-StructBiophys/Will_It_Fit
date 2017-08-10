double NanodiscWithEllipticCylinderAsMembraneProtein(double q, double HeightOfBelt, double HeightOfLipids, double HeightOfCore, double HeightOfMethyl,
                                                     double VolumeOfBelt, double ScatteringLengthOfCaps, double ScatteringLengthOfCore,
                                                     double ScatteringLengthOfBelt, double ScatteringLengthOfMethyl, double NumberOfLipids, double AggregationNumber,
                                                     double MajorAxisOfCore, double MinorAxisOfCore, double VolumeOfCore, double VolumeOfCaps,
                                                     double VolumeOfMethyl, double ThicknessOfBelt, double MinorAxisOfMembraneProtein,
                                                     double MajorAxisOfMembraneProtein, double HeightOfMembraneProtein, double ScatteringLengthOfMembraneProtein,
                                                     double VolumeOfMembraneProtein, double AngularOffset, double DisplacementOfMP)
{
    /// Declarations
    // Declarations of dummy variables used in computations
    double Dummy1;
    double Dummy2;
    int i;
    int j;

    // Declarations used to compute return value
    double SumOverAlpha;
    double PartialSumOfNanodisc;
    double RescaledSum;
    double ReturnValue;

    // Declarations of partial sums
    double SumOverBeta;

    // Declarations used to describe the formfactors
    double PsiCoreOuter;
    double PsiLipidsOuter;
    double PsiMethylOuter;

    double PsiCoreInner;
    double PsiLipidsInner;
    double PsiMethylInner;

    double PsiOuterBelt;
    double PsiInnerBelt;
    complex PsiMembraneProtein;

    // Declare volume weights
    double ComputedVolumeOfOuterBelt;
    double ComputedVolumeOfInnerBelt;

    double ComputedVolumeOfLipidsOuter;
    double ComputedVolumeOfLipidsInner;

    double ComputedVolumeOfCoreOuter;
    double ComputedVolumeOfCoreInner;

    double ComputedVolumeOfMethylOuter;
    double ComputedVolumeOfMethylInner;

    double ComputedVolumeOfLipids;
    double ComputedVolumeOfCore;
    double ComputedVolumeOfMethyl;

    // Declare radii
    double InnerRadiusOfBelt;
    double OuterRadiusOfBelt;
    double RadiusOfMembraneProtein;

    // Declarations used to compute amplitudes
    double FormfactorOfBelt;
    double FormfactorOfCaps;
    double FormfactorOfCore;
    double FormfactorOfMethyl;
    complex FormfactorOfMembraneProtein;

    double FormfactorOfLipidsWithHole;
    double FormfactorOfCoreWithHole;
    double FormfactorOfMethylWithHole;

    // Declarations of amplitudes
    double AmplitudeOfBelt;
    double AmplitudeOfCaps;
    double AmplitudeOfCore;
    double AmplitudeOfMethyl;
    complex AmplitudeOfDisc;
    complex AmplitudeOfMembraneProtein;

    // Introduce the parameters used in the integrations
    const int NumberOfStepsInAlpha = 50;
    const int NumberOfStepsInBeta = 20;

    const double AlphaStepSize = pi / (2.0f * NumberOfStepsInAlpha);
    const double BetaStepSize = pi / (2.0f * NumberOfStepsInBeta);

    double Alpha;
    double Beta;

    /// Begin integrations
    // Integrate over Alpha
    SumOverAlpha = 0.f;

    for (i = 1; i <= NumberOfStepsInAlpha; ++i) {
        Alpha = AlphaStepSize * i;
        SumOverBeta = 0.0f;

        // Integrate over Beta
        for (j = 1; j <= NumberOfStepsInBeta; ++j) {
            Beta = BetaStepSize * j;

            /// Compute characteristics of belt
            // Compute the inner radius of the belt
            Dummy1 = MinorAxisOfCore * cos(Beta);
            Dummy2 = MajorAxisOfCore * sin(Beta);

            InnerRadiusOfBelt = sqrt(pow(Dummy1, 2) + pow(Dummy2, 2));

            // Compute the outer radius of the belt
            Dummy1 = (MinorAxisOfCore + ThicknessOfBelt) * cos(Beta);
            Dummy2 = (MajorAxisOfCore + ThicknessOfBelt) * sin(Beta);

            OuterRadiusOfBelt = sqrt(pow(Dummy1, 2) + pow(Dummy2, 2));

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
            // Compute the radius of the membrane protein
            Dummy1 = MinorAxisOfMembraneProtein * cos(Beta + AngularOffset);
            Dummy2 = MajorAxisOfMembraneProtein * sin(Beta + AngularOffset);

            RadiusOfMembraneProtein = sqrt(pow(Dummy1, 2) + pow(Dummy2, 2));

            // Compute characteristics for the methylcore
            PsiLipidsOuter = FormfactorCylinder(q, InnerRadiusOfBelt, HeightOfLipids, Alpha);
            PsiCoreOuter = FormfactorCylinder(q, InnerRadiusOfBelt, HeightOfCore, Alpha);
            PsiMethylOuter = FormfactorCylinder(q, InnerRadiusOfBelt, HeightOfMethyl, Alpha);

            PsiLipidsInner = FormfactorCylinder(q, RadiusOfMembraneProtein, HeightOfLipids, Alpha);
            PsiCoreInner = FormfactorCylinder(q, RadiusOfMembraneProtein, HeightOfCore, Alpha);
            PsiMethylInner = FormfactorCylinder(q, RadiusOfMembraneProtein, HeightOfMethyl, Alpha);

            PsiMembraneProtein = FormfactorCylinderDisplaced(q, RadiusOfMembraneProtein, HeightOfMembraneProtein, DisplacementOfMP, Alpha);

            // Compute volumes
            ComputedVolumeOfLipidsOuter = pi * MinorAxisOfCore * MajorAxisOfCore * HeightOfLipids;
            ComputedVolumeOfLipidsInner = pi * MinorAxisOfMembraneProtein * MajorAxisOfMembraneProtein * HeightOfLipids;

            ComputedVolumeOfCoreOuter = pi * MinorAxisOfCore * MajorAxisOfCore * HeightOfCore;
            ComputedVolumeOfCoreInner = pi * MinorAxisOfMembraneProtein * MajorAxisOfMembraneProtein * HeightOfCore;

            ComputedVolumeOfMethylOuter = pi * MinorAxisOfCore * MajorAxisOfCore * HeightOfMethyl;
            ComputedVolumeOfMethylInner = pi * MinorAxisOfMembraneProtein * MajorAxisOfMembraneProtein * HeightOfMethyl;

            // Introduce functions to calculate the normalized amplitudes for the different parts of the disc
            FormfactorOfLipidsWithHole = (ComputedVolumeOfLipidsOuter * PsiLipidsOuter - ComputedVolumeOfLipidsInner * PsiLipidsInner) /
                                         (ComputedVolumeOfLipidsOuter - ComputedVolumeOfLipidsInner);

            FormfactorOfCoreWithHole = (ComputedVolumeOfCoreOuter * PsiCoreOuter - ComputedVolumeOfCoreInner * PsiCoreInner) /
                                       (ComputedVolumeOfCoreOuter - ComputedVolumeOfCoreInner);

            FormfactorOfMethylWithHole = (ComputedVolumeOfMethylOuter * PsiMethylOuter - ComputedVolumeOfMethylInner * PsiMethylInner) /
                                         (ComputedVolumeOfMethylOuter - ComputedVolumeOfMethylInner);

            // Computed ALL the volumes and formfactors
            ComputedVolumeOfLipids = ComputedVolumeOfLipidsOuter - ComputedVolumeOfLipidsInner;
            ComputedVolumeOfCore = ComputedVolumeOfCoreOuter - ComputedVolumeOfCoreInner;
            ComputedVolumeOfMethyl = ComputedVolumeOfMethylOuter - ComputedVolumeOfMethylInner;

            FormfactorOfCaps = (ComputedVolumeOfLipids * FormfactorOfLipidsWithHole - ComputedVolumeOfCore * FormfactorOfCoreWithHole) /
                               (ComputedVolumeOfLipids - ComputedVolumeOfCore);

            FormfactorOfCore = (ComputedVolumeOfCore * FormfactorOfCoreWithHole - ComputedVolumeOfMethyl * FormfactorOfMethylWithHole) /
                               (ComputedVolumeOfCore - ComputedVolumeOfMethyl);

            FormfactorOfMethyl = FormfactorOfMethylWithHole;

            FormfactorOfMembraneProtein = PsiMembraneProtein;

            /// Compute form factors of the different parts
            // Compute the different parts of the form factor
            AmplitudeOfBelt = ScatteringLengthOfBelt * VolumeOfBelt * FormfactorOfBelt;
            AmplitudeOfCaps = ScatteringLengthOfCaps * VolumeOfCaps * FormfactorOfCaps;
            AmplitudeOfCore = ScatteringLengthOfCore * VolumeOfCore * FormfactorOfCore;
            AmplitudeOfMethyl = ScatteringLengthOfMethyl * VolumeOfMethyl * FormfactorOfMethyl;
            AmplitudeOfMembraneProtein = ScatteringLengthOfMembraneProtein * VolumeOfMembraneProtein * FormfactorOfMembraneProtein;

            AmplitudeOfDisc = AmplitudeOfBelt + AmplitudeOfCore + AmplitudeOfCaps + AmplitudeOfMethyl + AmplitudeOfMembraneProtein;
            PartialSumOfNanodisc = pow(cabs(AmplitudeOfDisc), 2);

            SumOverBeta += PartialSumOfNanodisc;
        }

        SumOverAlpha += SumOverBeta * sin(Alpha);
    }

    /// Rescale results
    RescaledSum = SumOverAlpha * AggregationNumber;
    ReturnValue = RescaledSum * (2.0 / pi) * BetaStepSize * AlphaStepSize;

    return ReturnValue;
}
