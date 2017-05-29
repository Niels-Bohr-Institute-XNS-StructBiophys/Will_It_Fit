double AnisotropicMicelles(double q, double MinorRadiusOfCore, double MajorRadiusOfCore, double ScatteringLengthDensityOfHeads, double ScatteringLengthDensityOfTails,
                           double Concentration, double ThicknessOfShell, double ShellRoughness, double CoreRoughness, double AxisRatioOfCore, double AxisRatioOfShell)
{
    /// Declarations
    // Dummy variable
    int i;

    // Declarations used to describe the formfactors
    double PsiOuter;
    double PsiInner;

    double ProjectedRadiusOfCore;
    double ProjectedTotalRadius;

    // Declare volume weights
    double ComputedVolumeOfOuter;
    double ComputedVolumeOfInner;
    double ComputedVolumeOfShell;

    // Declarations used to compute amplitudes
    double FormfactorOfShell;
    double FormfactorOfInner;

    // Declarations of amplitudes
    double AmplitudeOfShell;
    double AmplitudeOfInner;
    double Intensity;

    // Introduce the parameters used in the integrations
    const int NumberOfStepsInAlpha = 50;
    const double AlphaStepSize = pi / (2.0 * NumberOfStepsInAlpha);
    double Alpha;
    double SumOverAlpha = 0.0;

    /// Begin computation
    // Derive volumes of these spheres
    ComputedVolumeOfInner = 4.0 / 3.0 * pi * AxisRatioOfCore * pow(MinorRadiusOfCore, 3);
    ComputedVolumeOfOuter = 4.0 / 3.0 * pi * AxisRatioOfShell * pow(MinorRadiusOfCore + ThicknessOfShell, 3);
    ComputedVolumeOfShell = ComputedVolumeOfOuter - ComputedVolumeOfInner;

    for (i = 1; i <= NumberOfStepsInAlpha; ++i) {
        Alpha = AlphaStepSize * i;

        // Derive projected radii
        ProjectedRadiusOfCore = MinorRadiusOfCore * sqrt(pow(sin(Alpha), 2) + pow(AxisRatioOfCore * cos(Alpha), 2));
        ProjectedTotalRadius  = (MinorRadiusOfCore + ThicknessOfShell) * sqrt(pow(sin(Alpha), 2) + pow(AxisRatioOfShell * cos(Alpha), 2));

        // Compute relevant formfactors
        PsiInner = FormfactorSphere(q * ProjectedRadiusOfCore);
        PsiOuter = FormfactorSphere(q * ProjectedTotalRadius);

        // Compute remaining formfactor
        FormfactorOfShell = (ComputedVolumeOfOuter * PsiOuter - ComputedVolumeOfInner * PsiInner) /
                            (ComputedVolumeOfOuter - ComputedVolumeOfInner);

        FormfactorOfInner = PsiInner;

        // Derive amplitudes
        AmplitudeOfShell = ComputedVolumeOfShell * FormfactorOfShell *
                           ScatteringLengthDensityOfHeads * exp(- 0.5 * pow(q * ShellRoughness, 2));

        AmplitudeOfInner = ComputedVolumeOfInner * FormfactorOfInner *
                           ScatteringLengthDensityOfTails * exp(- 0.5 * pow(q * CoreRoughness, 2));

        // Compute intensity of entire micelle
        Intensity = pow(AmplitudeOfShell + AmplitudeOfInner, 2);

        // Conclude integral step
        SumOverAlpha += Intensity * sin(Alpha);
    }

    return SumOverAlpha * Concentration * AlphaStepSize;
}
