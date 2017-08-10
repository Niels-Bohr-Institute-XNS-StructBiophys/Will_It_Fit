double MonodisperseMicelles(double q, double TotalRadius, double RadiusOfCore, double ScatteringLengthDensityOfHeads, double ScatteringLengthDensityOfTails,
                            double Concentration)
{
    /// Declarations
    // Declarations used to describe the formfactors
    double PsiOuterSphere;
    double PsiInnerSphere;

    // Declare volume weights
    double ComputedVolumeOfOuterSphere;
    double ComputedVolumeOfInnerSphere;
    double ComputedVolumeOfOuterLayer;

    // Declarations used to compute amplitudes
    double FormfactorOfOuterLayer;
    double FormfactorOfInnerSphere;

    // Declarations of amplitudes
    double AmplitudeOfOuterLayer;
    double AmplitudeOfInnerSphere;
    double Intensity;
    double ReturnValue;

    /// Begin computation
    // Compute formfactors of necessary spheres
    PsiOuterSphere = FormfactorSphere(q * TotalRadius);
    PsiInnerSphere = FormfactorSphere(q * RadiusOfCore);

    // Derive volumes of these spheres
    ComputedVolumeOfOuterSphere = 4.0 / 3.0 * pi * pow(TotalRadius, 3);
    ComputedVolumeOfInnerSphere = 4.0 / 3.0 * pi * pow(RadiusOfCore, 3);
    ComputedVolumeOfOuterLayer  = ComputedVolumeOfOuterSphere - ComputedVolumeOfInnerSphere;

    // Compute formfactors of the different shells
    FormfactorOfOuterLayer  = (ComputedVolumeOfOuterSphere * PsiOuterSphere - ComputedVolumeOfInnerSphere * PsiInnerSphere) /
                              (ComputedVolumeOfOuterSphere - ComputedVolumeOfInnerSphere);

    FormfactorOfInnerSphere = PsiInnerSphere;

    // Compute amplitudes of the different shells
    AmplitudeOfOuterLayer  = ComputedVolumeOfOuterLayer  * FormfactorOfOuterLayer  * ScatteringLengthDensityOfHeads;
    AmplitudeOfInnerSphere = ComputedVolumeOfInnerSphere * FormfactorOfInnerSphere * ScatteringLengthDensityOfTails;

    // Compute intensity of entire liposome
    Intensity = pow(AmplitudeOfOuterLayer + AmplitudeOfInnerSphere, 2);

    /// Rescale results
    ReturnValue = Intensity * Concentration;

    return ReturnValue;
}
