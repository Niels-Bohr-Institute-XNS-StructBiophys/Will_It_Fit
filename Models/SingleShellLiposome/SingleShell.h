double SingleShell(double q, double Radius, double Thickness, double ScatteringLength, double Concentration)
{
    /// Declarations
    // Declarations used to describe the formfactors
    double RadiusOuter;
    double RadiusInner;

    double PsiOuter;
    double PsiInner;

    // Declare volume weights
    double ComputedVolumeOfOuterSphere;
    double ComputedVolumeOfInnerSphere;

    double TotalVolume;

    // Declarations used to compute amplitudes
    double Formfactor;

    // Declarations of amplitudes
    double Amplitude;

    double Intensity;
    double ReturnValue;

    /// Begin computation
    // Compute formfactors of necessary spheres
    RadiusOuter = Radius + 0.5 * Thickness;
    PsiOuter = FormfactorSphere(q * RadiusOuter);

    RadiusInner = Radius - 0.5 * Thickness;
    PsiInner = FormfactorSphere(q * RadiusInner);

    // Derive volumes of these spheres
    ComputedVolumeOfOuterSphere = 4.0 / 3.0 * pi * pow(RadiusOuter, 3);
    ComputedVolumeOfInnerSphere = 4.0 / 3.0 * pi * pow(RadiusInner, 3);

    TotalVolume = ComputedVolumeOfOuterSphere - ComputedVolumeOfInnerSphere;

    // Compute formfactors of the different shells
    Formfactor = (ComputedVolumeOfOuterSphere * PsiOuter - ComputedVolumeOfInnerSphere * PsiInner) /
                 (ComputedVolumeOfOuterSphere - ComputedVolumeOfInnerSphere);

    // Compute amplitudes of the different shells
    Amplitude = TotalVolume * Formfactor * ScatteringLength;

    // Compute intensity of entire liposome
    Intensity = pow(Amplitude, 2);

    /// Rescale results
    ReturnValue = Intensity * Concentration;

    return ReturnValue;
}
