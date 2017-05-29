double TripleShellLiposomeWMLandSQ(double q, double Radius, double ThicknessOfTails, double ThicknessOfHeads,
                                   double ScatteringLengthOfHeads, double ScatteringLengthOfTails, double Concentration, double lipidvolume)
{
    /// Declarations
    // Declarations used to describe the formfactors
    double RadiusOuterSphere;
    double RadiusInnerSphere;
    double RadiusOuterTails;
    double RadiusInnerTails;

    double PsiOuterSphere;
    double PsiInnerSphere;
    double PsiOuterTails;
    double PsiInnerTails;

    // Declare volume weights
    double ComputedVolumeOfOuterSphere;
    double ComputedVolumeOfInnerSphere;
    double ComputedVolumeOuterTails;
    double ComputedVolumeInnerTails;

    double TotalVolumeOfOuterHeads;
    double TotalVolumeOfInnerHeads;
    double TotalVolumeOfTails;

    // Declarations used to compute amplitudes
    double FormfactorOfOuterHeads;
    double FormfactorOfInnerHeads;
    double FormfactorOfOfTails;

    // Declarations of amplitudes
    double AmplitudeOfOuterHeads;
    double AmplitudeOfInnerHeads;
    double AmplitudeOfTails;

    double Intensity;
    double ReturnValue;
    double V_tot;
    double b_tot;
    double rho_tot;

    /// Begin computation
    // Compute formfactors of necessary spheres
    RadiusOuterSphere = Radius + ThicknessOfHeads + ThicknessOfTails;
    PsiOuterSphere = FormfactorSphere(q * RadiusOuterSphere);

    RadiusOuterTails = Radius + ThicknessOfTails;
    PsiOuterTails = FormfactorSphere(q * RadiusOuterTails);

    RadiusInnerTails = Radius - ThicknessOfTails ;
    PsiInnerTails = FormfactorSphere(q * RadiusInnerTails);

    RadiusInnerSphere = Radius - ThicknessOfHeads - ThicknessOfTails;
    PsiInnerSphere = FormfactorSphere(q * RadiusInnerSphere);

    // Derive volumes of these spheres
    ComputedVolumeOfOuterSphere = 4.0 / 3.0 * pi * pow(RadiusOuterSphere, 3);
    ComputedVolumeOuterTails = 4.0 / 3.0 * pi * pow(RadiusOuterTails, 3);
    ComputedVolumeInnerTails = 4.0 / 3.0 * pi * pow(RadiusInnerTails, 3);
    ComputedVolumeOfInnerSphere = 4.0 / 3.0 * pi * pow(RadiusInnerSphere, 3);

    TotalVolumeOfOuterHeads = ComputedVolumeOfOuterSphere - ComputedVolumeOuterTails;
    TotalVolumeOfTails = ComputedVolumeOuterTails - ComputedVolumeInnerTails;
    TotalVolumeOfInnerHeads = ComputedVolumeInnerTails - ComputedVolumeOfInnerSphere;

    V_tot = TotalVolumeOfOuterHeads + TotalVolumeOfTails + TotalVolumeOfInnerHeads;
    b_tot = ScatteringLengthOfTails * TotalVolumeOfTails + ScatteringLengthOfHeads * (TotalVolumeOfOuterHeads + TotalVolumeOfInnerHeads);
    rho_tot = b_tot / V_tot;

    // Compute formfactors of the different shells
    FormfactorOfOuterHeads = (ComputedVolumeOfOuterSphere * PsiOuterSphere - ComputedVolumeOuterTails * PsiOuterTails) /
                             (ComputedVolumeOfOuterSphere - ComputedVolumeOuterTails);

    FormfactorOfOfTails = (ComputedVolumeOuterTails * PsiOuterTails - ComputedVolumeInnerTails * PsiInnerTails) /
                          (ComputedVolumeOuterTails - ComputedVolumeInnerTails);

    FormfactorOfInnerHeads= (ComputedVolumeInnerTails * PsiInnerTails - ComputedVolumeOfInnerSphere * PsiInnerSphere) /
                            (ComputedVolumeInnerTails - ComputedVolumeOfInnerSphere);

    // Compute amplitudes of the different shells
    AmplitudeOfOuterHeads = TotalVolumeOfOuterHeads * FormfactorOfOuterHeads * ScatteringLengthOfHeads;
    AmplitudeOfTails = TotalVolumeOfTails * FormfactorOfOfTails * ScatteringLengthOfTails;
    AmplitudeOfInnerHeads = TotalVolumeOfInnerHeads * FormfactorOfInnerHeads * ScatteringLengthOfHeads;

    // Compute intensity of entire liposome
    Intensity = pow((AmplitudeOfOuterHeads + AmplitudeOfInnerHeads + AmplitudeOfTails)/(rho_tot*V_tot), 2);

    /// Rescale results
    Intensity *= pow(V_tot * rho_tot, 2);
    Concentration /= V_tot / lipidvolume;
    ReturnValue = Intensity * Concentration;

    return ReturnValue;
}
