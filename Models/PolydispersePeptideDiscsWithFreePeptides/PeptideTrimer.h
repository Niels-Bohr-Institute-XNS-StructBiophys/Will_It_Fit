double PeptideTrimer(double q, double Radius, double Height, double ScatteringLengthDensity, double ConcentrationOfTrimers)
{
    /// Declarations
    // Dummies
    int i;
    int j;

    double Dummy1;
    double Dummy2;

    // Amplitudes and phases
    complex Phase1;
    complex Phase2;
    complex Phase3;

    complex Amplitude1;
    complex Amplitude2;
    complex Amplitude3;

    double Volume = pi * pow(Radius, 2) * Height;

    // Variables used in spherical average
    const int AlphaSteps = 50;
    double AlphaMin;
    double AlphaMax;
    double AlphaStepSize;
    double Alpha;

    const int BetaSteps = 50;
    double BetaMin;
    double BetaMax;
    double BetaStepSize;
    double Beta;

    // Sums used in integration
    double ReturnValue;
    double SumOverAlpha = 0.0;
    double SumOverBeta = 0.0;

    /// Integration
    // Begin integration over alpha
    AlphaMin = 0.0;
    AlphaMax = pi;
    AlphaStepSize = (AlphaMax - AlphaMin) / AlphaSteps;

    BetaMin = 0.0;
    BetaMax = 2.0 * pi;
    BetaStepSize = (BetaMax - BetaMin) / BetaSteps;

    for (i = 1; i <= AlphaSteps; ++i) {
        Alpha = i * AlphaStepSize + AlphaMin;

        SumOverBeta = 0.0;

        for (j = 1; j <= BetaSteps; ++j) {
            Beta = j * BetaStepSize + BetaMin;

            Dummy1 = 2.0 * j1(q * Radius * sin(Alpha)) / (q * Radius * sin(Alpha));
            Dummy2 = sin(cos(Alpha) * q * Height / 2.0) / (q * cos(Alpha) * Height / 2.0);

            Phase1 = 1.0;
            Phase2 = cexp(_Complex_I * q * 2 * Radius * sin(Alpha) * cos(Beta));
            Phase3 = cexp(_Complex_I * q * (2 * Radius * sin(Alpha) * cos(Beta) * cos(pi / 3.0) +
                                            2 * Radius * sin(Alpha) * sin(Beta) * sin(pi / 3.0)));

            Amplitude1 = Dummy1 * Dummy2 * Phase1 * ScatteringLengthDensity * Volume;
            Amplitude2 = Dummy1 * Dummy2 * Phase2 * ScatteringLengthDensity * Volume;
            Amplitude3 = Dummy1 * Dummy2 * Phase3 * ScatteringLengthDensity * Volume;

            SumOverBeta += pow(cabs(Amplitude1 + Amplitude2 + Amplitude3), 2);
        }

        SumOverAlpha += SumOverBeta * sin(Alpha);
    }

    /// Rescale and return
    SumOverAlpha *= ConcentrationOfTrimers;
    ReturnValue = 1.0 / (4.0 * pi) * SumOverAlpha * AlphaStepSize * BetaStepSize;

    return ReturnValue;
}
