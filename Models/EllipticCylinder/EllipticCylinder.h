double EllipticCylinder(double q, double MinorRadius, double MajorRadius, double Height, double Concentration, double Volume, double ScatteringLength)
{
    // Declarations of dummy variables used in computations
    double Dummy1;
    double Dummy2;
    int i;
    int j;

    // Declarations used to compute return value
    double SumOverAlpha;
    double Sum;
    double RescaledSum;
    double ReturnValue;
    double SumOverBeta;

    double Psi;
    double Amplitude;
    double ProjectedRadius;

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
            Dummy1 = MinorRadius * cos(Beta);
            Dummy2 = MajorRadius * sin(Beta);

            ProjectedRadius = sqrt(pow(Dummy1, 2) + pow(Dummy2, 2));

            Psi = FormfactorCylinder(q, ProjectedRadius, Height, Alpha);

            Amplitude = ScatteringLength * Volume * Psi;

            Sum = pow(Amplitude, 2);

            SumOverBeta += Sum;
        }

        SumOverAlpha += SumOverBeta * sin(Alpha);
    }


    /// Rescale results
    ReturnValue = SumOverAlpha * Concentration * (2.0 / pi) * BetaStepSize * AlphaStepSize;

    return ReturnValue;
}
