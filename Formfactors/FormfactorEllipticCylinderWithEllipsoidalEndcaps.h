/// This function is used to compute the structure of a cylinder with triaxial half ellipsoids as endcaps
double FormfactorCylinderWithEndcaps(double q, double Alpha, double Beta, double RadiusOne, double RadiusTwo, double Height,
                                     double ScaleFactorOfEndcaps, double VerticalAxisOfEndcaps)
{
    /// Declarations
    // Dummies
    int i;

    double Dummy1;
    double Dummy2;
    double Dummy3;

    double EquivalentRadiusInCylinder;
    double EquivalentRadiusInEndcaps;

    const int TSteps = 50;
    double TMin;
    double TMax;
    double TStepSize;
    double T;

    double ReturnValue;
    double SumOverT = 0.0;

    // Properties
    double RadiusOneOfEllipsoid = RadiusOne * ScaleFactorOfEndcaps;
    double RadiusTwoOfEllipsoid = RadiusTwo * ScaleFactorOfEndcaps;
    double ShiftOfEndcaps;
    double Volume;

    // Caps center inside cylidner
    ShiftOfEndcaps = - VerticalAxisOfEndcaps / RadiusTwoOfEllipsoid * sqrt(pow(RadiusTwoOfEllipsoid, 2) - pow(RadiusTwo, 2));

    // Compute volume
    Volume = 2 * 3.14159 * RadiusTwoOfEllipsoid * RadiusOneOfEllipsoid *
             (2.0 / 3.0 * VerticalAxisOfEndcaps - fabs(ShiftOfEndcaps) + pow(fabs(ShiftOfEndcaps), 3) / (3.0 * pow(VerticalAxisOfEndcaps, 2))) +
             3.14156 * RadiusTwo * RadiusOne * Height;

    // Compute amplitude
    EquivalentRadiusInCylinder = sqrt(pow(RadiusOne * cos(Beta), 2) + pow(RadiusTwo * sin(Beta), 2));
    EquivalentRadiusInEndcaps  = sqrt(pow(RadiusOneOfEllipsoid * cos(Beta), 2) + pow(RadiusTwoOfEllipsoid * sin(Beta), 2));

    TMin = - ShiftOfEndcaps / VerticalAxisOfEndcaps;
    TMax = 1.0;
    TStepSize = (TMax - TMin) / TSteps;

    Dummy1 = pi * RadiusOne * RadiusTwo * Height *
             2.0 * sin(q * Height / 2.0 * cos(Alpha)) / (q * Height / 2.0 * cos(Alpha)) *
             j1(q * EquivalentRadiusInCylinder * sin(Alpha)) / (q * EquivalentRadiusInCylinder * sin(Alpha));

    SumOverT = 0.0;

    for (i = 0; i < TSteps; ++i) {
        T = i * TStepSize + TMin;

        Dummy2 = cos(q * cos(Alpha) * (VerticalAxisOfEndcaps * T + ShiftOfEndcaps + Height / 2.0)) *
        (1 - pow(T, 2)) * j1(q * EquivalentRadiusInEndcaps * sin(Alpha) * sqrt(1.0 - pow(T, 2))) /
        (q * EquivalentRadiusInEndcaps * sin(Alpha) * sqrt(1.0 - pow(T, 2)));

        SumOverT += Dummy2 * TStepSize;
    }

    Dummy3 = 4.0 * pi * RadiusOneOfEllipsoid * RadiusTwoOfEllipsoid * VerticalAxisOfEndcaps * SumOverT;

    // Conclusion
    ReturnValue = (Dummy1 + Dummy3) / Volume;

    return ReturnValue;
}

