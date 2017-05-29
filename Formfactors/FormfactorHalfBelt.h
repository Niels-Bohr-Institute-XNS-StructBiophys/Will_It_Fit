/// This function computes the formfactor of a cylinder with a hole in the middle from AngularOffset to pi + AngularOffset.

complex FormfactorHalfBelt(double q, double Alpha, double Beta, double Height, double AngularOffset, double MinorAxis, double MajorAxis, double Width)
{
    /// Declarations
    // Dummy variables
    int i;

    complex Dummy1;
    complex Dummy2;

    // Variables to be used to compute the integral over phi
    double Phi;
    double PhiMin = AngularOffset;
    double PhiMax = AngularOffset + pi;
    const int NumberOfStepsInPhi = 200;
    double PhiStep = (PhiMax - PhiMin) / NumberOfStepsInPhi;
    complex SumOverPhi;

    double InnerRadius;
    double OuterRadius;

    // Partial sums
    complex ReturnValue;
    double VerticalSum;

    // Properties of structure
    double Volume = (pi * Height * (MajorAxis + Width) * (MinorAxis + Width) -
                    pi * Height * MajorAxis * MinorAxis) / 2.0;

    /// Computation
    VerticalSum = 2.0 / (q * cos(Alpha)) * sin(q * cos(Alpha) * Height / 2);

    // Begin integration over phi
    SumOverPhi = 0.0;

    for (i = 0; i < NumberOfStepsInPhi; ++i) {
        Phi = PhiMin + i * PhiStep;

        InnerRadius = MajorAxis * MinorAxis /
                      sqrt(pow(MajorAxis * sin(Phi), 2) + pow(MinorAxis * cos(Phi), 2));

        OuterRadius = (MajorAxis + Width) * (MinorAxis + Width) /
                      sqrt(pow((MajorAxis + Width) * sin(Phi), 2) + pow((MinorAxis + Width) * cos(Phi), 2));

        Dummy1 = _Complex_I * q * (cos(Beta) * sin(Alpha) * cos(Phi) + sin(Beta) * sin(Alpha) * sin(Phi));

        Dummy2 = 1.0 / Dummy1 * cexp(Dummy1 * OuterRadius) * (OuterRadius - 1.0 / Dummy1) -
                 1.0 / Dummy1 * cexp(Dummy1 * InnerRadius) * (InnerRadius - 1.0 / Dummy1);

        SumOverPhi += Dummy2 * PhiStep;
    }

    // Conclusion
    ReturnValue = VerticalSum * SumOverPhi / Volume;

    return ReturnValue;
}

