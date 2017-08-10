// Hard Sphere Structure Factor see p.410 Lindner and Zemb, "Neutrons, X-rays and light: Scattering Methods Applied to Soft Condenced Matter" 2002 North-Holland. (the "Bombannes book")
double LiposomeSQ(double q, double VolumeFraction, double Radius, double cRadiusHardSphere)
{
    double RadiusHardSphere = Radius*cRadiusHardSphere;
    double A = 2*q*RadiusHardSphere;
    double alpha = pow(1+2*VolumeFraction, 2)/pow(1-VolumeFraction,4);
    double beta = -6*VolumeFraction*pow(1+VolumeFraction/2, 2)/pow(1-VolumeFraction,4);
    double gamma = VolumeFraction*alpha/2;
    double GA = alpha*(sin(A)-A*cos(A))/pow(A, 2)+
    beta*(2*A*sin(A)+(2-pow(A, 2))*cos(A)-2)/pow(A, 3)+
    gamma*(-pow(A, 4)*cos(A)+4*((3*pow(A, 2)-6)*cos(A)+(pow(A, 3)-6*A)*sin(A)+6))/pow(A, 5);

    double ReturnValue;

    ReturnValue = 1/(1+24*VolumeFraction*GA/A);

    return ReturnValue;

}

double Model(double q, double * Parameters, double * Constraints, double Contrast, struct Protein ProteinStructure, struct UserDefined * UserDefinedStructure)
{
    /// Declarations
    // Parameters for model
    double Dummy;
    double Dummy1;
    double Dummy2;
    int i;

    // Variables describing the sample
    double Roughness = 0.0;

    double ScatteringLengthOfTails;
    double ScatteringLengthOfHeads;

    double ThicknessOfHeads;
    double ThicknessOfTails;

    double Radius;
    double ConcentrationValue;

    // Variables describing the polydispersity
    double RMin;
    double RMax;
    double RStep;

    const int NumberOfSteps = 200;

    double DistributionOfRadii[NumberOfSteps];
    double DistributionOfRadiiWeight[NumberOfSteps];
    double TotalDistribution;

    // Declarations of partial sums
    double Intensity = 0.0;
    double Scaling;
    double Background;
    double SQ;
    double lipidvolume;

    /// Get values from constraints
    ScatteringLengthOfHeads = Constraints[0];
    ScatteringLengthOfTails = Constraints[1];
    ThicknessOfHeads = Constraints[2];
    ThicknessOfTails = Constraints[3];
    ConcentrationValue = Constraints[4];
    Radius = Constraints[5];
    lipidvolume = Constraints[6] + Constraints[7];

    ConcentrationValue = Parameters[18] * ConcentrationValue;

    /// Compute

    /// Begin computation
    // Create gaussian distribution of radii
    RMin = Radius - 3.0 * Radius * Parameters[7];
    RMax = Radius + 3.0 * Radius * Parameters[7];

    if (RMin < 0.0) {
        RMin = 0.0;
    }

    RStep = (RMax - RMin) / (NumberOfSteps - 1.0);

    TotalDistribution = 0;

    for (i = 0; i < NumberOfSteps; ++i) {
        DistributionOfRadii[i] = RMin + i * RStep;
        Dummy1 = 0.5f * pow((DistributionOfRadii[i] - Radius) / (sqrt(2) * Radius * Parameters[7]), 2);

        if (Dummy1 > 40) {
            Dummy2 = 0;
        } else {
            Dummy2 = 1 / sqrt(2 * pi * pow(Radius * Parameters[7], 2)) * exp(-Dummy1);
        }

        DistributionOfRadiiWeight[i] = Dummy2 * RStep;

        TotalDistribution += DistributionOfRadiiWeight[i];
    }

    for (i = 0; i < NumberOfSteps; ++i) {
        Intensity += DistributionOfRadiiWeight[i] / TotalDistribution * TripleShellLiposomeWMLandSQ(q, DistributionOfRadii[i], ThicknessOfTails, ThicknessOfHeads,
                                                                                              ScatteringLengthOfHeads, ScatteringLengthOfTails,
                                                                                              ConcentrationValue, lipidvolume);
    }

    SQ = LiposomeSQ(q, Parameters[20], Radius, Parameters[21]);

    /// Rescale result and return
    double SecondTerm = fabs(Parameters[14]) * exp(- pow(q - Parameters[15], 2) / (2.0 * pow(Parameters[16], 2))) +
                        Parameters[22] * fabs(Parameters[14]) * exp(- pow(q - 2.0 * Parameters[15], 2) / (2.0 * pow(Parameters[16], 2)));

    if (Contrast >= 0.0 && Contrast <= 100.0) {
        Roughness  = exp(-pow(q * Parameters[13], 2));
        Scaling    = Parameters[18] * (Contrast / 100.0 * Parameters[0] + (100.0 - Contrast) / 100.0 * Parameters[1]);
        Background = Contrast / 100.0 * Parameters[8] + (100.0 - Contrast) / 100.0 * Parameters[9];
    } else {
        Roughness  = exp(-pow(q * Parameters[12], 2));
        Scaling    = Parameters[18] * Parameters[2];
        Background = Parameters[10];
    }

    Intensity = Intensity * Scaling * Roughness * SQ + SecondTerm - Background;

    return Intensity;
}
