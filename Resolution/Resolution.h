/*
 * For a detailed description of the manner in which, the
 * resolution effects have been implemented, please consult:
 *
 * Pedersen, Posselt, and Mortensen (1990)
 * J. Appl. Cryst 23, 321-333
 *
 * for details. Note that this function and the associated
 * ones computes a number of unused quantities (which could
 * be used to compute 2D-smearing.
 */

 // The modified Bessel-function of the first kind of order zero
double ModifiedBesselFunctionI0(double q, double DummyQ, double SigmaQ)
{
    double x;
    double y;

    double Weight1;
    double Weight2;
    double Result;

    const double P1 = 1.0;
    const double P2 = 3.5156229;
    const double P3 = 3.0899424;
    const double P4 = 1.2067492;
    const double P5 = 0.2659732;
    const double P6 = 0.360768e-1;
    const double P7 = 0.45813e-2;

    const double Q1 = 0.39894228;
    const double Q2 = 0.1328592e-1;
    const double Q3 = 0.225319e-2;
    const double Q4 = -0.157565e-2;
    const double Q5 = 0.916281e-2;
    const double Q6 = -0.2057706e-1;
    const double Q7 = 0.2635537e-1;
    const double Q8 = -0.1647633e-1;
    const double Q9 = 0.392377e-2;

    x       = q * DummyQ / pow(SigmaQ, 2);
    Weight1 = DummyQ / pow(SigmaQ, 2);

    if (fabs(x) < 3.75) {
        y       = pow(x / 3.75, 2);
        Weight2 = exp(- (pow(q, 2) + pow(DummyQ, 2)) / (2.0 * pow(SigmaQ, 2)));
        Result  = P1 + y * (P2 + y * (P3 + y * (P4 + y * (P5 + y * (P6 + y * P7)))));
    } else {
        y       = 3.75 / fabs(x);
        Weight2 = exp(- (pow(q, 2) + pow(DummyQ, 2)) / (2.0 * pow(SigmaQ, 2)) + x) / sqrt(fabs(x));
        Result  = Q1 + y * (Q2 + y * (Q3 + y * (Q4 + y * (Q5 + y * (Q6 + y * (Q7 + y * (Q8 + y * Q9)))))));
    }

    return Weight1 * Result * Weight2;
}

// Function computing resolution properties
int Resolution(struct Dataset * Data, int NumberOfSpectra, char ResolutionFile[256], int NumberOfSmearingFolds)
{
    /// Declarations
    // Dummy variables
    double DummyQ;
    double q;
    double Dummy1;
    double Dummy2;
    double Dummy3;
    double Dummy4;
    double Dummy5;
    double Dummy6;
    double Dummy7;
    double Dummy8;

    int DummyInt;
    int i;
    int j;
    int k;

    char Line[256];

    // Used in computation
    double Stepsize;
    double SigmaQ;
    double SumOfWeights;

    // File info
    FILE *Inputfile;

    // Resolution information
    double SampleToDetectorDistance[NumberOfSpectra];
    double CollimationLength[NumberOfSpectra];

    double RadiusOfFirstPinhole[NumberOfSpectra];
    double RadiusOfSecondPinhole[NumberOfSpectra];

    double Wavelength[NumberOfSpectra];
    double WavelengthSpread[NumberOfSpectra];

    double DetectorWidth[NumberOfSpectra];
    double BinningWidth[NumberOfSpectra];


    /// Obtain information from resolution info file
    // Open file
    Inputfile = fopen(ResolutionFile, "r");

    if (Inputfile == NULL){
        printf("An error occured when attempting to open resolution file.");
        return -1;
    }

    // Fetch values from file
    i = 0;

    printf("\n");
    printf(".res-file found! \n");

    while (fgets(Line, sizeof(Line), Inputfile) != NULL) {

        if (sscanf(Line, "%d %lf %lf %lf %lf %lf %lf %lf %lf", &DummyInt, &Dummy1, &Dummy2, &Dummy3, &Dummy4, &Dummy5, &Dummy6, &Dummy7, &Dummy8) == 9) {

            if (DummyInt != 0) {
                Data[i].IncludeResolutionEffects = true;
                printf("Including instrumental smearing in fits of dataset %d. \n", i + 1);
            }

            RadiusOfFirstPinhole[i]     = Dummy1;
            RadiusOfSecondPinhole[i]    = Dummy2;
            CollimationLength[i]        = Dummy3;
            SampleToDetectorDistance[i] = Dummy4;
            Wavelength[i]               = Dummy5;
            WavelengthSpread[i]         = Dummy6;
            DetectorWidth[i]            = Dummy7;
            BinningWidth[i]             = Dummy8;

            ++i;
        }
    }

    if (i == 0) {
        return -1;
    }

    if (i != NumberOfSpectra) {
        printf("The number of spectra in the .card-file does not match the number of spectra in the .res-file.");
        return -1;
    }

    // Close file
    fclose(Inputfile);

    /// Computations
    // Loop over all datasets
    for (i = 0; i < NumberOfSpectra; ++i) {

        if (Data[i].IncludeResolutionEffects == true) {

            // Loop over all points in these files
            for (j = Data[i].NMin; j < Data[i].NMax; ++j) {
                q = Data[i].QValues[j];

                SigmaQ = SigmaOfQ(q, RadiusOfFirstPinhole[i], RadiusOfSecondPinhole[i], CollimationLength[i], SampleToDetectorDistance[i], Wavelength[i],
                                  WavelengthSpread[i], DetectorWidth[i], BinningWidth[i]);

                Data[i].SigmaQValues[j] = SigmaQ;
                Stepsize  = 6.0 * SigmaQ / (1.0 * NumberOfSmearingFolds);
                SumOfWeights = 0.0;

                for (k = 0; k < NumberOfSmearingFolds; ++k) {
                    DummyQ = q + (k + 0.5 - NumberOfSmearingFolds / 2.0) * Stepsize;

                    if (DummyQ < 1e-5) {
                        DummyQ = 1e-5;
                    }

                    Data[i].ResolutionWeights[j][k] = ModifiedBesselFunctionI0(q, DummyQ, SigmaQ);

                    SumOfWeights += Data[i].ResolutionWeights[j][k];
                }

                for (k = 0; k < NumberOfSmearingFolds; ++k) {
                    Data[i].ResolutionWeights[j][k] /= SumOfWeights;
                }
            }
        }
    }

    return 0;
}
