double SigmaOfQ(double q, double RadiusOfFirstPinhole, double RadiusOfSecondPinhole, double CollimationLength, double SampleToDetectorDistance, double Wavelength,
                double WavelengthWidth, double DetectorWidth, double BinningWidth)
{
    /// Declarations
    // Physical quantities
    double Wavevector;
    double ScatteringAngle;

    // Computed deviations
    double SigmaDetector1;
    double SigmaDetector2;
    double SigmaWavelength;
    double SigmaCollimation1;
    double SigmaCollimation2;
    double SigmaBinning;

    // Dummies
    double AngularExtent1;
    double AngularExtent2;

    // Other
    double Beta1;
    double Beta2;
    double ReturnValue;
    double ConvertFWHMToSigma = 1.0 / (2.0 * sqrt(2.0 * log(2.0)));

    /// Computation
    // Physics (Scattering angle is ONE(!) theta)
    Wavevector = 2.0 * pi / Wavelength;
    ScatteringAngle = asin(q / (2.0 * Wavevector));

    // Detector
    SigmaDetector1 = Wavevector * cos(ScatteringAngle) * pow(cos(2 * ScatteringAngle), 2) * DetectorWidth / SampleToDetectorDistance * ConvertFWHMToSigma;
    SigmaDetector2 = Wavevector * cos(2 * ScatteringAngle) * DetectorWidth / SampleToDetectorDistance * ConvertFWHMToSigma;

    // Wavelength
    SigmaWavelength = q * WavelengthWidth * ConvertFWHMToSigma;

    // Collimation part 1
    AngularExtent1 = RadiusOfFirstPinhole / (CollimationLength + SampleToDetectorDistance / pow(cos(2.0 * ScatteringAngle), 2));
    AngularExtent2 = RadiusOfSecondPinhole * pow(cos(2.0 * ScatteringAngle), 2) / SampleToDetectorDistance;

    if (AngularExtent1 >= AngularExtent2) {
        Beta1 = 2.0 * RadiusOfFirstPinhole / CollimationLength - 0.5 * pow(RadiusOfSecondPinhole, 2) * pow(cos(2.0 * ScatteringAngle), 4) *
                pow(CollimationLength + SampleToDetectorDistance / pow(cos(2.0 * ScatteringAngle), 2), 2) /
                (RadiusOfFirstPinhole * pow(SampleToDetectorDistance, 2) * CollimationLength);
    } else {
        Beta1 = 2.0 * RadiusOfSecondPinhole * (1.0 / CollimationLength + pow(cos(2.0 * ScatteringAngle), 2) / SampleToDetectorDistance) -
                0.5 * pow(RadiusOfFirstPinhole, 2) * SampleToDetectorDistance / (RadiusOfSecondPinhole * CollimationLength *
                pow(cos(2.0 * ScatteringAngle), 2) * (CollimationLength + SampleToDetectorDistance / pow(cos(2.0 * ScatteringAngle), 2)));
    }

    // Collimation part 2
    AngularExtent1 = RadiusOfFirstPinhole / (CollimationLength + SampleToDetectorDistance / cos(2.0 * ScatteringAngle));
    AngularExtent2 = RadiusOfSecondPinhole * cos(2.0 * ScatteringAngle) / SampleToDetectorDistance;

    if (AngularExtent1 >= AngularExtent2) {
        Beta2 = 2.0 * RadiusOfFirstPinhole / CollimationLength - 0.5 * pow(RadiusOfSecondPinhole, 2) * pow(cos(2.0 * ScatteringAngle), 2) *
                pow(CollimationLength + SampleToDetectorDistance / cos(2.0 * ScatteringAngle), 2) /
                (RadiusOfFirstPinhole * pow(SampleToDetectorDistance, 2) * CollimationLength);
    } else {
        Beta2 = 2.0 * RadiusOfSecondPinhole * (1.0 / CollimationLength + cos(2.0 * ScatteringAngle) / SampleToDetectorDistance) -
                0.5 * pow(RadiusOfFirstPinhole, 2) * SampleToDetectorDistance / (RadiusOfSecondPinhole * CollimationLength *
                cos(2.0 * ScatteringAngle) * (CollimationLength + SampleToDetectorDistance / cos(2.0 * ScatteringAngle)));
    }

    // Collimation
    SigmaCollimation1 = Wavevector * cos(ScatteringAngle) * Beta1 * ConvertFWHMToSigma;
    SigmaCollimation2 = Wavevector * Beta2 * ConvertFWHMToSigma;

    // Binning kernel
    SigmaBinning = BinningWidth * ConvertFWHMToSigma;

    /// Combine results and return total deviation
    // For 1D-data
    ReturnValue = sqrt(pow(SigmaCollimation1, 2) + pow(SigmaDetector1, 2) + pow(SigmaWavelength, 2) + pow(SigmaBinning, 2));

    return ReturnValue;
}
