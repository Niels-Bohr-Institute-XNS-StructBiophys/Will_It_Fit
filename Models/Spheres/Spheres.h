double Spheres(double q, double Radius, double Concentration, double ScatteringLength)
{
    /// Declaration
    double Sum;
    double ReturnValue;
    double Volume;
    double Psi;
    double Amplitude;

    /// Begin computation
    Volume    = 4.0 / 3.0 * pi * pow(Radius, 3);
    Psi       = FormfactorSphere(q * Radius);
    Amplitude = ScatteringLength * Volume * Psi;
    Sum       = pow(Amplitude, 2);

    /// Rescale results
    ReturnValue = Sum * Concentration;

    return ReturnValue;
}
