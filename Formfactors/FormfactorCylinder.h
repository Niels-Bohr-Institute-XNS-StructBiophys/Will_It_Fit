/// This function is used to compute the formfactor of a cylinder.

double FormfactorCylinder(double q, double Radius, double Height, double Alpha)
{
    // Variables used in function
    double ReturnValue;

    // Computation
    ReturnValue = 2.0 * j1(q * Radius * sin(Alpha)) / (q * Radius * sin(Alpha)) * sin(q * Height * cos(Alpha) / 2.0) / (q * Height * cos(Alpha) / 2.0);

    return ReturnValue;
}
