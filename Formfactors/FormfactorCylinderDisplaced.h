/// This function computes the formfactor for a cylinder of a given height displaced a certain length away from the origin.

complex FormfactorCylinderDisplaced(double q, double Radius, double HeightOfCylinder, double Displacement, double Alpha)
{
    // Declarations
    double Dummy1;
    double Dummy2;
    complex Dummy3;
    complex ReturnValue;

    // Computation
    Dummy1 = j1(q * Radius * sin(Alpha)) / (q * Radius * sin(Alpha));
    Dummy2 = sin(q * HeightOfCylinder * cos(Alpha) / 2.0) / (q * HeightOfCylinder * cos(Alpha) / 2.0);
    Dummy3 = cexp(_Complex_I * q * cos(Alpha) * Displacement);

    ReturnValue = 2.0 * Dummy1 * Dummy2 * Dummy3;

    return ReturnValue;
}
