double FormfactorRodWithoutEndcaps(double q, double r, double l, double alpha)
{
    // Variables needed in function
    double ReturnValue;

    // Computation
    ReturnValue = j0(q * r * sin(alpha)) * sin(q * l * cos(alpha) / 2.f) / (q * l * cos(alpha) / 2.f);

    return ReturnValue;
}
