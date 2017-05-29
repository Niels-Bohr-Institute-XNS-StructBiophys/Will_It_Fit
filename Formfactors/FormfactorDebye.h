double FormfactorDebye(double q, double r)
{
    // Variables used in functions
    double ReturnValue;

    // Computation
    ReturnValue = 2 * ((exp(-(pow(q, 2) * pow(r, 2))) + (pow(q, 2) * pow(r, 2)) - 1) /
                       (pow(q, 4) * pow(r, 4)));

    return ReturnValue;
}
