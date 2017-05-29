double FormfactorHammouda(double x)
{
    // Variables used in function
    double ReturnValue;

    // Computation
    ReturnValue = (1.0 - exp(-pow(x, 2))) / pow(x, 2);

    return ReturnValue;
}
