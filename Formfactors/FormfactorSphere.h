/// This function computes the for factor for a sphere

double FormfactorSphere(double qR)
{
    // Variables used in function
    double ReturnValue;

    // Computation
    ReturnValue = (3.0 * (sin(qR) - qR * cos(qR))) / pow(qR, 3);

    return ReturnValue;
}
