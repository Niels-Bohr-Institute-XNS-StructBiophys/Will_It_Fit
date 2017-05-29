/************************************************
 *
 * Model name: PDB
 * Author    : Martin Cramer Pedersen
 *             Based on code by Søren Kynde
 * Email     : martin.pedersen@anu.edu.au
 *
 ***********************************************/

// List of headers to be included on compilation
#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_sf_bessel.h>

//Order of harmonics
#define NumberOfHarmonics 17

// More headers
#include "MathAndHelpfunctions.h"
#include "ComputeConstraints.h"
#include "Model.h"
#include "OutputData.h"
