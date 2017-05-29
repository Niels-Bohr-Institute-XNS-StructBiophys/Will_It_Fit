/************************************************
 *
 * Model name: Membrane Protein in nanodisc with curved bilayer
 * Author    : Søren Kynde
 * Email     : kynde@nbi.dk
 *
 ***********************************************/

// List of headers to be included on compilation

#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_sf_bessel.h>
#define Nh 17 //Order of harmonics
#define skip 2 //Skip uneaven harmonics (Because of mirror symmetry of disc)
#define Nphi ((Nh+1)*2)
#define Ntheta ((Nh+1)*2)

#include "MathAndHelpfunctions.h"
#include "ComputeConstraints.h"
#include "Model.h"
#include "OutputData.h"
