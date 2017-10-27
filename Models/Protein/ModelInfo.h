/************************************************
 *
 * Model name: Protein With HetAtoms
 * Author    : Søren Kynde and Nicholas Skar-Gislinge
 * Email     : kynde@nbi.ku.dk
 *
 ***********************************************/

 // List of headers to be included on compilation
#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_sf_bessel.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

 //Order of harmonics
#define NumberOfHarmonics 17
#include "MathAndHelpfunctions.h"
#include "Constraints.h"
//#include "SemiflexibleChain.h"
#include "Model.h"
#include "OutputData.h"
#include "OutputDataBash.h"

