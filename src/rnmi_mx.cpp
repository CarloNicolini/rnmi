#include <iostream>
#include <string.h>
#include <sstream>
#include <igraph.h>
#include <sys/time.h>

#include "benchm.h"
#include "set_parameters.h" // the LFR parameters
#include "../paco/igraph_utils.h" // to handle conversion from EigenMatrix to igraph object and then to mxArray

#ifdef __linux__
#include <mex.h>
#endif

#ifdef __APPLE__
#include "mex.h"
#endif

using namespace std;

void printUsage()
{
    mexPrintf("rNMI: Relative Normalized Mutual Information\n");
    mexPrintf("This program computes the relative normalized mutual information of two partitions. Membership vectors starts from 0.\n");
    mexPrintf("rnmi_value = rnmi(memb1, memb2)\n");

enum error_type
{
    NO_ERROR = 0,
    ERROR_TOO_MANY_OUTPUT_ARGS = 1,
    ERROR_NOT_ENOUGH_ARGS = 2,
    ERROR_ARG_VALUE = 3,
    ERROR_ARG_TYPE = 4,
    ERROR_MATRIX = 5,
    ERROR_ARG_EMPTY=6,
    ERROR_ARG_UNKNOWN=7
};

static const char *error_strings[] =
{
    "",
    "Too many output arguments.",
    "Not enough input arguments.",
    "Non valid argument value.",
    "Non valid argument type.",
    "Non valid input adjacency matrix. Must be symmetric real dense-type square matrix.",
    "Expected some argument value but empty found.",
    "Unkwown argument."
};



/**
 * @brief mexFunction
 * @param nOutputArgs
 * @param outputArgs
 * @param nInputArgs
 * @param inputArgs
 */
void mexFunction(int nOutputArgs, mxArray *outputArgs[], int nInputArgs, const mxArray * inputArgs[])
{

    // Finish the function
    return;
}
