#include <iostream>

#ifdef __linux__
#include <mex.h>
#endif

#ifdef __APPLE__
#include "mex.h"
#endif

#include <set>
#include "nmi.h"

using namespace std;

void printUsage()
{
    mexPrintf("rNMI: Relative Normalized Mutual Information\n");
    mexPrintf("This program computes the relative normalized mutual information of two partitions. Membership vectors starts from 0.\n");
    mexPrintf("rnmi_value = rnmi(memb1, memb2)\n");
}

/**
 * @brief mexFunction
 * @param nOutputArgs
 * @param outputArgs
 * @param nInputArgs
 * @param inputArgs
 */
void mexFunction(int nOutputArgs, mxArray *outputArgs[], int nInputArgs, const mxArray * inputArgs[])
{
    int randseed=int(time(NULL));
    ZRANDOMv3 rg(randseed);

    const mxArray* p1 = inputArgs[0];
    const mxArray* p2 = inputArgs[1];

    const size_t n1 = mxGetNumberOfElements(p1);
    const size_t n2 = mxGetNumberOfElements(p2);

    if (n1!=n2)
    {
        mexErrMsgTxt("The two partitions have different number of elements.");
        return;
    }

    vector<double> memb1(n1);
    vector<double> memb2(n2);

    // Copy the membership to vector<double>
    memcpy(memb1.data(),mxGetPr(p1),n1*sizeof(double));
    memcpy(memb2.data(),mxGetPr(p2),n2*sizeof(double));

    outputArgs[0] = mxCreateDoubleMatrix(1,1,mxREAL);
    outputArgs[1] = mxCreateDoubleMatrix(1,1,mxREAL);

    // Cast the membership to vector<int>
    std::vector<int> m1(memb1.begin(),memb1.end());
    std::vector<int> m2(memb2.begin(),memb2.end());

    // Compute number of communities of m1 and m2
    int q1 = std::set<int>(m1.begin(),m1.end()).size();
    int q2 = std::set<int>(m2.begin(),m2.end()).size();

    double the_nmi=0;

    bool agtb=true; //qa is greater than qb

    if(q1<q2)
        agtb=false;

    if(agtb)
        the_nmi=compute_nmi(m1,m2);
    else
        the_nmi=compute_nmi(m1,m2);

    // Correct the NMI by NMI of random partitions
    double tot_nmi=0;
    int nsample=100;
    for (int sample=0; sample<nsample; sample++)
    {
        double nmi=0;
        if (agtb)
        {
            shuffle_seq(m2,rg);
            nmi=compute_nmi(m1,m2);
        }
        else
        {
            shuffle_seq(m1,rg);
            nmi=compute_nmi(m2,m1);
        }
        tot_nmi += nmi;
    }
    tot_nmi /= nsample;

    // Compute rnmi
    *mxGetPr(outputArgs[0]) = the_nmi-tot_nmi;
    // Compute nmi
    *mxGetPr(outputArgs[1]) = the_nmi;

    return;
}
