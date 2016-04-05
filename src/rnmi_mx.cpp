/* This file is part of rnmi_mx a program to compute normalized mutual
* information and relative normalized mutual information of two clusterings.
*
*  Copyright (C) 2016 Carlo Nicolini <carlo.nicolini@iit.it>
*
*  rnmi_mx is free software; you can redistribute it and/or
*  modify it under the terms of the GNU Lesser General Public
*  License as published by the Free Software Foundation; either
*  version 3 of the License, or (at your option) any later version.
*
*  Alternatively, you can redistribute it and/or
*  modify it under the terms of the GNU General Public License as
*  published by the Free Software Foundation; either version 2 of
*  the License, or (at your option) any later version.
*
*  rnmi_mx is distributed in the hope that it will be useful, but WITHOUT ANY
*  WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
*  FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License or the
*  GNU General Public License for more details.
*
*  You should have received a copy of the GNU Lesser General Public
*  License and a copy of the GNU General Public License along with
*  rnmi_mx. If not, see <http://www.gnu.org/licenses/>.
*/

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
    mexPrintf("\t\nComputation of NMI and relative NMI of two partitions (Zhang,2015).Membership vectors must be in the [0,N] range.\n");
    mexPrintf("Usage:\n");
    mexPrintf("[rnmi_value, nmi_value] = rnmi(membA, membB, nparts)\n");
    mexPrintf("Input:\n");
    mexPrintf("\tmembA: membership vector of partition A\n");
    mexPrintf("\tmembB: membership vector of partition B\n");
    mexPrintf("\tnparts: the number of random partitions to correct against (default 100)\n");
    mexPrintf("Output:\n");
    mexPrintf("\trnmi_value: value of relative Normalized Mutual Information with 100 samples (default)\n");
    mexPrintf("\tnmi_value: value of Normalized Mutual Information\n");
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
    if (nInputArgs<2)
    {
        printUsage();
        mexErrMsgTxt("Not enough arguments");
        return;
    }
    int randseed=int(time(NULL));
    ZRANDOMv3 rg(randseed);

    const mxArray* p1 = inputArgs[0];
    const mxArray* p2 = inputArgs[1];


    int nsample=100;
    if (nInputArgs==3)
        nsample = int(*mxGetPr(inputArgs[2]));

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

    // Check if NaN exist
    double v1 = 0;
    for (int i=0; i<n1; ++i)
        v1 += *mxGetPr(p1) + *mxGetPr(p2);
    if ( mxIsNaN(v1) || mxIsInf(v1))
    {
        mexErrMsgTxt("Partitions contain NaN or Inf values. Only integers are allowed");
        return;
    }

    // Cast the membership to vector<int>
    std::vector<int> m1(memb1.begin(),memb1.end());
    std::vector<int> m2(memb2.begin(),memb2.end());

    // Compute number of communities of m1 and m2m
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
