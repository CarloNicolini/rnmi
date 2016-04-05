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

#ifndef __NMI__
#define __NMI__
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <vector>
#include <string>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <ctime>
#include <map>
#include <algorithm>
#include <assert.h>
#include <unistd.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <sstream>
#include "zrg.h"
#define FRANDOM (rg.rdflt())
const char deli[1024]="\t ";

//{{{long get_cpu_time(void)
long get_cpu_time(void)
{
    struct rusage ru;
    getrusage(RUSAGE_SELF, &ru);
    return ( ru.ru_utime.tv_sec*1000 +
	     ru.ru_utime.tv_usec/1000+
	     ru.ru_stime.tv_sec*1000 +
	     ru.ru_stime.tv_usec/1000 );
}
//}}}

int ps2p(vector<string> pas, vector<int> & pa)
{
	// convert string vector to int, return number of groups in the partition
	map <string,int> id2idx;
	pa.resize(pas.size());
	int idx=0;
	for(int i=0;i<pas.size();i++){
		string id=pas[i];
		if(id2idx.count(id) == 0) id2idx[id] = idx++;
		pa[i] = id2idx[id];
	}
	return idx;
}

//{{{shuffle_seq(vector <int> &sequence)
void shuffle_seq(vector <int> &sequence,ZRANDOMv3 &rg)
{
	//random shuffle the sequence
	int N=sequence.size();
	//for(int i=0;i<N;i++) sequence[i]=i;
	int tmp;
	for(int i=0;i<N;i++){
		int tmpindex=int(FRANDOM*N);
		tmp=sequence[tmpindex];
		sequence[tmpindex]=sequence[N-i-1];
		sequence[N-i-1]=tmp;
	}
}
//}}}


double compute_nmi(vector<int> pa, vector<int> pb)
{
	assert(pa.size() == pb.size() && "Two partitions have different number of nodes !");
	int n=pa.size();
	int qa=-1,qb=-1;
	vector <int > ga;//group a
	vector <int > gb;//group b
	for(int i=0;i<n;i++){
		if(qa<pa[i]) qa=pa[i];
		if(qb<pb[i]) qb=pb[i];
	}
	qa++;
	qb++;
	if(qa==1 && qb==1) return 0.0;
	ga.resize(qa);
	for(int q=0;q<qa;q++) ga[q]=0;
	gb.resize(qb);
	for(int q=0;q<qb;q++) gb[q]=0;

	vector< vector<int> > A;
	vector< vector<int> > B;
	A.resize(qa); //existing structure
	B.resize(qa); //counting structure
	for(int i=0;i<n;i++){
		int q=pa[i];
		int t=pb[i];
		ga[q]++;
		gb[t]++;
		int idx=-1;
		for(int j=0;j<A[q].size();j++){
			if(A[q][j] == t) {
				idx=j;
				break;
			}
		}
		if(idx == -1){//pair [x y] did not show up
			A[q].push_back(t);
			B[q].push_back(1);
		}else{// [x y] is there
			B[q][idx] += 1;
		}
	}
	double Ha=0;
	for(int q=0;q<qa;q++){
		if(ga[q]==0) continue;
		double prob=1.0*ga[q]/n;
		Ha += prob*log(prob);
	}
	double Hb=0;
	for(int q=0;q<qb;q++){
		if(gb[q]==0) continue;
		double prob=1.0*gb[q]/n;
		Hb += prob*log(prob);
	}
	double Iab=0;
	for(int q=0;q<qa;q++){
		for(int idx=0;idx<A[q].size();idx++){
			double prob=1.0*B[q][idx]/n;
			int t=A[q][idx];
			Iab += prob*log(prob/ ( 1.0*ga[q]/n*gb[t]/n ));
		}
	}
	return -2.0*Iab/(Ha+Hb);
}


double compute_rnmi(vector<int> pa, vector<int> pb)
{
	assert(pa.size() == pb.size() && "Two partitions have different number of nodes !");
	int n=pa.size();
	int qa=-1,qb=-1;
	vector <int > ga;//group a
	vector <int > gb;//group b
	for(int i=0;i<n;i++){
		if(qa<pa[i]) qa=pa[i];
		if(qb<pb[i]) qb=pb[i];
	}
	qa++;
	qb++;
	ga.resize(qa);
	for(int q=0;q<qa;q++) ga[q]=0;
	gb.resize(qb);
	for(int q=0;q<qb;q++) gb[q]=0;

	vector< vector<int> > A;
	vector< vector<int> > B;
	A.resize(qa); //existing structure
	B.resize(qa); //counting structure
	for(int i=0;i<n;i++){
		int q=pa[i];
		int t=pb[i];
		ga[q]++;
		gb[t]++;
		int idx=-1;
		for(int j=0;j<A[q].size();j++){
			if(A[q][j] == t) {
				idx=j;
				break;
			}
		}
		if(idx == -1){//pair [x y] did not show up
			A[q].push_back(t);
			B[q].push_back(1);
		}else{// [x y] is there
			B[q][idx] += 1;
		}
	}
	double Ha=0;
	for(int q=0;q<qa;q++){
		if(ga[q]==0) continue;
		double prob=1.0*ga[q]/n;
		Ha += prob*log(prob);
	}
	Ha *= -1;
	double Hb=0;
	for(int q=0;q<qb;q++){
		if(gb[q]==0) continue;
		double prob=1.0*gb[q]/n;
		Hb += prob*log(prob);
	}
	Hb *= -1;
	double Iab=0;
	for(int q=0;q<qa;q++){
		for(int idx=0;idx<A[q].size();idx++){
			double prob=1.0*B[q][idx]/n;
            //int t=A[q][idx];
			Iab += prob*log(prob);
		}
	}
	Iab *=-1;
	Iab = Ha+Hb-Iab;
	double corr = (qa*qb-qa-qb+1.0)/2/n;
	return 2.0*(Iab-corr)/(Ha+Hb);
}











#endif
