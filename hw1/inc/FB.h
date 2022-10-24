//----------------------------------------------------------------------
// File:       FB.h
// Author:     Yu-Hao Yeh
// Synopsis:   Realization of Forward-Backward Algorithm (Baum-Welch Algorithm)
// Date:       2022/10/25
//----------------------------------------------------------------------
#ifndef _FB_H_
#define _FB_H_

#include "bits/stdc++.h"
#include "../inc/hmm.h"

#ifndef defLINE
#define dLINE 10000
#endif

#ifndef defTIME
#define dTIME 50
#endif

#ifndef defDIGIT
#define dDIGIT 100000
#endif

using namespace std;

typedef struct ForwardBackward
{
   vector<vector<int>> seq;
   vector<vector<vector<double>>> a;         // alph[line][time][state]
   vector<vector<vector<double>>> b;         // beta[line][time][state]
   vector<vector<vector<double>>> g;         // gamma[line][time][state]
   vector<vector<vector<vector<double>>>> x; // xi[line][time][state][state]

   vector<double> newI;         // new Initial probability set[state]
   vector<vector<double>> newT; // new Transition matrix[state][state]
   vector<vector<double>> newO; // new Observation matrix[observation][state]

} FB;

class FBAlg
{
private:
   int iteration;
   int state;
   int obnum; // observation number
   FB fb;
   HMM hmm; // observation[observation][state] : obeservation i comes from state j
            // transition[state i][state j] : state i to state j

   // Calculate Variables; called by CalVar()
   void CalAlph();
   void CalBeta();
   void CalGamma();
   void CalXi();

   // Update data set; called by Update()
   void UpdateInitial();
   void UpdateTransitionA();
   void UpdateObservationB();
   void UpdateHMM();

public:
   FBAlg(HMM hmm_a, int iter) //: hmm(hmm), iteration(iteration)
   {
      state = hmm_a.state_num;
      obnum = hmm_a.observ_num;
      hmm = hmm_a;
      iteration = iter;
      fb.a.assign(dLINE, vector<vector<double>>(dTIME, vector<double>(state, 0.0)));
      fb.b.assign(dLINE, vector<vector<double>>(dTIME, vector<double>(state, 0.0)));
      fb.g.assign(dLINE, vector<vector<double>>(dTIME, vector<double>(state, 0.0)));
      fb.x.assign(dLINE, vector<vector<vector<double>>>(dTIME, vector<vector<double>>(state, vector<double>(state, 0.0))));
      fb.seq.assign(dLINE, vector<int>(dTIME, 0));
      fb.newI.assign(state, 0.0);
      fb.newT.assign(state, vector<double>(state, 0.0));
      fb.newO.assign(obnum, vector<double>(state, 0.0));
   }
   ~FBAlg(){};

   // Input sequence
   void GetSeq(string);

   // Calculate Variables
   void CalVar();

   // Update data set
   void Update();

   // Print
   void PrintHMM();

   // Write File
   void WriteHMM(string);
};

#endif