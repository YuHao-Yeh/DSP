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

#ifndef defDIGIT
#define dDIGIT 100000
#endif

using namespace std;

typedef struct ForwardBackward
{
   vector<vector<vector<double>>> a;         // alph[line][time][state]
   vector<vector<vector<double>>> b;         // beta[line][time][state]
   vector<vector<vector<double>>> g;         // gamma[line][time][state]
   vector<vector<vector<vector<double>>>> x; // xi[line][time][state][state]

   vector<double> newI;         // new Initial probability set[state]
   vector<vector<double>> newT; // new Transition matrix[state][state]
   vector<vector<double>> newO; // new Observation matrix[observation][state]

   ~ForwardBackward()
   {
      a.clear();
      b.clear();
      g.clear();
      x.clear();
      newI.clear();
      newT.clear();
      newO.clear();
   }
} FB;

class FBAlg
{
private:
   vector<vector<int>> seq;
   vector<short> seq_size;
   int line;
   int time;
   int stnum;
   int obnum; // observation number
   FB fb;
   HMM hmm; // observation[observation][state] : obeservation i comes from state j
            // transition[state i][state j] : state i to state j

   // Re-construct fb after read sequences
   void ConstructFB();

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
   FBAlg(HMM hmm_a);
   ~FBAlg();

   // Input sequence
   void ReadSeq(string);

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