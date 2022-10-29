//----------------------------------------------------------------------
// File:       FB.h
// Author:     Yu-Hao Yeh
// Synopsis:   Implementation of Forward-Backward Algorithm
//             (Baum-Welch Algorithm)
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

   ~ForwardBackward()
   {
      a.clear();
      b.clear();
      g.clear();
      x.clear();
   }
} FB;

class FBAlg
{
private:
   vector<vector<int>> seq;
   int line;  // total lines of training sequences
   int time;  // longest length of training sequences
   int stnum; // total state number
   int obnum; // total observation number
   FB fb;     // forward-backward parameters
   HMM hmm;   // observation[observation][state] : obeservation i comes from state j
              // transition[state i][state j] : state i to state j

public:
   FBAlg(HMM);
   ~FBAlg();

   // Read training sequences
   void ReadSeq(string);

   // Start calculating forward-backward parameters
   void StartForBack();

   // Update HMM
   void UpdateHMM();

   // Print
   void PrintP();   // Print P(O | lambda) according to Alph_T
   void PrintHMM(); // Print info of HMM

   // Write File
   void WriteHMM(string);
};

#endif