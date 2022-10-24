//----------------------------------------------------------------------
// File:       Viterbi.h
// Author:     Yu-Hao Yeh
// Synopsis:   Realization of Viterbi Algorithm
// Date:       2022/10/25
//----------------------------------------------------------------------
#ifndef _VITERBI_H_
#define _VITERBI_H_

#include "bits/stdc++.h"
#include "../inc/hmm.h"

#ifndef defMAX_NUM
#define dMAX_NUM 5
#endif

#ifndef defLINE
#define dLINE 2500
#endif

#ifndef defTIME
#define dTIME 50
#endif

using namespace std;

typedef struct ViterbiAlg
{
   vector<vector<vector<double>>> d; // delta[line][time][state]
   vector<vector<vector<short>>> p;  // psi[line][time][state]
} Vit;

class Viterbi
{
private:
   short state;
   short obnum;
   int modnum;
   vector<vector<int>> seq;
   vector<pair<int, double>> max; // // result probabilit and path (model, probability)
   HMM hmm[dMAX_NUM];
   Vit vit[dMAX_NUM];

   // Viterbi main process
   void CalDelta(int);
   void FindMax();

public:
   Viterbi(HMM (&hmm_a)[dMAX_NUM], int modnum_a)
   {
      modnum = modnum_a;
      state = hmm_a[0].state_num;
      obnum = hmm_a[0].state_num;
      for (int i = 0; i < modnum; i++)
      {
         hmm[i] = hmm_a[i];
         vit[i].d.assign(dLINE, vector<vector<double>>(dTIME, vector<double>(state, 0.0)));
         vit[i].p.assign(dLINE, vector<vector<short>>(dTIME, vector<short>(state, 0)));
      }
      seq.assign(dLINE, vector<int>(dTIME, 0));
      max.assign(dLINE, make_pair(0, 0.0));
   }
   ~Viterbi(){};

   // Test sequence
   void GetSeq(string);

   // Viterbi main process
   void Process();

   // Write File
   void WriteViterbi(string);
   void WriteAccuracy(); // Compared with open case at "./data/test_lbl.txt"
};

#endif