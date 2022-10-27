//----------------------------------------------------------------------
// File:       Viterbi.h
// Author:     Yu-Hao Yeh
// Synopsis:   Realization of Viterbi Algorithm
// Date:       2022/10/25
//----------------------------------------------------------------------
#ifndef _VITERBI_H_
#define _VITERBI_H_

#include "bits/stdc++.h"
#include <thread>
#include "../inc/hmm.h"

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

   ~ViterbiAlg()
   {
      d.clear();
      p.clear();
   }
} Vit;

class Viterbi
{
private:
   int line;
   int time;
   vector<short> state_num;
   vector<short> observe_num;
   short state;
   short obnum;
   int modnum;
   vector<vector<int>> seq;
   vector<short> seq_size;
   vector<pair<int, double>> max; // // result probability and path (model, probability)
   HMM hmm[dMAX_LINE];
   Vit vit[dMAX_LINE];

   // Viterbi main process
   void CalDelta(int);
   void FindMax();

public:
   Viterbi(int modnum_a);
   ~Viterbi();

   // Test sequence
   void RecvSeq(string);

   // Receive HMMs
   void RecvHMM(vector<string>);

   // Start Viterbi Algorithm
   void StartVit();

   // Write File
   void WriteViterbi(string);
   void WriteAccuracy(); // Compared with open case at "./data/test_lbl.txt"
};

#endif