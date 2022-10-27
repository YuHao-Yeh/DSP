//----------------------------------------------------------------------
// File:       train.cpp
// Author:     Yu-Hao Yeh
// Synopsis:   HMM training via Forward-Backward Algorithm
// Date:       2022/10/25
//----------------------------------------------------------------------
#include "bits/stdc++.h"
#include "../inc/hmm.h"
#include "../inc/FB.h"
#include "../inc/tm_usage.h"

using namespace std;

void Help_message()
{
   puts("Usage: ./train <iter> <model_init_path> <seq_path> <output_model_path>");
}

int main(int argc, char *argv[])
{
   if (argc != 5)
   {
      Help_message();
      exit(1);
   }
   CommonNs::TmUsage tmusg;
   CommonNs::TmStat stat;

   //-------------------------------------------------------------------
   // Read Input File
   //-------------------------------------------------------------------
   int iteration = atoi(argv[1]); // iteration number
   string model_init = argv[2];   // model initial path
   string fseq = argv[3];         // sequence path
   string fout = argv[4];         // output model path

   HMM hmm;
   loadHMM(&hmm, model_init.c_str());
   FBAlg forbackward(hmm);
   forbackward.ReadSeq(fseq);
   //-------------------------------------------------------------------
   // Main Part
   //-------------------------------------------------------------------
   short iter = 0;
   tmusg.periodStart();

   while (iter < iteration)
   {
      forbackward.CalVar();
      forbackward.Update();

      iter++;
      if (!(iter % 20))
         cout << iter << " iterations are done.\n";
   }

   forbackward.WriteHMM(fout);

   tmusg.getPeriodUsage(stat);
   cout << "The total CPU time: " << (stat.uTime + stat.sTime) / 1000.0 / 1000.0 << "s" << endl;
   cout << "memory: " << stat.vmPeak / 1024.0 << "MB" << endl; // print peak memory
   puts("----------------------------");

   return 0;
}
