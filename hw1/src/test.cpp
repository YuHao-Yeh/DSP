//----------------------------------------------------------------------
// File:       test.cpp
// Author:     Yu-Hao Yeh
// Synopsis:   HMM finding most possible path via Viterbi Algorithm
// Date:       2022/10/25
//----------------------------------------------------------------------
#include "bits/stdc++.h"
#include "../inc/hmm.h"
#include "../inc/Viterbi.h"
#include "../inc/tm_usage.h"

using namespace std;

void Help_message()
{
   puts("Usage: ./test <models_list_path> <seq_path> <output_result_path>");
}

int main(int argc, char *argv[])
{
   if (argc != 4)
   {
      Help_message();
      exit(1);
   }
   CommonNs::TmUsage tmusg;
   CommonNs::TmStat stat;

   //-------------------------------------------------------------------
   // Read Input File
   //-------------------------------------------------------------------
   string flist = argv[1]; // model list filepath
   string fseq = argv[2];  // test sequence filepath
   string fout = argv[3];  // output filepath

   tmusg.periodStart();

   ifstream ifs(flist, ios::in);
   vector<string> modname;
   short modnum = 0;
   char tmp[dMAX_LINE + 1];
   while (ifs.getline((char *)&tmp, dMAX_LINE + 1))
   {
      modname.push_back(tmp);
      modnum++;
   }

   ifs.close();

   //-------------------------------------------------------------------
   // Start Viterbi Algorithm
   //-------------------------------------------------------------------
   Viterbi vit(modnum);
   vit.RecvSeq(fseq);
   vit.RecvHMM(modname);
   vit.StartVit();

   //-------------------------------------------------------------------
   // Write File
   //-------------------------------------------------------------------
   vit.WriteViterbi(fout);

   vit.WriteAccuracy();

   // dump_models(hmm, 5);

   tmusg.getPeriodUsage(stat);
   cout << "The total CPU time: " << (stat.uTime + stat.sTime) / 1000.0 / 1000.0 << "s" << endl;
   cout << "memory: " << stat.vmPeak / 1024.0 << "MB" << endl; // print peak memory

   return 0;
}