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
   string flist = argv[1];
   string fseq = argv[2];
   string fout = argv[3];

   //-------------------------------------------------------------------
   // Start Viterbi Algorithm
   //-------------------------------------------------------------------
   HMM hmm[5];
   // puts("A");
   load_models(flist.c_str(), hmm, dMAX_NUM);
   // puts("B");
   Viterbi vit(hmm, dMAX_NUM);
   // puts("C");
   vit.GetSeq(fseq);
   // puts("D");
   vit.Process();

   //-------------------------------------------------------------------
   // Write File
   //-------------------------------------------------------------------
   // puts("E");
   vit.WriteViterbi(fout);

   // puts("F");
   vit.PrintAccuracy();

   // puts("G");
   tmusg.getPeriodUsage(stat);
   cout << "The total CPU time: " << (stat.uTime + stat.sTime) / 1000.0 / 1000.0 << "s" << endl;
   cout << "memory: " << stat.vmPeak / 1024.0 << "MB" << endl; // print peak memory

   // puts("H");
   return 0;
}