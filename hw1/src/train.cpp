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
   CommonNs::TmUsage tmusg2;
   CommonNs::TmStat stat2;

   //-------------------------------------------------------------------
   // Read Input File
   //-------------------------------------------------------------------
   int iteration = atoi(argv[1]); // iteration number
   string model_init = argv[2];   // model initial path
   string fseq = argv[3];         // sequence path
   string fout = argv[4];         // output model path

   HMM hmm;
   loadHMM(&hmm, model_init.c_str());
   FBAlg forbackward(hmm, iteration);
   forbackward.GetSeq(fseq);
   //-------------------------------------------------------------------
   // Main Part
   //-------------------------------------------------------------------
   short iter = 0;
   tmusg.periodStart();
   tmusg2.periodStart();

   while (iter < iteration)
   {
      // ofstream ofs;
      forbackward.CalVar();
      forbackward.Update();
      cout << "iter " << iter << " : ";

      tmusg2.getPeriodUsage(stat2);
      cout << (stat2.uTime + stat2.sTime) / 1000.0 / 1000.0 << "s" << endl;
      if (!((iter + 1) % 10))
      {
         cout << "Iteration: " << iter + 1 << " complete.\n";
         tmusg.getPeriodUsage(stat);
         cout << "The total CPU time: " << (stat.uTime + stat.sTime) / 1000.0 / 1000.0 << "s" << endl;
         cout << "memory: " << stat.vmPeak / 1024.0 << "MB" << endl; // print peak memory
         puts("---------------------------------");
      }

      iter++;
   }

   forbackward.WriteHMM(fout);

   return 0;
}
