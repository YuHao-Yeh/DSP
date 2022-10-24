#include "bits/stdc++.h"
#include "../inc/hmm.h"
#include "FB.h"
#include "tm_usage.h"

using namespace std;

void Help_message()
{
   puts("Usage: ./train <iter> <model_init_path> <seq_path> <output_model_path>");
}

void PrintHMM_(HMM hmm)
{
   cout << "This is main file." << endl;
   cout << "filename: " << hmm.model_name << endl;
   cout << "initial: " << hmm.state_num << endl;
   for (int i = 0; i < hmm.state_num; i++)
      cout << hmm.initial[i] << " ";
   puts("");
   cout << "transition: " << hmm.state_num << endl;
   for (int i = 0; i < hmm.state_num; i++)
   {
      for (int j = 0; j < hmm.state_num; j++)
         cout << hmm.transition[i][j] << " ";
      puts("");
   }
   puts("");
   cout << "observation: " << hmm.state_num << endl;
   for (int i = 0; i < hmm.observ_num; i++)
   {
      for (int j = 0; j < hmm.state_num; j++)
         cout << hmm.observation[i][j] << " ";
      puts("");
   }
   puts("");
   puts("----------------------------");
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
   PrintHMM_(hmm);
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

      // if (iter == iteration)
      //    ofs.open(fout, ios::out);
      // else
      //    ofs.open(model_init, ios::out);

      // vector<float> NewI = forbackward.GetNewI();
      // vector<vector<float>> NewT = forbackward.GetNewT();
      // vector<vector<float>> NewO = forbackward.GetNewO();

      // for (int i = 0; i < hmm.state_num; i++)
      //    hmm.initial[i] = NewI.at(i);
      // for (int i = 0; i < hmm.state_num; i++)
      //    for (int j = 0; j < hmm.state_num; j++)
      //       hmm.transition[i][j] = NewT.at(i)[j];
      // for (int i = 0; i < hmm.observ_num; i++)
      //    for (int j = 0; j < hmm.state_num; j++)
      //       hmm.observation[i][j] = NewO.at(i)[j];
      // dumpHMM(fp, &hmm);

      // Write File
      // if (!ofs.is_open())
      // {
      //    perror(model_init.c_str());
      //    exit(1);
      // }

      // ofs << "initial: " << hmm.state_num << endl;
      // for (int i = 0; i < hmm.state_num; i++)
      //    ofs << NewI.at(i) << " ";

      // ofs << "\n\ntransition: " << hmm.state_num << endl;
      // for (int i = 0; i < hmm.state_num; i++)
      // {
      //    for (int j = 0; j < hmm.state_num; j++)
      //       ofs << NewT.at(i)[j] << " ";
      //    ofs << endl;
      // }

      // ofs << "\n\nobservation: " << hmm.observ_num << endl;
      // for (int i = 0; i < hmm.state_num; i++)
      // {
      //    for (int j = 0; j < hmm.observ_num; j++)
      //       ofs << NewO.at(i)[j] << " ";
      //    ofs << endl;
      // }

      // ofs.close();
      tmusg2.getPeriodUsage(stat2);
      cout << (stat2.uTime + stat2.sTime) / 1000.0 << "ms" << endl;
      if (!((iter + 1) % 10))
      {
         cout << "Iteration: " << iter + 1 << " complete.\n";
         tmusg.getPeriodUsage(stat);
         cout << "The total CPU time: " << (stat.uTime + stat.sTime) / 1000.0 << "ms" << endl;
         cout << "memory: " << stat.vmPeak << "KB" << endl; // print peak memory
         puts("---------------------------------");
      }

      iter++;
   }

   forbackward.WriteHMM(fout);

   return 0;
}
