#include "Viterbi.h"

void Viterbi::GetSeq(string filepath)
{
   ifstream ifs(filepath, ios::in);
   if (!ifs.is_open())
   {
      perror(filepath.c_str());
      exit(1);
   }

   int i = 0;
   char tmp[dTIME + 1];
   while (ifs.getline((char *)&tmp, dTIME + 1))
   {
      for (int j = 0; j < dTIME; j++)
         seq.at(i)[j] = int(tmp[j] - 'A');
      i++;
   }
   ifs.close();
}

void Viterbi::Process()
{
   for (int i = 0; i < modnum; i++)
   {
      // cout << "Process " << i << endl;
      CalDelta(i);
   }
   // puts("Finish process");
   FindMax();
   // puts("Max finded");
}

void Viterbi::CalDelta(int mod)
{
   for (int l = 0; l < dLINE; l++)
      for (int i = 0; i < state; i++)
      {
         vit[mod].d.at(l)[0][i] = hmm[mod].initial[i] * hmm[mod].transition[seq.at(l)[0]][0];
         vit[mod].p.at(l)[0][i] = 0;
      }

   for (int l = 0; l < dLINE; l++)
      for (int t = 1; t < dTIME; t++)
         for (int j = 0; j < state; j++)
         {
            vector<double> tmp(dTIME, 0.0);
            for (int i = 0; i < state; i++)
               tmp.at(i) = vit[mod].d.at(l)[t - 1][i];

            double max = *max_element(tmp.begin(), tmp.end());
            int index = max_element(tmp.begin(), tmp.end()) - tmp.begin();

            vit[mod].d.at(l)[t][j] = max * hmm[mod].observation[seq.at(l)[t]][j];
            vit[mod].p.at(l)[t][j] = index;
         }
}

void Viterbi::FindMax()
{
   vector<double> tmp(modnum, 0.0);
   for (int l = 0; l < dLINE; l++)
   {
      // cout << "current : " << l << endl;
      // puts("F1");
      for (int i = 0; i < modnum; i++)
         tmp.at(i) = vit[i].d.at(l)[dTIME - 1][seq.at(l)[dTIME - 1]];
      // puts("F2");
      max[l].second = *max_element(tmp.begin(), tmp.end());
      // puts("F3");
      max[l].first = max_element(tmp.begin(), tmp.end()) - tmp.begin();
      // tmp.clear();
   }
}

void Viterbi::PrintHMM(int h)
{
   cout << "filename: " << hmm[h].model_name << endl;
   cout << "initial: " << hmm[h].state_num << endl;
   for (int i = 0; i < hmm[h].state_num; i++)
      cout << hmm[h].initial[i] << " ";
   puts("");

   cout << "transition: " << hmm[h].state_num << endl;
   for (int i = 0; i < hmm[h].state_num; i++)
   {
      for (int j = 0; j < hmm[h].state_num; j++)
         cout << hmm[h].transition[i][j] << " ";
      puts("");
   }
   puts("");

   cout << "observation: " << hmm[h].state_num << endl;
   for (int i = 0; i < hmm[h].observ_num; i++)
   {
      for (int j = 0; j < hmm[h].state_num; j++)
         cout << hmm[h].observation[i][j] << " ";
      puts("");
   }
   puts("");
   puts("----------------------------");
}

void Viterbi::PrintAccuracy()
{
   string file = "./data/test_lbl.txt";
   ifstream ifs(file, ios::in);
   if (!ifs.is_open())
   {
      perror(file.c_str());
      exit(1);
   }
   vector<int> ans(dLINE, 0);
   char tmp[20];
   int i = 0;
   while (ifs.getline(tmp, 20))
   {
      ans.at(i) = int(tmp[7] - '0');
   }
   int count = 0;
   for (i = 0; i < dLINE; i++)
      if (ans.at(i) == max.at(i).first)
         count++;
   cout << "Accuracy = " << 1.0 * count / dLINE * 100.0 << "%.\n";
}

void Viterbi::WriteViterbi(string filepath)
{
   ofstream ofs(filepath, ios::out);
   for (int l = 0; l < dLINE; l++)
      ofs << "model_0" << max.at(l).first << ".txt " << max.at(l).second << "\n";
}