//----------------------------------------------------------------------
// File:       Viterbi.cpp
// Author:     Yu-Hao Yeh
// Synopsis:   Realization of Viterbi Algorithm (Baum-Welch Algorithm)
// Date:       2022/10/25
//----------------------------------------------------------------------
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
      CalDelta(i);

   FindMax();
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
               tmp.at(i) = vit[mod].d.at(l)[t - 1][i] * hmm[mod].transition[i][j];

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
      for (int i = 0; i < modnum; i++)
         tmp.at(i) = vit[i].d.at(l)[dTIME - 1][seq.at(l)[dTIME - 1]];
      max[l].second = *max_element(tmp.begin(), tmp.end());
      max[l].first = max_element(tmp.begin(), tmp.end()) - tmp.begin() + 1;
   }
}

void Viterbi::WriteAccuracy()
{
   string fin = "./data/test_lbl.txt";
   string fout = "./data/test_accuracy.txt";
   ifstream ifs(fin, ios::in);
   if (!ifs.is_open())
   {
      perror(fin.c_str());
      exit(1);
   }
   ofstream ofs(fout, ios::out);
   if (!ofs.is_open())
   {
      perror(fout.c_str());
      exit(1);
   }

   vector<int> ans(dLINE, 0);
   char tmp[13];
   int i = 0;
   while (ifs.getline(tmp, 13))
   {
      ans.at(i) = int(tmp[7] - '0');
      i++;
   }
   int count = 0;
   for (int l = 0; l < dLINE; l++)
      if (ans.at(l) == this->max.at(l).first)
         count++;

   ofs << count << " cases are corrected.\n";
   ofs << "Accuracy = " << 1.0 * count / dLINE * 100.0 << "%.\n";

   ifs.close();
   ofs.close();
}

void Viterbi::WriteViterbi(string filepath)
{
   ofstream ofs(filepath, ios::out);
   for (int l = 0; l < dLINE; l++)
      ofs << "model_0" << max.at(l).first << ".txt " << max.at(l).second << "\n";
}