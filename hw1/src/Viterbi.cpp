//----------------------------------------------------------------------
// File:       Viterbi.cpp
// Author:     Yu-Hao Yeh
// Synopsis:   Realization of Viterbi Algorithm (Baum-Welch Algorithm)
// Date:       2022/10/25
//----------------------------------------------------------------------
#include "../inc/Viterbi.h"

Viterbi::Viterbi(int modnum_a)
{
   modnum = modnum_a;
   line = 0;
   time = 0;
   state_num.assign(modnum, 0);
   observe_num.assign(modnum, 0);
   seq.assign(dLINE + 1, vector<int>(dMAX_SEQ, 0));
   seq_size.assign(dLINE, 0);
}

Viterbi::~Viterbi()
{
   state_num.clear();
   observe_num.clear();
   seq.clear();
   seq_size.clear();
   max.clear();
}

void Viterbi::RecvSeq(string filepath)
{
   ifstream ifs(filepath, ios::in);
   if (!ifs.is_open())
   {
      perror(filepath.c_str());
      exit(1);
   }

   char tmp[dMAX_SEQ + 1];
   while (ifs.getline((char *)&tmp, dMAX_SEQ + 1))
   {
      for (int j = 0; j < ifs.gcount(); j++)
         seq.at(line)[j] = int(tmp[j] - 'A');
      seq_size.at(line) = ifs.gcount() - 1;
      line++;
   }
   time = *max_element(seq_size.begin(), seq_size.end());
   seq.resize(line, vector<int>(time));
   seq_size.resize(line);

   ifs.close();
}

void Viterbi::RecvHMM(vector<string> modlist)
{
   ifstream ifs;
   for (int i = 0; i < modnum; i++)
   {
      loadHMM(&hmm[i], modlist.at(i).c_str());
      this->state_num.at(i) = hmm[i].state_num;
      this->observe_num.at(i) = hmm[i].state_num;
   }
   this->state = *max_element(state_num.begin(), state_num.end());
   this->obnum = hmm[0].state_num;
}

void Viterbi::StartVit()
{
   for (int i = 0; i < modnum; i++)
   {
      vit[i].d.assign(line, vector<vector<double>>(time, vector<double>(state, 0.0)));
      vit[i].p.assign(line, vector<vector<short>>(time, vector<short>(state, 0)));
      CalDelta(i);
   }
   max.assign(line, make_pair(0, 0.0));
   FindMax();
}

void Viterbi::CalDelta(int mod)
{
   for (int l = 0; l < line; l++)
      for (int i = 0; i < state; i++)
      {
         vit[mod].d.at(l)[0][i] = hmm[mod].initial[i] * hmm[mod].transition[seq.at(l)[0]][0];
         // vit[mod].p.at(l)[0][i] = 0; // It has been assigned 0 at initialization
      }

   for (int l = 0; l < line; l++)
      for (int t = 1; t < time; t++)
      {
         vector<double> tmp(time, 0.0);
         for (int j = 0; j < state; j++)
         {
            for (int i = 0; i < state; i++)
               tmp.at(i) = vit[mod].d.at(l)[t - 1][i] * hmm[mod].transition[i][j];

            double max = *max_element(tmp.begin(), tmp.end());
            int index = max_element(tmp.begin(), tmp.end()) - tmp.begin();

            vit[mod].d.at(l)[t][j] = max * hmm[mod].observation[seq.at(l)[t]][j];
            vit[mod].p.at(l)[t][j] = index;
         }
      }
}

void Viterbi::FindMax()
{
   vector<double> tmp(modnum, 0.0);
   for (int l = 0; l < line; l++)
   {
      for (int i = 0; i < modnum; i++)
         tmp.at(i) = vit[i].d.at(l)[time - 1][seq.at(l)[time - 1]];
      max[l].second = *max_element(tmp.begin(), tmp.end());
      max[l].first = max_element(tmp.begin(), tmp.end()) - tmp.begin() + 1;
   }
}

void Viterbi::WriteViterbi(string filepath)
{
   ofstream ofs(filepath, ios::out);
   for (int l = 0; l < line; l++)
      ofs << "model_0" << max.at(l).first << ".txt " << max.at(l).second << "\n";
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

   vector<int> ans(line, 0);
   char tmp[13];
   int i = 0;
   while (ifs.getline(tmp, 13))
   {
      ans.at(i) = int(tmp[7] - '0');
      i++;
   }
   int count = 0;
   for (int l = 0; l < line; l++)
      if (ans.at(l) == this->max.at(l).first)
         count++;

   ofs << count << " cases are correct.\n";
   ofs << "Accuracy = " << 1.0 * count / line * 100.0 << "%.\n";

   cout << count << " cases are correct.\n";
   cout << "Accuracy = " << 1.0 * count / line * 100.0 << "%.\n";

   ifs.close();
   ofs.close();
}