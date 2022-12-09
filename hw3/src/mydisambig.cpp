//----------------------------------------------------------------------
// File:       mydisambig.cpp
// Author:     Yu-Hao Yeh
// Synopsis:   Convert ZhuYin to Mandarin, and find the most possible
//             path via bigram HMM and Viterbi Algorithm
// Date:       2022/12/07
// Reference:  https://desilinguist.org/pdf/langmodel_code.pdf
//----------------------------------------------------------------------
// #ifdef _MYDISAMBIG_

#include "bits/stdc++.h"
// #include <iostream>
// #include <fstream>
// #include <string>
// #include <algorithm>
// #include <vector>
// #include <map>
// #include <set>
// #include <File.h>
#include "Ngram.h"
#include "../inc/tm_usage.h"
#include "../inc/mydisambig_log.h"

using namespace std;

const double dMIN = -1e8;
const int Ngram_Order = 2;
Vocab voc;
Ngram lm(voc, Ngram_Order);
map<string, set<string>> mapping; // record the mapping in ZhuYin-Big5
// string flog = "mydisambig_log.txt";

void Help_message(string s = "")
{
   if (!s.empty())
      cerr << "Fail to open" << s << endl;
   puts("Usage: ./mydisambig <segmented_file_path> <ZhuYin-Big5_mapping_path> <language_model_path> <output_file_path>.");
   exit(1);
}

// Get P(W2 | W1) -- bigram
double getBigramProb(const char *w1, const char *w2, Vocab &voc, Ngram &lm)
{
   VocabIndex wid1 = voc.getIndex(w1); // test whether the word is in the corpus
   VocabIndex wid2 = voc.getIndex(w2);

   if (wid1 == Vocab_None) // OOV
      wid1 = voc.getIndex(Vocab_Unknown);
   if (wid2 == Vocab_None) // OOV
      wid2 = voc.getIndex(Vocab_Unknown);

   VocabIndex context[] = {wid1, Vocab_None};
   return lm.wordProb(wid2, context);
}

int main(int argc, char *argv[])
{
   if (argc != 5)
      Help_message();

   CommonNs::TmUsage tmusg;
   CommonNs::TmStat stat;
   tmusg.periodStart();

   string fin = argv[1];  // segmented filepath
   string fmap = argv[2]; // ZhuYin-Big5 filepath
   string flm = argv[3];  // language model filepath
   string fout = argv[4]; // output filepath
   string line;

   //-------------------------------------------------------------------
   // Load language model to Ngram model
   //-------------------------------------------------------------------
   File lmFile(flm.c_str(), "r");
   if (!lmFile)
      Help_message("language model file");
   lm.read(lmFile);
   lmFile.close();

   //-------------------------------------------------------------------
   // Get mapping
   //-------------------------------------------------------------------
   ifstream ifs(fmap, ios::in | ios::binary);
   if (!ifs.is_open())
      Help_message("ZhuYin-Big5 mapping file");

   while (getline(ifs, line))
   {
      string word, symbol = line.substr(0, 2);
      set<string> &current = mapping[symbol];
      line = line.substr(3);
      line.erase(remove(line.begin(), line.end(), ' '), line.end());
      for (int i = 0; i < (line.size() - 1); i += 2)
         current.insert(line.substr(i, 2));
   }
   ifs.close();
   mapping[Vocab_SentStart].insert(Vocab_SentStart);
   mapping[Vocab_SentEnd].insert(Vocab_SentEnd);

   //-------------------------------------------------------------------
   // Get test data
   //-------------------------------------------------------------------
   vector<string> test; // record the lines in segmented file
   ifs.open(fin, ios::in | ios::binary);
   if (!ifs.is_open())
      Help_message("segmented file");

   while (getline(ifs, line))
      test.push_back(line);

   ifs.close();
   //-------------------------------------------------------------------
   // Viterbi algorithm : max P(q_i|q_j) * �__{t?1}(q_j)
   //-------------------------------------------------------------------
   ofstream ofs(fout, ios::out | ios::binary);
   if (!ofs.is_open())
      Help_message("output file");

   for (int l = 0; l < test.size(); l++)
   {
      line = test.at(l);
      vector<string> words;

      // Split the sentence into words
      line.erase(remove(line.begin(), line.end(), ' '), line.end());
      for (int i = 0; i < line.length() / 2; i++)
         words.push_back(line.substr(2 * i, 2));
      words.push_back(Vocab_SentEnd);

      int len = words.size();
      vector<vector<string>> str(len);
      vector<vector<int>> path(len);
      vector<vector<double>> prob(len);

      // Record all candidate route
      for (int t = 0; t < words.size(); t++)
         if (words.at(t)[0] != -93)
            str.at(t).push_back(words.at(t));
         else
            for (auto it : mapping[words.at(t)])
               str.at(t).push_back(it);

      // Initialization
      for (auto it : str.at(0))
         prob.at(0).push_back(getBigramProb(Vocab_SentStart, it.c_str(), voc, lm));
      path.at(0).assign(str.at(0).size(), 0);

      // Main part
      for (int t = 1; t < len; t++)
      {
         prob.at(t).assign(str.at(t).size(), dMIN);
         path.at(t).assign(str.at(t).size(), 0);
         for (int curr = 0; curr < str.at(t).size(); curr++)
         {
            if (str.at(t - 1).size() == 1)
            {
               prob.at(t)[curr] = prob.at(t - 1)[0] + getBigramProb(str.at(t - 1)[0].c_str(), str.at(t)[curr].c_str(), voc, lm);
               continue;
            }
            vector<double> tmp(str.at(t - 1).size(), 0.0);
            for (int prev = 0; prev < str.at(t - 1).size(); prev++)
               tmp.at(prev) = prob.at(t - 1)[prev] + getBigramProb(str.at(t - 1)[prev].c_str(), str.at(t)[curr].c_str(), voc, lm);

            double max = *max_element(tmp.begin(), tmp.end());
            int index = max_element(tmp.begin(), tmp.end()) - tmp.begin();

            prob.at(t)[curr] = max;
            path.at(t)[curr] = index;
         }
      }

      // Backtrack the most possible path
      vector<string> answer{Vocab_SentStart, Vocab_SentEnd};
      int index = path.at(str.size() - 1)[0];
      for (int t = str.size() - 1; t > 0; t--)
      {
         answer.insert(answer.begin() + 1, str.at(t - 1)[index]);
         index = path.at(t - 1)[index];
      }

      // Write file
      for (int t = 0; t < answer.size(); t++)
         ofs << answer.at(t) << " ";
      ofs << endl;
   }

   ofs.close();

   tmusg.getPeriodUsage(stat);
   cout << "The total CPU time: " << (stat.uTime + stat.sTime) / 1000.0 / 1000.0 << "s" << endl;
   cout << "memory: " << stat.vmPeak / 1024.0 << "MB" << endl; // print peak memory

   // logging
   Log log(Ngram_Order);
   log.WriteTMusage(fin, (stat.uTime + stat.sTime) / 1000.0 / 1000.0, stat.vmPeak / 1024.0);

   return 0;
}

// #endif