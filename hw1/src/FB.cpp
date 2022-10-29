//----------------------------------------------------------------------
// File:       FB.cpp
// Author:     Yu-Hao Yeh
// Synopsis:   Implementation of Forward-Backward Algorithm
// Date:       2022/10/25
//----------------------------------------------------------------------
#include "../inc/FB.h"

FBAlg::FBAlg(HMM hmm_a)
{
   hmm = hmm_a;
   line = 0;
   time = 0; // fixed sequence length in this hw
   stnum = hmm_a.state_num;
   obnum = hmm_a.observ_num;
   seq.assign(dLINE, vector<int>(dMAX_SEQ, 0));
}

FBAlg::~FBAlg()
{
   seq.clear();
};

void FBAlg::ReadSeq(string filepath)
{
   ifstream ifs(filepath, ios::in);
   if (!ifs.is_open())
   {
      perror(filepath.c_str());
      exit(1);
   }

   vector<short> seq_size(dLINE, 0);
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

   // Construct fb
   fb.a.assign(line, vector<vector<double>>(time, vector<double>(stnum, 0.0)));
   fb.b.assign(line, vector<vector<double>>(time, vector<double>(stnum, 0.0)));
   fb.g.assign(line, vector<vector<double>>(time, vector<double>(stnum, 0.0)));
   fb.x.assign(line, vector<vector<vector<double>>>(time, vector<vector<double>>(stnum, vector<double>(stnum, 0.0))));

   ifs.close();
}

void FBAlg::StartForBack()
{
   //----------------------------------------------------------------
   // Base case :
   // 1. aplh_1(i) = pi_i * beta_i(obsevation_1)
   // 2. beta_T(i) = 1
   //----------------------------------------------------------------
   int T = time - 1, tmpt;
   double sum;
   for (int l = 0; l < line; l++)
      for (int i = 0; i < stnum; i++)
      {
         fb.a.at(l)[0][i] = hmm.initial[i] * hmm.observation[seq.at(l)[0]][i];
         fb.b.at(l)[T][i] = 1.0;
      }

   //----------------------------------------------------------------
   // Parameter : Alph
   // Inductive step:
   // alph_{t+1}(i) = sigma_{i=1}^{i=N} {alph_t(i) * tran_a_ij} * ob_b_j(o_{t+1})
   // t = [1, T-1], j = [1, N]
   //----------------------------------------------------------------
   for (int l = 0; l < line; l++)
      for (int t = 0; t < T; t++)
      {
         tmpt = t + 1;
         for (int j = 0; j < stnum; j++)
         {
            for (int i = 0; i < stnum; i++)
               fb.a.at(l)[tmpt][j] += fb.a.at(l)[t][i] * hmm.transition[i][j];
            fb.a.at(l)[tmpt][j] *= hmm.observation[seq.at(l)[tmpt]][j];
         }
      }
   //----------------------------------------------------------------
   // Parameter : Beta
   // Inductive step:
   // beta_t(i) = sigma_{j=1}^N tran_a_ij * ob_b_j(o_{t+1}) * beta_{t+1}(j)
   // t = T-1, T-2, ..., 1
   // 1 <= i <= N
   //----------------------------------------------------------------
   for (int l = 0; l < line; l++)
      for (int t = time - 2; t >= 0; t--)
      {
         tmpt = t + 1;
         for (int i = 0; i < stnum; i++)
            for (int j = 0; j < stnum; j++)
               fb.b.at(l)[t][i] += hmm.transition[i][j] * hmm.observation[seq.at(l)[tmpt]][j] * fb.b.at(l)[tmpt][j];
      }
   //----------------------------------------------------------------
   // Parameter : Xi
   // Inductive step:
   // Xi_t(i, j) = P(q_t = i, q_{t+1} = j | O, lambda) = n/d
   // n = alph_t(i) * tran_a_ij * ob_b_j(0_{t+1}) * beta_{t+1}(j)
   // d = sigma_{i=1}^{i=N} {alph_t(i) * beta_{t}(i)}}
   // t = 1, 2, ... T-1
   //----------------------------------------------------------------
   // Parameter : Gamma
   // Inductive step:
   // Gamma_t(i) = P(q_t = i | O, lambda) = n/d
   //            = sigma_t{j=1}^{j=N} xi_t(i, j), t=[1, T-1]
   //              + alph_T(i) * beta_T(j) / sigma_{j=1}^{j=N} alph_T(i) * beta_{t+1}(j)
   // n = alph_t(i) * beta_t(j)
   // d = sigma_{j=1}^{j=N} alph_t(i) * beta_{t+1}(j)
   //----------------------------------------------------------------
   for (int l = 0; l < line; l++)
      for (int t = 0; t < time - 1; t++)
      {
         sum = 0;
         // sigma_{i=1}^{i=N} alph_t(i) * beta_t(i)
         for (int i = 0; i < stnum; i++)
            sum += fb.a.at(l)[t][i] * fb.b.at(l)[t][i];

         for (int i = 0; i < stnum; i++)
         {
            // Gamma
            fb.g.at(l)[t][i] = fb.a.at(l)[t][i] * fb.b.at(l)[t][i] / sum;
            // Xi
            for (int j = 0; j < stnum; j++)
               fb.x.at(l)[t][i][j] = fb.a.at(l)[t][i] * hmm.transition[i][j] * hmm.observation[seq.at(l)[t + 1]][j] * fb.b.at(l)[t + 1][j] / sum;
         }
      }
   // Gamma_T(i)
   for (int l = 0; l < line; l++)
   {
      sum = 0;
      for (int i = 0; i < stnum; i++)
         sum += fb.a.at(l)[T][i] * fb.b.at(l)[T][i];
      for (int i = 0; i < stnum; i++)
         fb.g.at(l)[T][i] = fb.a.at(l)[T][i] * fb.b.at(l)[T][i] / sum;
   }
}

void FBAlg::UpdateHMM()
{
   int T = time - 1;
   //----------------------------------------------------------------
   // Sum of gamma
   //----------------------------------------------------------------
   vector<double> sum_g(stnum, 0.0);
   for (int l = 0; l < line; l++)
      for (int i = 0; i < stnum; i++)
         for (int t = 0; t < time; t++)
            sum_g[i] += fb.g.at(l)[t][i];

   //----------------------------------------------------------------
   // Probability : Pi
   // Inductive step:
   //----------------------------------------------------------------
   double check_sum = 0;
   for (int i = 0; i < stnum; i++)
   {
      for (int l = 0; l < line; l++)
         hmm.initial[i] += fb.g.at(l)[0][i];
      hmm.initial[i] /= line;
      check_sum += hmm.initial[i];
   }

   // Check Initial Matrix : Sum of each probability should equal to 1
   if (check_sum != 1.0)
      for (int i = 0; i < stnum; i++)
         hmm.initial[i] /= check_sum;

   //----------------------------------------------------------------
   // Probability : Observation
   //----------------------------------------------------------------
   double sum;

   for (int i = 0; i < stnum; i++)
   {
      for (int k = 0; k < obnum; k++)
      {
         sum = 0;
         for (int l = 0; l < line; l++)
            for (int t = 0; t < time; t++)
               if (k == seq.at(l)[t])
                  sum += fb.g.at(l)[t][k];
         hmm.observation[k][i] = sum / sum_g[i];
      }
   }

   // Check Observation Matrix : Each column's sum = 1
   for (int j = 0; j < obnum; j++)
   {
      check_sum = 0.0;
      for (int i = 0; i < stnum; i++)
         check_sum += hmm.observation[i][j];
      if (check_sum != 1.0)
         for (int i = 0; i < obnum; i++)
            hmm.observation[i][j] /= check_sum;
   }

   //----------------------------------------------------------------
   // Probability : Transmission
   //----------------------------------------------------------------
   vector<double> tmp(stnum, 0);
   for (int l = 0; l < line; l++)
      for (int i = 0; i < stnum; i++)
         tmp[i] = sum_g[i] - fb.g.at(l)[T][i];
   for (int i = 0; i < stnum; i++)
      for (int j = 0; j < stnum; j++)
      {
         for (int l = 0; l < line; l++)
            for (int t = 0; t < T; t++)
               hmm.transition[i][j] += fb.x.at(l)[t][i][j];

         hmm.transition[i][j] /= tmp[i];
      }

   // Check Trasition Matrix : Each row's sum = 1
   for (int i = 0; i < stnum; i++)
   {
      check_sum = 0.0;
      for (int j = 0; j < stnum; j++)
         check_sum += hmm.transition[i][j];
      if (check_sum != 1.0)
         for (int j = 0; j < stnum; j++)
            hmm.transition[i][j] /= check_sum;
   }
}

void FBAlg::PrintHMM()
{
   cout << "Initial filename: " << hmm.model_name << endl
        << endl;
   cout << "initial: " << hmm.state_num << endl;
   for (int i = 0; i < hmm.state_num; i++)
      cout << setw(11) << hmm.initial[i] << " ";
   puts("");
   puts("");
   cout << "transition: " << hmm.state_num << endl;
   for (int i = 0; i < hmm.state_num; i++)
   {
      for (int j = 0; j < hmm.state_num; j++)
         cout << setw(11) << hmm.transition[i][j] << " ";
      puts("");
   }
   puts("");
   cout << "observation: " << hmm.state_num << endl;
   for (int i = 0; i < hmm.observ_num; i++)
   {
      for (int j = 0; j < hmm.state_num; j++)
         cout << setw(11) << hmm.observation[i][j] << " ";
      puts("");
   }
   puts("----------------------------");
}

void FBAlg::PrintP()
{
   double sum = 0, T = time - 1;
   for (int l = 0; l < line; l++)
      for (int i = 0; i < stnum; i++)
         sum += fb.a.at(l)[T][i];
   sum /= line;
   cout << "P(O | lambda) = " << sum << endl;
}

void FBAlg::WriteHMM(string filepath)
{
   puts("----------------------------");
   cout << "Output filepath: " << filepath << endl;

   for (int i = 0; i < stnum; i++)
      hmm.initial[i] = round(dDIGIT * hmm.initial[i]) / dDIGIT;

   for (int i = 0; i < stnum; i++)
   {
      for (int j = 0; j < stnum; j++)
         hmm.transition[i][j] = round(dDIGIT * hmm.transition[i][j]) / dDIGIT;
      for (int j = 0; j < obnum; j++)
         hmm.observation[i][j] = round(dDIGIT * hmm.observation[i][j]) / dDIGIT;
   }

   FILE *fp;
   fp = fopen(filepath.c_str(), "w");
   if (open_or_die(filepath.c_str(), "w"))
      dumpHMM(fp, &hmm);
   PrintHMM();
}