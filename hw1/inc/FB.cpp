//----------------------------------------------------------------------
// File:       FB.cpp
// Author:     Yu-Hao Yeh
// Synopsis:   Realization of Forward-Backward Algorithm
// Date:       2022/10/25
//----------------------------------------------------------------------
#include "FB.h"

void FBAlg::GetSeq(string filepath)
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
         fb.seq.at(i)[j] = int(tmp[j] - 'A');
      i++;
   }

   ifs.close();
}

void FBAlg::CalVar()
{
   CalAlph();
   CalBeta();
   CalGamma();
   CalXi();
}

void FBAlg::CalAlph()
{
   // Base case : aplh_1(i) = pi_i * beta_i(obsevation_1)
   for (int l = 0; l < dLINE; l++)
      for (int j = 0; j < state; j++)
         fb.a.at(l)[0][j] = hmm.initial[j] * hmm.observation[fb.seq.at(l)[0]][j];

   //----------------------------------------------------------------
   // Inductive step:
   // alph_{t+1}(i) = sigma_{i=1}^N {alph_t(i) * tran_a_ij} * ob_b_j(o_{t+1})
   // 1 <= t <= T-1
   // 1 <= j <= N
   //----------------------------------------------------------------
   for (int l = 0; l < dLINE; l++)
      for (int t = 1; t < dTIME; t++)
      {
         double *tmp;
         tmp = (double *)calloc(state, sizeof(double));

         for (int j = 0; j < state; j++)
         {
            for (int i = 0; i < state; i++)
               tmp[j] += fb.a.at(l)[t - 1][i] * hmm.transition[i][j];
            fb.a.at(l)[t][j] = tmp[j] * hmm.observation[fb.seq.at(l)[t]][j];
         }
         free(tmp);
      }
}

void FBAlg::CalBeta()
{
   // Base case : beta_T(i) = 1
   for (int l = 0; l < dLINE; l++)
      for (int i = 0; i < state; i++)
         fb.b.at(l)[dTIME - 1][i] = 1.0;

   //----------------------------------------------------------------
   // Inductive step:
   // beta_t(i) = sigma_{j=1}^N tran_a_ij * ob_b_j(o_{t+1}) * beta_{t+1}(j)
   // t = T-1, T-2, ..., 1
   // 1 <= i <= N
   //----------------------------------------------------------------
   for (int l = 0; l < dLINE; l++)
      for (int t = dTIME - 2; t >= 0; t--)
         for (int i = 0; i < state; i++)
            for (int j = 0; j < state; j++)
               fb.b.at(l)[t][i] += hmm.transition[i][j] * hmm.observation[fb.seq.at(l)[t + 1]][j] * fb.b.at(l)[t + 1][j];
}

void FBAlg::CalGamma()
{
   for (int l = 0; l < dLINE; l++)
      for (int t = 0; t < dTIME; t++)
      {
         double sum = 0;
         double *tmp;
         tmp = (double *)calloc(state, sizeof(double));

         for (int j = 0; j < state; j++)
         {
            tmp[j] = fb.a.at(l)[t][j] * fb.b.at(l)[t][j];
            sum += tmp[j];
         }

         for (int i = 0; i < state; i++)
            fb.g.at(l)[t][i] = tmp[i] / sum;
      }
}

void FBAlg::CalXi()
{
   for (int l = 0; l < dLINE; l++)
      for (int t = 0; t < dTIME - 1; t++)
      {
         double sum = 0, tmp;

         for (int j = 0; j < state; j++)
         {
            for (int i = 0; i < state; i++)
            {
               tmp = fb.a.at(l)[t][i];
               tmp *= hmm.transition[i][j];
               tmp *= hmm.observation[fb.seq.at(l)[t + 1]][j];
               tmp *= fb.b.at(l)[t + 1][j];
               sum += tmp;
            }
         }

         for (int i = 0; i < state; i++)
            for (int j = 0; j < state; j++)
            {
               tmp = fb.a.at(l)[t][i];
               tmp *= hmm.transition[i][j];
               tmp *= hmm.observation[fb.seq.at(l)[t + 1]][j];
               tmp *= fb.b.at(l)[t + 1][j];
               fb.x.at(l)[t][i][j] = tmp / sum;
            }
      }
}

void FBAlg::Update()
{
   UpdateInitial();
   UpdateTransitionA();
   UpdateObservationB();
   UpdateHMM();
}

void FBAlg::UpdateInitial()
{
   double sum;
   for (int i = 0; i < state; i++)
   {
      sum = 0;
      for (int l = 0; l < dLINE; l++)
         sum += fb.g.at(l)[0][i];
      fb.newI.at(i) = sum / dLINE;
   }

   // Check Initial Matrix : Sum of each probability should equal to 1
   sum = 0.0;
   for (int i = 0; i < state; i++)
      sum += fb.newI.at(i);
   if (sum != 1.0)
      for (int i = 0; i < state; i++)
         fb.newI.at(i) /= sum;
   // fb.newI.at(i) = round(dDIGIT * (fb.newI.at(i) / sum)) / dDIGIT;
}

void FBAlg::UpdateTransitionA()
{
   double sum, sum_g, sum_x;

   for (int i = 0; i < state; i++)
      for (int j = 0; j < state; j++)
      {
         sum_g = 0.0;
         sum_x = 0.0;

         for (int l = 0; l < dLINE; l++)
            for (int t = 0; t < dTIME - 1; t++)
            {
               sum_x += fb.x.at(l)[t][i][j];
               sum_g += fb.g.at(l)[t][i];
            }

         fb.newT.at(i)[j] = sum_x / sum_g;
      }
   // Check Trasition Matrix : Sum of each row should equal to 1
   for (int i = 0; i < state; i++)
   {
      sum = 0.0;
      for (int j = 0; j < state; j++)
         sum += fb.newT.at(i)[j];
      if (sum != 1.0)
         for (int j = 0; j < state; j++)
            fb.newT.at(i)[j] /= sum;
      // fb.newT.at(i)[j] = round(dDIGIT * (fb.newT.at(i)[j] / sum)) / dDIGIT;
   }
}

void FBAlg::UpdateObservationB()
{
   double sum, sum_o;

   for (int j = 0; j < state; j++)
   {
      for (int k = 0; k < obnum; k++)
      {
         sum = 0.0;
         sum_o = 0.0;
         for (int l = 0; l < dLINE; l++)
            for (int t = 0; t < dTIME; t++)
            {
               if (k == fb.seq.at(l)[t])
                  sum_o += fb.g[l][t][j];
               sum += fb.g[l][t][j];
            }
         fb.newO.at(k)[j] = sum_o / sum;
      }
   }
   // Check Observation Matrix : Sum of each column should equal to 1
   for (int j = 0; j < obnum; j++)
   {
      sum = 0.0;
      for (int i = 0; i < state; i++)
         sum += fb.newT.at(i)[j];
      if (sum != 1.0)
         for (int i = 0; i < state; i++)
            fb.newO.at(i)[j] /= sum;
      // fb.newO.at(i)[j] = round(dDIGIT * (fb.newO.at(i)[j] / sum)) / dDIGIT;
   }
}

void FBAlg::UpdateHMM()
{
   for (int i = 0; i < state; i++)
      hmm.initial[i] = fb.newI.at(i);
   for (int i = 0; i < state; i++)
      for (int j = 0; j < state; j++)
         hmm.transition[i][j] = fb.newT.at(i)[j];
   for (int i = 0; i < obnum; i++)
      for (int j = 0; j < state; j++)
         hmm.observation[i][j] = fb.newO.at(i)[j];
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

void FBAlg::WriteHMM(string filepath)
{
   puts("----------------------------");
   cout << "Output filepath: " << filepath << endl;

   for (int i = 0; i < state; i++)
      fb.newI.at(i) = round(dDIGIT * fb.newI.at(i)) / dDIGIT;

   for (int i = 0; i < state; i++)
      for (int j = 0; j < state; j++)
         fb.newT.at(i)[j] = round(dDIGIT * fb.newT.at(i)[j]) / dDIGIT;

   for (int j = 0; j < obnum; j++)
      for (int i = 0; i < state; i++)
         fb.newO.at(i)[j] = round(dDIGIT * fb.newO.at(i)[j]) / dDIGIT;

   FILE *fp;
   fp = fopen(filepath.c_str(), "w");
   if (open_or_die(filepath.c_str(), "w"))
      dumpHMM(fp, &hmm);
   PrintHMM();
}