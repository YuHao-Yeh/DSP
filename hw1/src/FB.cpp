//----------------------------------------------------------------------
// File:       FB.cpp
// Author:     Yu-Hao Yeh
// Synopsis:   Realization of Forward-Backward Algorithm
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
   seq_size.assign(dLINE, 0);
}

FBAlg::~FBAlg()
{
   seq.clear();
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

   ConstructFB();

   ifs.close();
}

void FBAlg::ConstructFB()
{
   fb.a.assign(line, vector<vector<double>>(time, vector<double>(stnum, 0.0)));
   fb.b.assign(line, vector<vector<double>>(time, vector<double>(stnum, 0.0)));
   fb.g.assign(line, vector<vector<double>>(time, vector<double>(stnum, 0.0)));
   fb.x.assign(line, vector<vector<vector<double>>>(time, vector<vector<double>>(stnum, vector<double>(stnum, 0.0))));
   fb.newI.assign(stnum, 0.0);
   fb.newT.assign(stnum, vector<double>(stnum, 0.0));
   fb.newO.assign(obnum, vector<double>(stnum, 0.0));
}

void FBAlg::CalVar()
{
   CalAlph();
   CalAlph();
   CalBeta();
   CalGamma();
   CalXi();
}

void FBAlg::CalAlph()
{
   // Base case : aplh_1(i) = pi_i * beta_i(obsevation_1)
   for (int l = 0; l < line; l++)
      for (int j = 0; j < stnum; j++)
         fb.a.at(l)[0][j] = hmm.initial[j] * hmm.observation[seq.at(l)[0]][j];

   //----------------------------------------------------------------
   // Inductive step:
   // alph_{t+1}(i) = sigma_{i=1}^N {alph_t(i) * tran_a_ij} * ob_b_j(o_{t+1})
   // 1 <= t <= T-1
   // 1 <= j <= N
   //----------------------------------------------------------------
   double *tmp;
   for (int l = 0; l < line; l++)
      for (int t = 1; t < time; t++)
      {

         tmp = (double *)calloc(stnum, sizeof(double));

         for (int j = 0; j < stnum; j++)
         {
            for (int i = 0; i < stnum; i++)
               tmp[j] += fb.a.at(l)[t - 1][i] * hmm.transition[i][j];
            fb.a.at(l)[t][j] = tmp[j] * hmm.observation[seq.at(l)[t]][j];
         }
         free(tmp);
      }
}

void FBAlg::CalBeta()
{
   // Base case : beta_T(i) = 1
   for (int l = 0; l < line; l++)
      for (int i = 0; i < stnum; i++)
         fb.b.at(l)[time - 1][i] = 1.0;

   //----------------------------------------------------------------
   // Inductive step:
   // beta_t(i) = sigma_{j=1}^N tran_a_ij * ob_b_j(o_{t+1}) * beta_{t+1}(j)
   // t = T-1, T-2, ..., 1
   // 1 <= i <= N
   //----------------------------------------------------------------
   for (int l = 0; l < line; l++)
      for (int t = time - 2; t >= 0; t--)
         for (int i = 0; i < stnum; i++)
            for (int j = 0; j < stnum; j++)
               fb.b.at(l)[t][i] += hmm.transition[i][j] * hmm.observation[seq.at(l)[t + 1]][j] * fb.b.at(l)[t + 1][j];
}

void FBAlg::CalGamma()
{
   double sum;
   double *tmp;
   for (int l = 0; l < line; l++)
      for (int t = 0; t < time; t++)
      {
         sum = 0;
         tmp = (double *)calloc(stnum, sizeof(double));

         for (int j = 0; j < stnum; j++)
         {
            tmp[j] = fb.a.at(l)[t][j] * fb.b.at(l)[t][j];
            sum += tmp[j];
         }

         for (int i = 0; i < stnum; i++)
            fb.g.at(l)[t][i] = tmp[i] / sum;
         free(tmp);
      }
}

void FBAlg::CalXi()
{
   double sum;
   vector<vector<double>> tmp(stnum, vector<double>(stnum));

   for (int l = 0; l < line; l++)
      for (int t = 0; t < time - 1; t++)
      {
         sum = 0;
         for (int j = 0; j < stnum; j++)
         {
            for (int i = 0; i < stnum; i++)
            {
               tmp.at(i)[j] = fb.a.at(l)[t][i];
               tmp.at(i)[j] *= hmm.transition[i][j];
               tmp.at(i)[j] *= hmm.observation[seq.at(l)[t + 1]][j];
               tmp.at(i)[j] *= fb.b.at(l)[t + 1][j];
               sum += tmp.at(i)[j];
            }
         }

         for (int i = 0; i < stnum; i++)
            for (int j = 0; j < stnum; j++)
               fb.x.at(l)[t][i][j] = tmp.at(i)[j] / sum;
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
   for (int i = 0; i < stnum; i++)
   {
      sum = 0;
      for (int l = 0; l < line; l++)
         sum += fb.g.at(l)[0][i];
      fb.newI.at(i) = sum / line;
   }

   // Check Initial Matrix : Sum of each probability should equal to 1
   sum = 0.0;
   for (int i = 0; i < stnum; i++)
      sum += fb.newI.at(i);
   if (sum != 1.0)
      for (int i = 0; i < stnum; i++)
         fb.newI.at(i) /= sum;
   // fb.newI.at(i) = round(dDIGIT * (fb.newI.at(i) / sum)) / dDIGIT;
}

void FBAlg::UpdateTransitionA()
{
   double sum, sum_g, sum_x;

   for (int i = 0; i < stnum; i++)
      for (int j = 0; j < stnum; j++)
      {
         sum_g = 0.0;
         sum_x = 0.0;

         for (int l = 0; l < line; l++)
            for (int t = 0; t < time - 1; t++)
            {
               sum_x += fb.x.at(l)[t][i][j];
               sum_g += fb.g.at(l)[t][i];
            }

         fb.newT.at(i)[j] = sum_x / sum_g;
      }
   // Check Trasition Matrix : Sum of each row should equal to 1
   for (int i = 0; i < stnum; i++)
   {
      sum = 0.0;
      for (int j = 0; j < stnum; j++)
         sum += fb.newT.at(i)[j];
      if (sum != 1.0)
         for (int j = 0; j < stnum; j++)
            fb.newT.at(i)[j] /= sum;
      // fb.newT.at(i)[j] = round(dDIGIT * (fb.newT.at(i)[j] / sum)) / dDIGIT;
   }
}

void FBAlg::UpdateObservationB()
{
   double sum, sum_o;

   for (int j = 0; j < stnum; j++)
   {
      for (int k = 0; k < obnum; k++)
      {
         sum = 0.0;
         sum_o = 0.0;
         for (int l = 0; l < line; l++)
            for (int t = 0; t < time; t++)
            {
               if (k == seq.at(l)[t])
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
      for (int i = 0; i < stnum; i++)
         sum += fb.newT.at(i)[j];
      if (sum != 1.0)
         for (int i = 0; i < stnum; i++)
            fb.newO.at(i)[j] /= sum;
      // fb.newO.at(i)[j] = round(dDIGIT * (fb.newO.at(i)[j] / sum)) / dDIGIT;
   }
}

void FBAlg::UpdateHMM()
{
   for (int i = 0; i < stnum; i++)
      hmm.initial[i] = fb.newI.at(i);
   for (int i = 0; i < stnum; i++)
      for (int j = 0; j < stnum; j++)
         hmm.transition[i][j] = fb.newT.at(i)[j];
   for (int i = 0; i < obnum; i++)
      for (int j = 0; j < stnum; j++)
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

   for (int i = 0; i < stnum; i++)
      fb.newI.at(i) = round(dDIGIT * fb.newI.at(i)) / dDIGIT;

   for (int i = 0; i < stnum; i++)
      for (int j = 0; j < stnum; j++)
         fb.newT.at(i)[j] = round(dDIGIT * fb.newT.at(i)[j]) / dDIGIT;

   for (int j = 0; j < obnum; j++)
      for (int i = 0; i < stnum; i++)
         fb.newO.at(i)[j] = round(dDIGIT * fb.newO.at(i)[j]) / dDIGIT;

   FILE *fp;
   fp = fopen(filepath.c_str(), "w");
   if (open_or_die(filepath.c_str(), "w"))
      dumpHMM(fp, &hmm);
   PrintHMM();
}