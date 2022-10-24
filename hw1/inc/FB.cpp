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
   // puts("AL");
   CalAlph();
   // puts("BE");
   CalBeta();
   // puts("GA");
   CalGamma();
   // puts("XI");
   CalXi();
}

void FBAlg::CalAlph()
{
   // Base case : aplh_1(i) = pi_i * beta_i(obsevation_1)
   // puts("A0");
   for (int l = 0; l < dLINE; l++)
      for (int j = 0; j < state; j++)
         fb.a.at(l)[0][j] = hmm.initial[j] * hmm.observation[fb.seq.at(l)[0]][j];
   // cout << "a[0][0][0] should be " << hmm.initial[0] * hmm.observation[fb.seq->at(0)[0]][0] << endl;
   // for (int l = 0; l < 5; l++)
   // {
   //    for (int j = 0; j < state; j++)
   //       cout << fb.a.at(l)[0][j] << " ";
   //    puts("");
   // }
   //----------------------------------------------------------------
   // Inductive step:
   // alph_{t+1}(i) = sigma_{i=1}^N {alph_t(i) * tran_a_ij} * ob_b_j(o_{t+1})
   // 1 <= t <= T-1
   // 1 <= j <= N
   //----------------------------------------------------------------
   // puts("A1");
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
   // puts("A2");
}

void FBAlg::CalBeta()
{
   // Base case : beta_T(i) = 1
   // puts("B1");
   for (int l = 0; l < dLINE; l++)
      for (int i = 0; i < state; i++)
         fb.b.at(l)[dTIME - 1][i] = 1.0;

   //----------------------------------------------------------------
   // Inductive step:
   // beta_t(i) = sigma_{j=1}^N tran_a_ij * ob_b_j(o_{t+1}) * beta_{t+1}(j)
   // t = T-1, T-2, ..., 1
   // 1 <= i <= N
   //----------------------------------------------------------------
   // puts("B2");

   for (int l = 0; l < dLINE; l++)
      for (int t = dTIME - 2; t >= 0; t--)
         for (int i = 0; i < state; i++)
            for (int j = 0; j < state; j++)
               fb.b.at(l)[t][i] += hmm.transition[i][j] * hmm.observation[fb.seq.at(l)[t + 1]][j] * fb.b.at(l)[t + 1][j];

   // puts("B3");
}

void FBAlg::CalGamma()
{
   // puts("G1");
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
   // puts("G2");
}

void FBAlg::CalXi()
{
   // puts("X1");
   for (int l = 0; l < dLINE; l++)
   {
      // puts("X2");
      for (int t = 0; t < dTIME - 1; t++)
      {
         double sum = 0, tmp;

         // puts("X3");
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

         // puts("X4");
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
   // puts("X5");
}

void FBAlg::Update()
{
   // puts("Update1");
   for (int i = 0; i < state; i++)
      fb.newI.at(i) = 0.0;
   // puts("Update2");
   for (int i = 0; i < state; i++)
      for (int j = 0; j < state; j++)
         fb.newT.at(i)[j] = 0.0;
   // puts("Update3");
   for (int i = 0; i < obnum; i++)
      for (int j = 0; j < state; j++)
         fb.newT.at(i)[j] = 0.0;
   // puts("Update4");
   UpdateInitial();
   // puts("Update5");
   UpdateTransitionA();
   // puts("Update6");
   UpdateObservationB();
   UpdateHMM();
   // PrintHMM();
   Print();
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
}

void FBAlg::UpdateTransitionA()
{
   double sum_g, sum_x;

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
}

void FBAlg::UpdateHMM()
{
   int digit = 100000;
   for (int i = 0; i < state; i++)
      hmm.initial[i] = round(digit * fb.newI.at(i)) / digit;
   for (int i = 0; i < state; i++)
      for (int j = 0; j < state; j++)
         hmm.transition[i][j] = round(digit * fb.newT.at(i)[j]) / digit;
   for (int i = 0; i < obnum; i++)
      for (int j = 0; j < state; j++)
         hmm.observation[i][j] = round(digit * fb.newO.at(i)[j]) / digit;
}

void FBAlg::Print()
{
   // for (int i = 0; i < dTIME; i++)
   //    cout << fb.seq->at(0)[i] << " ";
   // puts("");

   // puts(".....................");

   // for (int i = 0; i < state; i++)
   // {
   //    for (int j = 0; j < dTIME; j++)
   //       cout << fb.a.at(0)[j][i] << " ";
   //    puts("");
   // }

   // puts("--------------------");

   // for (int i = 0; i < 5; i++)
   // {
   //    for (int j = 0; j < dTIME; j++)
   //       cout << fb.b.at(i)[j][0] << " ";
   //    puts("");
   // }

   // puts(".....................");

   // for (int i = 0; i < 5; i++)
   // {
   //    for (int j = 0; j < dTIME; j++)
   //       cout << fb.g.at(i)[j][0] << " ";
   //    puts("");
   // }

   // puts(".....................");
   double sum = 0.0;
   for (int i = 0; i < state; i++)
      sum += fb.a.at(dLINE - 1)[dTIME - 1][i];
   cout << "P(O | lambda) = " << sum << endl;
}

void FBAlg::PrintHMM()
{
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

void FBAlg::WriteHMM(string filepath)
{
   cout << "Output filepath: " << filepath << endl;
   FILE *fp;
   fp = fopen(filepath.c_str(), "w");
   if (open_or_die(filepath.c_str(), "w"))
      dumpHMM(fp, &hmm);
   PrintHMM();
}