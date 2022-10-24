#include "bits/stdc++.h"
using namespace std;

typedef struct ForwardBackward
{
   vector<vector<int>> *seq;
   vector<vector<vector<double>>> a;         // alph[ttttttt][state]
   vector<vector<vector<double>>> b;         // beta[ttttttt][state]
   vector<vector<vector<double>>> g;         // gamma[ttttttt][state]
   vector<vector<vector<vector<double>>>> x; // xi[ttttttt][state][state]

   vector<double> newI;         // new Initial probability set[state]
   vector<vector<double>> newT; // new Transition matrix[state][state]
   vector<vector<double>> newO; // new Observation matrix[observation][state]

} FB;

int TIME = 50;
int LINE = 10000;

int main()
{
   FB fb;
   puts("Get in...");
   cout << "size : " << fb.seq->size() << endl;
   fb.seq = new vector<vector<int>>(LINE, vector<int>(TIME, 0));

   string filepath = "./train_seq_01.txt";
   puts("Opening file...");
   ifstream ifs(filepath, ios::in);

   if (!ifs.is_open())
   {
      puts("Failed");
      perror(filepath.c_str());
      exit(1);
   }
   puts("Success");
   puts("Reading file...");

   int i = 0;
   char tmp[TIME + 1];
   // char = (char *)calloc(TIME + 1, sizeof(char));
   while (ifs.getline(tmp, TIME + 1))
   {
      for (int j = 0; j < TIME; j++)
      {
         // cout << tmp[j] << " ";
         fb.seq->at(i)[j] = int(tmp[j] - 'A');
      }
      i++;
   }
   puts("");

   puts("Start");
   for (i = 0; i < TIME; i++)
      cout << fb.seq->at(0)[i] << " ";
   puts("");

   free(tmp);
   // tmp = nullptr;
   ifs.close();

   system("pause");
   return 0;
}
