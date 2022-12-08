// #include "bits/stdc++.h"
#include <iostream>
#include <fstream>
#include <string>
#include <chrono>
#include <ctime>

using namespace std;

class Log
{
private:
   string flog = "mydisambig_log.txt";
   int order = 2;
   ofstream ofs;
   time_t rawtime;
   struct tm *info;
   char buffer[80];

public:
   Log(string, int);
   Log(int);
   Log();
   ~Log();

   void Help_message(string);
   bool OpenLog(string);
   void WriteTMusage(string, double, double);
   void WriteError(string);
};