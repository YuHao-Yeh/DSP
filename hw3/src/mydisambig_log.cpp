#include "../inc/mydisambig_log.h"

Log::Log(string log_file, int ngram_order = 2)
{
   this->flog = log_file;
   this->order = ngram_order;
}

Log::Log() {}

Log::~Log()
{
   this->ofs.clear();
}

void Log::Help_message(string e = "")
{
   if (!e.empty())
      cerr << "Failed to open " << e << " .\n";
   puts("Usage: ./mydisambig <segmented_file_path> <ZhuYin-Big5_mapping_path> <language_model_path> <output_file_path>.");
   exit(1);
}

bool Log::OpenLog(string log_file)
{
   this->flog = log_file;
   ofs.open(flog, ios::out | ios::app | ios::binary);
   if (!ofs.is_open())
   {
      ofs.close();
      return false;
   }
   ofs.close();
   return true;
}

void Log::WriteTMusage(string fin, double cpu_time, double peak_memory)
{
   ofs.open(flog, ios::out | ios::app | ios::binary);
   if (!ofs.is_open())
      Help_message(flog);

   time(&this->rawtime);
   this->info = localtime(&this->rawtime);
   strftime(buffer, 80, "%Y-%m-%d %a %H:%M:%S", info);

   ofs << "[" << buffer << "]\t";
   switch (this->order)
   {
   case 1:
      ofs << "Unigram\t";
      break;
   case 2:
      ofs << "Bbigram\t";
      break;
   case 3:
      ofs << "Triigram\t";
      break;
   default:
      break;
   }
   ofs << fin << "\t"
       << "CPU time: " << cpu_time << "s \t"
       << "Peak memory: " << peak_memory << "MB" << endl;

   ofs.close();
}

void Log::WriteError(string e)
{
   ofs.open(flog, ios::out | ios::app | ios::binary);
   if (!ofs.is_open())
      Help_message(flog);

   time(&this->rawtime);
   this->info = localtime(&this->rawtime);
   strftime(buffer, 80, "%Y-%m-%d %a %H:%M:%S", info);

   ofs << "[" << buffer << "]\t" << e << endl;
}