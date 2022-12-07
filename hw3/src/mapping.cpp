#include "bits/stdc++.h"

using namespace std;

void Help_message()
{
   puts("Usage: ./map <Big5-ZhuYin_mapping_file_path> <ZhuYin-Big5_mapping_file_path>");
}

int main(int argc, char *argv[])
{
   if (argc != 3)
   {
      Help_message();
      exit(1);
   }

   string fin = argv[1];        // Big5-ZhuYin mapping file
   string fout = argv[2];       // ZhuYin-Big5 mapping file
   vector<string> phonetic[37]; // 37 ZhuYin Symbols
   vector<string> dictionary;   // Record every Mandarin dictionary

   //-------------------------------------------------------------------
   // Read Input File and Mapping
   //-------------------------------------------------------------------
   ifstream ifs;
   ifs.open(fin, ios::in | ios::binary);
   if (!ifs.is_open())
   {
      cerr << "Failed to open input file.\n";
      Help_message();
      exit(1);
   }

   string word, ZhuYin;
   while (ifs >> word >> ZhuYin)
   {
      dictionary.push_back(word);

      // First Zhu-Yin
      if (ZhuYin[0] == -93)
         if (ZhuYin[1] > -96 && ZhuYin[1] < -69)
            phonetic[ZhuYin[1] + 106].push_back(word);
         else if (ZhuYin[1] > 115 && ZhuYin[1] < 127)
            phonetic[ZhuYin[1] - 116].push_back(word);

      // Check Polyphones
      for (int i = 2; i < ZhuYin.length() - 2; i++)
      {
         if (ZhuYin[i] != 0x2f || ZhuYin[i + 1] != -93) // /
            continue;

         if (ZhuYin[i + 2] > -96 && ZhuYin[i + 2] < -69) // ㄐ-ㄩ
         {
            if (ZhuYin[i + 2] != ZhuYin[1])
               phonetic[ZhuYin[i + 2] + 106].push_back(word);
         }
         else if (ZhuYin[i + 2] > 0x73 && ZhuYin[i + 2] < 0x7f) // ㄅ-ㄏ
            if (ZhuYin[i + 2] != ZhuYin[1])
               phonetic[ZhuYin[i + 2] - 0x74].push_back(word);
      }
   }
   ifs.close();

   //-------------------------------------------------------------------
   // Write Output File
   //-------------------------------------------------------------------
   ofstream ofs;
   ofs.open(fout, ios::out | ios::binary);
   if (!ofs.is_open())
   {
      cerr << "Failed to open output file.\n";
      Help_message();
      exit(1);
   }

   // ㄅ-ㄦ
   for (int i = 0; i < 37; i++)
   {
      if (phonetic[i].empty())
         continue;

      if (i < 11)
         ofs << (char)0xa3 << (char)(0x74 + i) << "  ";
      else
         ofs << (char)0xa3 << (char)(0x96 + i) << "  ";

      for (int j = 0; j < phonetic[i].size(); j++)
         ofs << phonetic[i][j] << " ";
      ofs << endl;
   }

   // Every words
   for (int i = 0; i < dictionary.size(); i++)
      ofs << dictionary[i] << "  " << dictionary[i] << endl;

   ofs.close();

   return 0;
}

/*
   注音 Big5 mapping:
   ㄅ~ㄏ: A374~A37E
   ㄐ~ㄦ: A3A1~A3B7
   ㄧ~ㄩ: A3B8~A3BA
   ˉ~ˋ  : A3BC~A3BF
*/