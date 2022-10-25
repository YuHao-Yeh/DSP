This is README file for DSP HW#1
Author: 葉育豪
Date: 2022/10/25

========================
SYNOPSIS:

train :
	./train <iteration_number> <model_init_path> <seq_path> <output_model_path>

test :
	./test <models_list_path> <seq_path> <output_result_path>

========================
DIRECTORY:

bin/	   executable binary
data/    input data (train sequences and test sequences)
doc/	   reports
inc/	   include library (*.h files)
outputs/ test result
src/ 	   source C++ codes (*.cpp files)
model_01~5.txt
Makefile
README

========================
HOW TO CIMPILE:

To compile the demo, simply follow the following steps
	make

========================
HOW TO RUN:

	./train <iteration_number> <model_init_path> <seq_path> <output_model_path>

	./test <models_list_path> <seq_path> <output_result_path>

	For example,
	under hw1
	./bin/train 100 ./model_init.txt ./data/train_seq_01.txt ./model_01.txt
	./bin/train 100 ./model_init.txt ./data/seq_model_02.txt ./model_02.txt
	./bin/train 100 ./model_init.txt ./data/seq_model_03.txt ./model_03.txt
 	./bin/train 100 ./model_init.txt ./data/seq_model_04.txt ./model_04.txt
 	./bin/train 100 ./model_init.txt ./data/seq_model_05.txt ./model_05.txt
	
	./bin/test ./modellist.txt ./data/test_seq.txt ./outputs/result.txt
========================
OTHER NOTICE:
     hw1.zip/
     +- hw1/
      +- inc/
      |  +- hmm.h, FB.h, Viterbi.h
      +- src/
      |  +- train.cpp, FB.cpp      : train model by forward-backward algorithm.
      |  +- test.cpp, Viterbi.cpp   : test model by Viterbi algorithm
      +- model_init.txt                     : initial training model
      +- report.pdf
      +- Makefile
