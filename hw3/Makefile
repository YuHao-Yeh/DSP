SRIPATH ?= /root/srilm-1.5.10
MACHINE_TYPE ?= i686-m64
INC_PATH ?= inc
SRC_PATH ?= src

# Language model
CORPUS ?= corpus.txt
CORPUS_SEP ?= corpus_sep.txt
CORPUS_NGRAM ?= corpus_ngram.txt
ORDER ?= 2
LM ?= bigram.lm
LM3 ?= trigram.lm

CXX = g++
CXXFLAGS = -O2 -I$(SRIPATH)/include -I$(INC_PATH) # -Wall
OPTFLAGS = -g
LDFLAGS = -lpthread
SRLIBFLAGS = -loolm -ldstruct -lmisc
vpath lib%.a $(SRIPATH)/lib/$(MACHINE_TYPE)
vpath %.c $(SRC_PATH)
vpath %.cpp $(SRC_PATH)

ifeq ($(TARGET), $(MYDISAMBIG))
	# CXXFLAGS := $(filter-out -D _MYDISAMBIG_TRIGRAM_, $(CXXFLAGS))
   CXXFLAGS += -D _MYDISAMBIG_
endif

ifeq ($(TARGET), $(MYDISAMBIG_TRIGRAM))
	# CXXFLAGS := $(filter-out -D _MYDISAMBIG_, $(CXXFLAGS))
   CXXFLAGS += -D _MYDISAMBIG_TRIGRAM_
endif

MYDISAMBIG = mydisambig
MYDISAMBIG_TRIGRAM = mydisambig_trigram
TARGET: $(MYDISAMBIG) $(MYDISAMBIG_TRIGRAM)
INC = I$(INC_PATH)
SRC = $(wildcard $(SRC_PATH)/*.cpp)
OBJ = $(SRC: .cpp = .o)
FROM ?= Big5-ZhuYin.map
TO ?= ZhuYin-Big5.map

.PHONY: all mclean map 

all: $(TARGET)

# ./mydisambig test_data/example.txt ZhuYin-Big5.map separator_big5.pl test_data/result1.txt
# $(TARGET): $(OBJ) -loolm -ldstruct -lmisc
# 	$(CXX) $(LDFLAGS) -o $@ $^

$(MYDISAMBIG): tm_usage.o mydisambig_log.o mydisambig.o $(SRLIBFLAGS)
	@echo "Building mydisambig"
	$(CXX) $(CXXFLAGS) $(OPTFLAGS) $(LDFLAGS) $^ -o mydisambig -D _$(shell echo $(MYDISAMBIG) | tr '[a-z]' '[A-Z]')_
	@echo "Finished"

$(MYDISAMBIG_TRIGRAM): tm_usage.o mydisambig_log.o mydisambig_trigram.o $(SRLIBFLAGS)
	@echo "Building mydisambig_trigram"
	$(CXX) $(CXXFLAGS) $(OPTFLAGS) $(LDFLAGS) $^ -o mydisambig_trigram -D _$(shell echo $(MYDISAMBIG_TRIGRAM) | tr '[a-z]' '[A-Z]')_
	@echo "Finished"

%.o : %.cpp
	@echo ">> compiling: $<"
	$(CXX) $(CXXFLAGS) -lpthread -c $< -o $@

run_bi:
	mkdir -p result_bi
	@echo "Start mydisambig"
	@for i in $(shell seq 1 10) ; do \
		./mydisambig test_data/$$i.txt $(TO) $(LM) result_bi/$$i.txt | iconv -f big5; \
		echo "Finish $$i.txt"; \
		echo "----------------------------"; \
	done;
	@rm -rf *.o

run_tri:
	mkdir -p result_tri
	@echo "Start mydisambig_trigram"
	@for i in $(shell seq 1 10) ; do \
		./mydisambig_trigram test_data/$$i.txt $(TO) $(LM3) result_tri/$$i.txt | iconv -f big5; \
		echo "Finish $$i.txt"; \
		echo "----------------------------"; \
	done;
	@rm -rf *.o

# make map Big5-ZhuYin.map ZhuYin-Big5.map
map: $(SRC_PATH)/mapping.cpp
	@echo "Start mapping"
	$(CXX) $(CXXFLAGS) -o $@ $^
	./map $(FROM) $(TO) | iconv -f big5
	@#python3 $(SRC_PATH)/mapping.py $(FROM) $(TO)
	@echo "Finish mapping"

clean:
	$(RM) $(TARGET)
	@rm -rf map mydisambig mydisambig_trigram
	@rm -rf *.o src/*.o inc/*.o int/*.gch

test_bi:
	mkdir -p result_bi
	./mydisambig test_data/example.txt $(TO) $(LM) result_bi/example_result.txt
	@echo "Finish testing example"

test_tri:
	mkdir -p result_tri
	./mydisambig_trigram test_data/example.txt $(TO) $(LM3) result_tri/example_result.txt
	@echo "Finish testing example"

ans:
	mkdir -p  answer_2
	mkdir -p answer_3
	@for i in $(shell seq 1 10) ; do \
		disambig -text test_data/$$i.txt -map $(TO) -lm $(LM) -order 2 -keep-unk > answer_2/$$i.txt | iconv -f big5; \
		disambig -text test_data/$$i.txt -map $(TO) -lm $(LM) -order 3 -keep-unk > answer_3/$$i.txt | iconv -f big5; \
		echo "Finish $$i.txt"; \
	done;
	disambig -text test_data/example.txt -map $(TO) -lm $(LM) -order 2 -keep-unk > answer_2/example.txt | iconv -f big5;
	disambig -text test_data/example.txt -map $(TO) -lm $(LM) -order 3 -keep-unk > answer_3/example.txt | iconv -f big5;
	@echo "Get all answers"	

lm:
	bash ./setup.sh;
	@echo "Start building bigram language model"
	perl separator_big5.pl  $(CORPUS) > $(CORPUS_SEP);
	ngram-count -text $(CORPUS_SEP) -write $(CORPUS_NGRAM) -order $(ORDER)
	ngram-count -read $(CORPUS_NGRAM) -lm $(LM) -order $(ORDER) -unk | iconv -f big5
	@echo "Finished"
	
	@echo "Start building trigram language model"
	perl separator_big5.pl  $(CORPUS) > $(CORPUS_SEP);
	ngram-count -text $(CORPUS_SEP) -write $(CORPUS_NGRAM) -order 3 | iconv -f big5
	ngram-count -read $(CORPUS_NGRAM) -lm $(LM3) -order 3 -unk
	rm $(CORPUS_SEP) $(CORPUS_NGRAM)
	@echo "Finished"
