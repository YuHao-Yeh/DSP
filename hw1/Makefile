.PHONY: all clean run

CXX = g++
CFLAGS = -O2 -Wall -std=c++11 -Wno-psabi
OPTFLAGS = -g
LDFLAGS = -lm
DBGFLAGS = -g -D_DEBUG_ON_
TRAIN_ITER = 100

SRC_DIR = src/
INC_DIR = inc/
INC_DIR = src/

INC = -I$(INC_DIR)
SRC = $(wildcard $(SRC_DIR)*.cpp)
OBJ = $(SRC: %.cpp = %.o)

TARGET: train test

all: $(TARGET)
	@echo -n ""

train: src/train.o src/FB.o src/tm_usage.o
	$(CXX) $(CFLAGS) $(INC) $(OPTFLAGS) $^ -o train
	@echo "FINISH train"

test: src/test.o src/Viterbi.o src/tm_usage.o
	$(CXX) $(CFLAGS) $(INC) $(OPTFLAGS) $^ -o test
	@echo "FINISH test"

%.o : %.cpp
	@echo ">> compiling: $<"
	$(CXX) $(CFLAGS) $(INC) -lpthread -c $< -o $@

clean:
	rm -f $(OBJ)
	rm -rf *.o src/*.o inc/*.o src/*.o