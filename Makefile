.PHONY: all clean run

CXX = g++
CFLAGS = -O2 -Wall -std=c++11 -Wno-psabi
LDFLAGS = -lm
DBGFLAGS = -g -D_DEBUG_ON_
TRAIN_ITER = 100

SRC_DIR = src/
INC_DIR = inc/
INC_DIR = src/

INC = -I$(INC_DIR)
SRC = $(wildcard $(SRC_DIR)*.cpp)
OBJ = $(SRC: %.cpp = %.o)

TARGET: train #test

all: $(TARGET)
	@echo -n ""
train: src/train.o src/FB.o
	$(CXX) src/train.o src/FB.o $(CFLAGS) $(INC) -o $@
%.o : %.cpp inc/hmm.h
	${CXX} $< ${CFLAGS} $(INC) -lpthread -c
# src/train.o : src/train.cpp inc/hmm.h
# 	$(CXX) src/train.cpp $(CFLAGS) $(INC) -c
# src/FB.o : src/FB.h src/FB.cpp inc/hmm.h
# 	$(CXX) src/FB.cpp $(CFLAGS) $(INC) -c

clean:
	@rm -rf *.o src/*.o