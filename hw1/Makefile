# .PHONY: all clean run
# CXX = g++
# CC = gcc
# CFLAGS = -Wall -std=c++11 -Wno-psabi
# OPTFLAGS = -g
# OPTFLAGS = -O3
# LDFLAGS = -lm
# DBGFLAGS = -g -D_DEBUG_ON_
# TRAIN_ITER = 100

# SRC_DIR = src/
# INC_DIR = inc/
# INC_DIR = src/

# TARGET = training # test
# INC = -I$(INC_DIR)
# SRC = $(wildcard $(SRC_DIR*.cpp))
# OBJ = $(SRC:%.cpp = %.o)

# all: $(TARGET)
# 	@echo -n ""

# train : src/train_opt.o src/FB_opt.o inc
# 	$(CXX) $(OPTFLAGS) train_opt.o FB_opt.o -hmm -Llib -o train
# train_opt.o : src/train.cpp inc/hmm.h
# 	$(CC) $(CFLAGS) $< -Ilib -o $@
# FB_opt.o : src/FB.cpp src/FB.h
# 	$(CC) $(CFLAGS) $(OPTFLAGS) $< -o $@
# #  $(CXX) $(CFLAGS) $(INC) $(OPTFLAGS) -o $@ $^  $(LDFLAGS) -Iinc
# # @echo "train FINISH "

# # test: src/test.cpp
# # 	$(CXX) -o $@ $^ $(CFLAGS) $(LDFLAGS) -Iinc
# # 	@echo "test FINISH"

# clean:
# 	rm -f $(TARGET)


CXX = g++
CFLAGS = -Wall -std=c++11 -Wno-psabi
OPTFLAGS = -g
OPTFLAGS = -O3

SRC_DIR = src/
INC_DIR = src/
INC_DIR = inc/

EXEC = train
BIN	= $(EXEC)
INC = -I$(INC_DIR)
SRC = $(wildcard $(SRC_DIR)*.cpp)
OBJ = $(SRC:%.cpp=%.o)

.PHONY: all clean

all: $(BIN)

$(BIN): $(OBJ)
	$(CXX) $(CFLAGS) $(INC) $(OPTFLAGS) $^ -o $(BIN)
	@echo "FINISH"

%.o : %.cpp
	@echo ">> compiling: $<"
	$(CXX) $(CFLAGS) -c $< -o $@
clean:
	rm -f $(BIN) $(OBJ)
	rm -rf *.o lib/*.a lib/*.o bin/*
 