#!/usr/bin/make

all_proxy: all

VERSION=0.0.1

# Toolchain programs.
CC=g++
LD=ld
AR=ar
LN=ln


# Directories.
SRC=src
INC=src
DATA=data

# Source and header files.
SRC_OBJS=main Groups online trillionDTW util deque TimeSeries OnlineSession Grouping Algo
INC_OBJS=     Groups online trillionDTW util deque TimeSeries OnlineSession Grouping Algo

OBJ_PATHS=$(addprefix $(SRC)/, $(addsuffix .o, $(SRC_OBJS)))
SRC_PATHS=$(addprefix $(SRC)/, $(addsuffix .cpp, $(SRC_OBJS)))
INC_PATHS=$(addprefix $(INC)/, $(addsuffix .h, $(INC_OBJS)))



# Compiler flags.
CCFLAGS=-Isrc -std=c++0x -Wall -Wextra -Wno-unused-parameter -Wformat -Wpedantic -O2
DEBUG=#-pg -g

# Default rule for compiling object files.
%.o: %.cpp $(INC_PATHS)
	$(CC) $(CCFLAGS) $(DEBUG) -c $< -o $@

$(SRC)/util.o:        $(INC)/util.h
$(SRC)/trillionDTW.o: $(INC)/util.h $(INC)/deque.h $(INC)/trillionDTW.h
$(SRC)/online.o:      $(INC)/util.h $(INC)/deque.h $(INC)/online.h
$(SRC)/main.o:        $(INC)/trillionDTW.h $(INC)/Groups.h $(INC)/online.h
$(SRC)/Groups.o:      $(INC)/Groups.h
$(SRC)/deque.o:       $(INC)/deque.h


# Rule for compiling main executable.
onex: $(OBJ_PATHS)
	$(CC) $(DEBUG) $^ -o $@

run: onex
	./onex

# Main rule.
all: onex

# Clean repo.
clean:
	rm -f $(OBJ_PATHS) onex

.PHONY: all_proxy all clean
