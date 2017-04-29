# Linux Makefile

GPP = g++
FLAGS = -g -Wall -D_REENTRANT -std=c++0x -pthread
INCLUDE = -isystem include -I /usr/include
LIBS = -lconfig++ -lpng

DIR_PATHS = bin obj

CPP_FILES = $(wildcard src/*.cpp)

SRC_FILES = $(CPP_FILES)
OBJ_FILES = $(addprefix obj/, $(notdir $(SRC_FILES:.cpp=.o)))

EXE_NAME = bin/furry_circle

all: $(EXE_NAME)

$(EXE_NAME): $(OBJ_FILES)
	$(GPP) $(FLAGS) -o $@ $(INCLUDE) $^ $(LIBS)

dir:
	mkdir -p $(DIR_PATHS)

obj/%.o: src/%.cpp | dir
	$(GPP) $(FLAGS) -c -o $@ $(INCLUDE) $< 

clean:
	rm -rf bin/* obj/* $(EXE_NAME)
	rm -f src/*~ scenes/*~
	rm -f *~

.PHONY: clean
