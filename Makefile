#Conditionally declare CC and CXX so they can be overwritten from command line
CC=clang
CXX=clang++
CFLAGS=-Wall -march=native -Ofast
CXXFLAGS=$(CFLAGS) -std=c++14
LDFLAGS=-lboost_program_options
#CPPFILES := $(wildcard src/*.cpp)
#OBJFILES := $(addprefix obj/,$(notdir $(CPP_FILES:.cpp=.o)))

all: FastIdx PepTree PepteamMap PepteamProfile

FastIdx: obj/FastIdx.o obj/FastIdx_drv.o | bindir
	$(CXX) $(LDFLAGS) obj/FastIdx.o obj/FastIdx_drv.o -o bin/FastIdx

PepTree: obj/FastIdx.o obj/PepTree.o obj/PepTree_drv.o | bindir
	$(CXX) $(LDFLAGS) obj/FastIdx.o obj/PepTree.o obj/PepTree_drv.o -o bin/PepTree

PepteamMap: obj/FastIdx.o obj/PepTree.o obj/PepteamMap.o | bindir
	$(CXX) $(LDFLAGS) obj/FastIdx.o obj/PepTree.o obj/PepteamMap.o -o bin/PepteamMap

PepteamProfile: obj/FastIdx.o obj/PepTree.o obj/PepteamProfile.o | bindir
	$(CXX) $(LDFLAGS) obj/FastIdx.o obj/PepTree.o obj/PepteamProfile.o -o bin/PepteamProfile

obj/%.o: src/%.cpp | objdir
	$(CXX) $(CXXFLAGS) -c -o $@ $<

objdir:
	mkdir -p obj

bindir:
	mkdir -p bin

clean:
	rm -rf obj bin
