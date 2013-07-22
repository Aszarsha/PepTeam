#Conditionally declare CC and CXX so they can be overwritten from command line
CC?=gcc
CXX?=g++
CFLAGS=-Wall -march=native -Ofast
CXXFLAGS=$(CFLAGS) -std=c++1y
LDFLAGS=
#CPPFILES := $(wildcard src/*.cpp)
#OBJFILES := $(addprefix obj/,$(notdir $(CPP_FILES:.cpp=.o)))

all: FastIdx PepTree PepteamMap PepteamProfile

FastIdx: bindir obj/FastIdx.o obj/FastIdx_drv.o
	$(CXX) $(LDFLAGS) obj/FastIdx.o obj/FastIdx_drv.o -o bin/FastIdx

PepTree: bindir obj/FastIdx.o obj/PepTree.o obj/PepTree_drv.o
	$(CXX) $(LDFLAGS) obj/FastIdx.o obj/PepTree.o obj/PepTree_drv.o -o bin/PepTree

PepteamMap: bindir obj/FastIdx.o obj/PepTree.o obj/PepteamMap.o
	$(CXX) $(LDFLAGS) obj/FastIdx.o obj/PepTree.o obj/PepteamMap.o -o bin/PepteamMap

PepteamProfile: bindir obj/FastIdx.o obj/PepTree.o obj/PepteamProfile.o
	$(CXX) $(LDFLAGS) obj/FastIdx.o obj/PepTree.o obj/PepteamProfile.o -o bin/PepteamProfile

obj/%.o: src/%.cpp objdir
	$(CXX) $(CXXFLAGS) -c -o $@ $<

objdir:
	mkdir -p obj

bindir:
	mkdir -p bin

clean:
	rm -rf obj bin
