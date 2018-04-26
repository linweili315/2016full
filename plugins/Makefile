OS            := $(shell uname)

ifeq ($(OS), Darwin)
CXX           := clang++
DEBUG         :=
endif

ifeq ($(OS), Linux)
CXX           := g++
DEBUG         := -g
endif

CXXFLAGS      := $(DEBUG) -ansi -Wall -Wextra -m64 -O3

SOURCEDIR     := ../src
SOURCES       := $(filter-out $(SOURCEDIR)/B0KstMuMu.cc,$(wildcard $(SOURCEDIR)/*.cc)  $(SOURCEDIR)/RooDataHist.cxx)

INCLUDEDIR    := ../interface
INCLUDES      := -I $(INCLUDEDIR)

XERCESINFLAGS := -I $(HOME)/xerces/include
BOOSTINFLAGS  := -I $(HOME)/boost
GSLINCFLAGS   := -I $(HOME)/gsl/include

XERCESFLAGS   := -L $(HOME)/xerces/lib -lxerces-c
BOOSTFLFLAGS  := -L $(HOME)/boost/stage/lib
GSLFLAGS      := -L $(HOME)/gsl/lib -lgsl

ROOTFLAGS     := `root-config --cflags --glibs` -lMinuit
ROOFITFLAGS   := $(ROOTFLAGS) -lRooFit -lRooFitCore -lFoam


ExtractYield: $(SOURCEDIR) $(INCLUDEDIR) ExtractYield.cc
	$(CXX) $(CXXFLAGS) -o $@ $(ROOFITFLAGS) $(INCLUDES) $(SOURCES) -D"ROOFIT" $@.cc

LOCALCC := $(wildcard *.cc)
clean:
	$(RM) $(LOCALCC:%.cc=%) Macros_C.*
