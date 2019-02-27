OS            := $(shell uname)

#ifeq ($(OS), Linux)
CXX           := g++
DEBUG         := -g
#endif

CXXFLAGS      := $(DEBUG) -ansi -Wall -Wextra -m64 -O3 -std=c++11
#DEBUGFLAGS := -O3 -Wall -std=c++0x
#CXXFLAGS := $(DEBUGFLAGS)

SOURCEDIR     := ../src
SOURCES       := $(filter-out $(SOURCEDIR)/B0KstMuMu.cc,$(wildcard $(SOURCEDIR)/*.cc))

INCLUDEDIR    := ../interface
INCLUDES      := -I $(INCLUDEDIR)

XERCESINFLAGS := -I $(HOME)/xerces/include
BOOSTINFLAGS  := -I $(HOME)/boost
GSLINCFLAGS   := -I $(HOME)/gsl/include

XERCESFLAGS   := -L $(HOME)/xerces/lib -lxerces-c
BOOSTFLFLAGS  := -L $(HOME)/boost/stage/lib
GSLFLAGS      := -L $(HOME)/gsl/lib -lgsl

ROOTFLAGS     := $(shell root-config --glibs)
ROOFITFLAGS   := $(shell root-config --cflags --libs) -lRooFit -lRooFitCore
ROOFITFLAGS1  := $(shell root-config --cflags --libs) -lMinuit -lRooFit -lRooFitCore -lFoam
#compile class
ROOTCINT := $(shell which rootcint)
CLASS2       := RooBernsteinEffi
CLASSDICT2   := $(CLASS2)Dictionary.cxx

LIBS2 := $(CLASS2).cxx  $(CLASSDICT2)

$(CLASSDICT2): $(CLASS2).h $(CLASS2)LinkDef.h
	$(ROOTCINT) -f $@ -c $^


ExtractYield : $(SOURCEDIR) $(INCLUDEDIR) ExtractYield.cc
	$(CXX) $(CXXFLAGS) -o $@  $(ROOFITFLAGS) $(INCLUDES) $(SOURCES)  -D"ROOFIT" $@.cc
Moment : Moment.cc
	$(CXX) $(CXXFLAG) -o $@ $(ROOTFLAGS) $(ROOFITFLAGS1)  -D"ROOFIT" $@.cc

LOCALCC := $(wildcard *.cc)
clean:
	$(RM) $(LOCALCC:%.cc=%)  $(CLASSDICT2)  Macros_C.*
