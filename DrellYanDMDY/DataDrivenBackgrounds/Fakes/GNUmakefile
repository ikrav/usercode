ARCH         := $(shell $(ROOTSYS)/bin/root-config --arch)

ROOTCFLAGS    = $(shell ${ROOTSYS}/bin/root-config --cflags)
ROOTLIBS      = $(shell ${ROOTSYS}/bin/root-config --libs)
ROOTGLIBS     = $(shell ${ROOTSYS}/bin/root-config --glibs)
CERNLIB       = /cern/pro/lib

CXX = g++
LD  = g++

CXXFLAGS      = -g -O -Wall -fPIC
#CXXFLAGS      = -O -Wall -fPIC
LDFLAGS       = -g
SOFLAGS       = -shared

CXX           += -I/afs/cern.ch/cms/sw/slc5_ia32_gcc434/lcg/roofit/5.26.00-cms5/include
CXXFLAGS      += $(ROOTCFLAGS)

LIBS           = $(ROOTLIBS) 

NGLIBS         = $(ROOTGLIBS) 
NGLIBS        += -lMinuit -lHtml -lTreePlayer -lEG



GLIBS          = $(filter-out -lNew, $(NGLIBS))
CLIBS       = -L$(CERNLIB) -lpacklib -lmathlib -lkernlib -lpawlib -lg2c -lnsl

.SUFFIXES: .cc,.C


# ================================================================================
MakePlots :  MakePlots.o  
# -------------------------

	$(LD) $(LDFLAGS) -o $@ $^ $(GLIBS)
MakePlots.o:  DrawPlots.cc  
# -------------------------
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -rf MakePlots.o  MakePlots 



