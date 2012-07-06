ARCH         := $(shell $(ROOTSYS)/bin/root-config --arch)

ROOTCFLAGS    = $(shell ${ROOTSYS}/bin/root-config --cflags)
ROOTLIBS      = $(shell ${ROOTSYS}/bin/root-config --libs)
ROOTGLIBS     = $(shell ${ROOTSYS}/bin/root-config --glibs)
CERNLIB       = /cern/pro/lib

CXX = g++
LD  = g++

CXXFLAGS      = -g -O2 -Wall -fPIC -pg -D__noRooFit
#CXXFLAGS      = -O2 -Wall -fPIC
LDFLAGS       = -g 
SOFLAGS       = -shared

INC           = ../../Include/
CXX           += -I${INC} 
CXX     	+= -I/usr/include/python2.6
CXXFLAGS      += $(ROOTCFLAGS)

SOFLAGS 	= -shared
ADDROOTLIBS	= -lHistPainter -lASImage

LIBS           = $(ROOTLIBS)
SLIB 		= -L/usr/lib64/ -lboost_python
NGLIBS         = $(ROOTGLIBS)
#NGLIBS        += -lMinuit -lHtml -lTreePlayer -lEG

  

GLIBS          = $(filter-out -lNew, $(NGLIBS))
CLIBS       = -L$(CERNLIB) -lpacklib -lmathlib -lkernlib -lpawlib -lg2c -lnsl

#OBJDIR		= ../../obj
#SRCDIR		= ../../src

#vpath %.cc ../../src
#vpath %.hh ../../include

auxOBJS=CPlot.o PUReweight.o

default: eMuBkgExe

# ================================================================================
run:	bkg.so eMuBkgExe
# -------------------------
eMuBkgExe: main.o eMu.o	CmdLineOpts.o ${auxOBJS}
	$(LD) $(LDFLAGS) -o $@ $^ $(GLIBS) 
bkg.so: eMu_pyConf.o eMu.o
	$(LD) $(LDFLAGS) $(SOFLAGS) $(ADDROOTLIBS) $(SLIB) $(GLIBS) -o $@ $^
%.o: %.cc %.hh
	$(CXX) $(CXXFLAGS) -c $<
${OBJDIR}/%.o:	%.cc
	$(CXX) $(CXXFLAGS) -c $< -o $@
clean:
	rm -rf eMuBkgExe eMu.o bkg.so main.o eMu_pyConf.o CmdLineOpts.o


PUReweight.o: ${INC}PUReweight.cc ${INC}PUReweight.hh ${INC}MyTools.hh CPlot.o
	${CXX} ${CXXFLAGS} -c $< -o $@

CPlot.o: ${INC}CPlot.cc ${INC}CPlot.hh
	${CXX} ${CXXFLAGS} -c $< -o $@
