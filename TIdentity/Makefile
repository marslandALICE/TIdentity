##### Temporary makefile
##### A. Rustamov
PREFIX     = ./
NAME       = TIdentity
SHLIB      = $(PREFIX)lib/lib$(NAME).so

ROOTLIBS   = $(shell root-config --libs  | sed 's/-pthread //') -lMinuit -lMatrix
ROOTGLIBS  = $(shell root-config --glibs | sed 's/-pthread //') -lMinuit -lMatrix
LIBS       = $(ROOTLIBS)

INCLUDE    = -I$(ROOTSYS)/include \
             -I$(CINTSYSDIR) \
             -I$(PREFIX)                  

INCDIR     = $(PREFIX)
SRCDIR     = $(PREFIX)
VPATH      = $(SRCDIR) $(INCDIR) 

DICT       = $(NAME)Dict
LINKDEF    = $(NAME)LinkDef

CXX        = g++
LD         = g++

ROOTCFLAGS = $(shell root-config --cflags | sed 's/-pthread //') 
CXXFLAGS   = -g -Wall -fPIC
CXXFLAGS  += $(ROOTCFLAGS)
LDFLAGS    = -g 

SOURCES    = $(wildcard *.cc)
OBJECTS    = $(filter-out $(DICT).o,$(SOURCES:.cc=.o)) $(DICT).o 
HEADER     = TIdentity2D.h    \
             TIdentityFunctions.h \
             $(LINKDEF).h        

.PHONY: all clean 

all: 	
	$(MAKE) Makefile $(SHLIB)

$(SHLIB): $(OBJECTS) $(HEADER)
	$(CXX) -shared $(LDFLAGS) -o $(SHLIB) $(OBJECTS) $(INCLUDE) 

$(DICT).C: $(HEADER)
	@echo "Generating dictionary $(NAME)Dict.C"
	$(ROOTSYS)/bin/rootcint -f $(DICT).C -c $(INCLUDE) $(HEADER)

clean:
	rm -f *.o *Dict.* $(SHLIB) core *~

.c.o:
	cc $(CFLAGS) $(ANSI) $(INCLUDE) -c  $*.c

.cpp.o:
	$(CXX)   $(CXXFLAGS) $(INCLUDE) -c  $*.C

.C.o:
	$(CXX)   $(CXXFLAGS) $(INCLUDE) -c  $*.C

.cc.o:
	$(CXX)   $(CXXFLAGS) $(INCLUDE) -c  $*.cc

.f.o:
	$(FF)    $(FFLAGS)   $(INCLUDE) -c  $*.f












