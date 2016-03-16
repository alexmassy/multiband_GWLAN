#------------------------------------------------------------
#
# When you adapt this makefile to compile your CPLEX programs
# please copy this makefile and set CPLEXDIR and CONCERTDIR to
# the directories where CPLEX and CONCERT are installed.
#
#------------------------------------------------------------

ILOGSTUDIODIR = /home/alessio/CPLEX
CPLEXDIR      = $(ILOGSTUDIODIR)/cplex
CONCERTDIR    = $(ILOGSTUDIODIR)/concert
# ---------------------------------------------------------------------
# Compiler selection 
# ---------------------------------------------------------------------

CCC = g++ -O0

# ---------------------------------------------------------------------
# Link options and libraries
# ---------------------------------------------------------------------

CPLEXBINDIR   = $(CPLEXDIR)/bin/$(BINDIST)
CPLEXLIBDIR   = $(CPLEXDIR)/lib/x86-64_sles10_4.1/static_pic
CONCERTLIBDIR = $(CONCERTDIR)/lib/x86-64_sles10_4.1/static_pic

CCLNFLAGS = -L$(CPLEXLIBDIR) -lilocplex -lcplex -L$(CONCERTLIBDIR) -lconcert -lm -pthread


all:
	make -f Makefile all_cpp

execute: all
	make -f Makefile execute_cpp

CONCERTINCDIR = $(CONCERTDIR)/include
CPLEXINCDIR   = $(CPLEXDIR)/include

EXDIR         = $(pwd)
EXINC         = $(EXDIR) #/include
EXDATA        = $(EXDIR) #/data
EXSRCCPP      = $(EXDIR) #/src/cpp

CCFLAGS = -g2 $(CCOPT) -I$(CPLEXINCDIR) -I$(CONCERTINCDIR) 


#------------------------------------------------------------
#  make all      : to compile the examples. 
#  make execute  : to compile and execute the examples.
#------------------------------------------------------------



CPP_EX = clean gwlan_ex inst_mak 

all_cpp: $(CPP_EX)


execute_cpp: $(CPP_EX)
		./gwlan_ex


# ------------------------------------------------------------

clean :
	rm -rf *.o *~ *.class
	rm -rf $(CPP_EX)
	rm -rf *.mps *.ord *.sos *.lp *.sav *.net *.msg *.log *.clp

# ------------------------------------------------------------


gwlan_ex: algo_main.o gwlan_lib.o
	$(CCC) $(CCFLAGS) algo_main.o gwlan_lib.o -o gwlan_ex $(CCLNFLAGS)
algo_main.o: $(EXSRCCPP)algo_main.cc 
	$(CCC) -c $(CCFLAGS) $(EXSRCCPP)algo_main.cc 
gwlan_lib.o: $(EXSRCCPP)gwlan_lib.cc
	g++ -std=c++11 -c $(CCFLAGS) $(EXSRCCPP)gwlan_lib.cc 

inst_mak: instance-maker.cc
	g++ -std=c++0x -o inst-maker -I. ConfigFile.cpp instance-maker.cc

