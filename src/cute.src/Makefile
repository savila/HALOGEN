
include ../../Makefile.defs

# COMPILER AND OPTIONS
COMPCPU = gcc
OPTCPU = -Wall -O3 -fopenmp -g

#INCLUDES AND LIBRARIES
LIBCPU = -lm

#.o FILES
#CUTE
DEF = define.o
COM = common.o
IO = io.o
NEIGH = neighbors.o
CORR = correlator.o

OFILES = $(DEF) $(COM) $(IO) $(NEIGH) $(CORR) #$(MAIN)

#RULES
cute : $(OFILES)

#RULE TO MAKE .o's FROM .c's
$(DEF) : define.c
	$(COMPCPU) $(OPTCPU) $(LIBCPU) -c $< -o $@
$(COM) : common.c
	$(COMPCPU) $(OPTCPU) $(LIBCPU) -c $< -o $@
$(CORR) : correlator.c
	$(COMPCPU) $(OPTCPU) $(LIBCPU) -c $< -o $@
$(IO) : io.c
	$(COMPCPU) $(OPTCPU) $(LIBCPU) -c $< -o $@
$(NEIGH) : neighbors.c
	$(COMPCPU) $(OPTCPU) $(LIBCPU) -c $< -o $@
