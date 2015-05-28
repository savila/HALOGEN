include ../Makefile.defs
EXEC   = 2LPT-HALOGEN

###### DEFINE OBJECT CATEGORIES
CUTE_OBJ    = $(addprefix cute.src/,common.o define.o io.o neighbors.o correlator.o)

# objects that can be either 2lpt-only or not
2LPT_ONLY_OBJ = read_param.o read_glass.o save.o place_halos.o
2LPT_OBJ   = allvars.o power.o $(2LPT_ONLY_OBJ) $(addprefix nrsrc/,nrutil.o qromb.o polint.o trapzd.o) 

HLGN_OBJ   = read_snapshot.o populate_mass_function.o
HLGN_ONLY_OBJ   = populate_mass_function.o
FIT_OBJ    = $(HLGN_OBJ) correlate.o
ALL_OBJ    = $(CUTE_OBJ) $(2LPT_OBJ) $(FIT_OBJ)

####### SPECIFY ALL LIBRARIES
FFTW_LIB =  $(FFTW_LIBS)  -ldrfftw_mpi -ldfftw_mpi -ldrfftw -ldfftw
GSL_LIB = $(GSL_LIBS) -lgsl  -lgslcblas
OTHER_LIB   =   -lm  $(MPICHLIB)


2LPT_INCL   = allvars.h proto.h  nrsrc/nrutil.h
ALL_INCL = $(FFTW_INCL) $(GSL_INCL)
ALL_LIB  = $(FFTW_LIB) $(GSL_LIB) $(OTHER_LIB)

###### RULES
all: $(EXEC) halogen fit
	-rm $(2LPT_ONLY_OBJ)
	$(MAKE) 2LPT

$(EXEC): main.c $(2LPT_OBJ) $(HLGN_OBJ) 
	$(CCMPI) $(CFLAGS) $(DEFS) $(2LPT_OBJ) $(HLGN_ONLY_OBJ) main.c $(ALL_INCL) $(ALL_LIB)  -o  $(EXEC)  

$(2LPT_OBJ): %.o : %.c $(2LPT_INCL)
	 $(CCMPI) $(CFLAGS) $(DEFS) $(2LPTOPTS) $(ALL_INCL) -o $@ -c $< $(ALL_LIB) $(OTHER_LIB) 

2LPT: DEFS+=-DONLY_2LPT
2LPT: 2LPT.c $(2LPT_OBJ) $(2LPT_INCL)
	$(CCMPI) $(CFLAGS) $(DEFS) $(2LPTOPTS) $(2LPT_OBJ) 2LPT.c $(ALL_INCL) $(ALL_LIB)  -o  2LPT  
	-rm $(2LPT_ONLY_OBJ)

halogen:	halogen.c $(HLGN_OBJ) $(HLGN_OBJ:o=h)
	$(CC) $(CFLAGS) $(DEFS) halogen.c $(HLGN_OBJ) $(OTHER_LIB) -o halogen

fit:            fit.c halogen.c $(FIT_OBJ) $(FIT_OBJ:o=h)
	$(CC) $(CFLAGS) $(DEFS) $(GSL_INCL) -o fit $(CUTE_OBJ) $(FIT_OBJ) fit.c $(GSL_LIB) $(OTHER_LIB)

clean:
	-rm -f $(ALL_OBJ) $(EXEC) halogen fit 2LPT

$(HLGN_ONLY_OBJ) : %.o : %.c %.h
	$(CCMPI) $(CFLAGS) $(DEFS) -c $< $(ALL_INCL) $(OTHER_LIB) 

correlate.o: %.o : cute %.c %.h
	$(CC) $(CFLAGS) $(DEFS) -c correlate.c $(OTHER_LIB)

cute:                   
	$(MAKE) -C cute.src cute

