include ../Makefile.defs

EXEC   = 2LPT-HALOGEN

###### DEFINE OBJECT CATEGORIES
CUTE_OBJ    = $(addprefix cute.src/,common.o define.o io.o neighbors.o correlator.o)
2LPT_OBJ   = power.o allvars.o save.o read_param.o  read_glass.o \
	$(addprefix nrsrc/,nrutil.o qromb.o polint.o trapzd.o) 
HLGN_OBJ   = read_snapshot.o place_halos.o populate_mass_function.o
FIT_OBJ    = $(HLGN_OBJ) correlate.o
ALL_OBJ    = $(CUTE_OBJ) $(2LPT_OBJ) $(FIT_OBJ)

####### SPECIFY ALL LIBRARIES
FFTW_LIB =  $(FFTW_LIBS)  -ldrfftw_mpi -ldfftw_mpi -ldrfftw -ldfftw
GSL_LIB = $(GSL_LIBS) -lgsl  -lgslcblas
OTHER_LIB   =   -lm  $(MPICHLIB)


2LPT_INCL   = allvars.h proto.h  nrsrc/nrutil.h
ALL_INCL = $(2LPT_INCL) $(FFTW_INCL) $(GSL_INCL)
ALL_LIB  = $(FFTW_LIB) $(GSL_LIB) $(OTHER_LIB)

###### RULES
$(EXEC): main.c $(2LPT_OBJ) $(HLGN_OBJ) 
	$(CCMPI) $(CFLAGS) $(DEFS) main.c $(2LPT_OBJ) $(HLGN_OBJ) $(ALL_INCL) $(ALL_LIB)  -o  $(EXEC)  

$(2LPT_OBJ): %.o : %.c $(2LPT_INCL)
	 $(CCMPI) $(CFLAGS) -DONLY_2LPT $(FFTW_INCL) -o $@ -c $< $(FFTW_LIB) $(OTHER_LIB) 

2LPT: 2LPT.c $(2LPT_OBJ) $(2LPT_INCL)
	$(CCMPI) $(CFLAGS) -DONLY_2LPT 2LPT.c $(2LPT_OBJ) $(FFTW_INCL) $(FFTW_LIB) $(OTHER_LIB)  -o  2LPT  


halogen:	halogen.c $(HLGN_OBJ) $(HLGN_OBJ:o=h)
	$(CC) $(CFLAGS) $(DEFS) halogen.c $(HLGN_OBJ) $(OTHER_LIB)

fit:            fit.c halogen.c $(FIT_OBJ) $(FIT_OBJ:o=h)
	$(CC) $(CFLAGS) $(DEFS) $(GSL_INCL) -o fit $(CUTE_OBJ) $(FIT_OBJ) fit.c $(GSL_LIB) $(OTHER_LIB)

clean:
	-rm -f $(ALL_OBJ) $(EXEC) halogen fit 2LPT

$(HLGN_OBJ) : %.o : %.c %.h
	$(CC) $(CFLAGS) $(DEFS) -c $< $(OTHER_LIB)

correlate.o: %.o : cute %.c %.h
	$(CC) $(CFLAGS) $(DEFS) -c correlate.c correlate.h $(OTHER_LIB)

cute:                   
	$(MAKE) -C cute.src cute

