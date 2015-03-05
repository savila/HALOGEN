include ../../Makefile.defs

OBJ = define.o common.o io.o neighbors.o correlator.o

#RULES
cute : $(OBJ)

#RULE TO MAKE .o's FROM .c's
$(OBJ) : %.o : %.c
	$(CC) $(CFLAGS) -c $< -o $@ -lm
