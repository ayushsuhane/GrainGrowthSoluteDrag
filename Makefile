IDIR = lib

ODIR = obj
SDIR = src

#CC=gcc -W -Wall -ggdb3 -g
CC=gcc -O3 -fopenmp -std=c99 -w -g
CFLAGS=-I ${IDIR}/ -I ${HOME}/gsl/include -L ${HOME}/gsl/lib/$
LFLAGS=-lgsl -lgslcblas

LIBS= -lm

_DEPS = Global.h functions.h 
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS)) 

_OBJ = main.o setup.o postprocess.o phasefield.o misorientation.o Global.o
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))

#all :
#	@echo '$(ODIR)/%.o: ${SDIR}/%.c $(DEPS) $(CC) -c -o ${ODIR}/$@ $< $(CFLAGS)'

$(ODIR)/%.o: ${SDIR}/%.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS) $(LIBS) 

GG: $(OBJ)
	$(CC) -o $@ $^ $(LIBS) $(CFLAGS)  $(LFLAGS)

.PHONY: clean

clean:
	rm -f $(ODIR)/*.o *~ core $(INCDIR)/*~ GG

