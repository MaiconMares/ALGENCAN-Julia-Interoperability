LDFLAGS := -lalgencan -lcutest -lhsl -lblas
SIF := ELFUN EXTER GROUP RANGE
MAINFILE := $(shell ./get.type $(PROBNAME))

%.o: %.f
	$(FC) $(FFLAGS) -I$(ALGSRC)/inc -I$(HSLSRC)/inc -c -o $@ $<

%.o: %.f90
	$(FC) $(FFLAGS) -I$(ALGSRC)/inc -I$(HSLSRC)/inc -c -o $@ $<

all: sifdec algencan

sifdec:
	$(SIFDECODE)/bin/sifdecoder $(MASTSIF)/$(PROBNAME).SIF
	$(FC) -c $(FFLAGS) $(addsuffix .f,$(SIF))
	rm -f AUTOMAT.d

algencan: $(MAINFILE).o $(addsuffix .o,$(SIF))
	$(FC) $^ -L$(CUTEST)/objects/$(MYARCH)/double -L$(ALGSRC)/lib -L$(BLASSRC)/lib -L$(HSLSRC)/lib -o $@ $(LDFLAGS)
	mkdir -p $(BIN)/cutest
	mv OUTSDIF.d $(BIN)/cutest
	mv $@ $(BIN)/cutest

clean:
	rm -f *.o
	rm -f $(addsuffix .f,$(SIF))

distclean: clean
	rm -f $(BIN)/cutest/OUTSDIF.d
	rm -f $(BIN)/cutest/algencan
