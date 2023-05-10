fortran: algencanma_f.o
	$(FC) $(FFLAGS) -L$(BLASSRC)/lib -L$(HSLSRC)/lib -L$(ALGSRC)/lib -lalgencan -lhsl -lblas -shared -o libalgencanma.so algencanma_f.o
	mkdir -p $(BIN)/examples
	# mv algencanma $(BIN)/examples

c: algencanma_c.o
	gcc -g $(FFLAGS) algencanma_c.o -L$(BLASSRC)/lib -L$(HSLSRC)/lib -L$(ALGSRC)/lib -lalgencan -lhsl -lblas -lgfortran -lm -ldl -o algencanma_c
	mkdir -p $(BIN)/examples
	mv algencanma_c $(BIN)/examples

%_f.o: %.f90
	$(FC) -fpic -c $(FFLAGS) -I$(ALGSRC)/inc -o $@ $<

%_c.o: %.c
	gcc -g $(CFLAGS) -I$(ALGSRC)/inc -c -o $@ $< -lm

clean:
	rm -f algencanma_f.o algencanma_c.o

distclean: clean
	rm -f $(BIN)/examples/algencanma
	rm -f $(BIN)/examples/algencanma_c
