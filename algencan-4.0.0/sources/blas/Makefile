BLAS := dgemm.o dgemv.o dtpmv.o dtpsv.o idamax.o lsame.o xerbla.o

ifneq ("$(wildcard ./lib/libblas.a)","")

all:

else

all: $(BLAS)
	mkdir -p ./lib/
	$(AR) rcs ./lib/libblas.a $(BLAS)

endif

%.o: %.f
	$(FC) $(FFLAGS) -c -o $@ $<

clean:
	rm -f *.o

distclean: clean
	rm -f lib/libblas.a
