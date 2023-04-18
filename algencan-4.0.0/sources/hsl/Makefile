HSL := fakemetis.o hsl_zd11d.o mc21ad.o mc47ad.o mc64ad.o hsl_ma57d.o ma57ad.o mc34ad.o mc59ad.o mc71ad.o

%.o: %.f
	$(FC) $(FFLAGS) -c -o $@ $<

%.o: %.f90
	$(FC) $(FFLAGS) -c -o $@ $<

ifneq ("$(wildcard ./lib/libhsl.a)","")

all:

else

all: $(HSL)
	mkdir -p ./inc/
	mkdir -p ./lib/
	$(AR) rcl ./lib/libhsl.a $(HSL)
	mv hsl_ma57_double.mod ./inc/

endif

fakemetis.o: fakemetis.f

hsl_ma57d.o: hsl_ma57d.f90 ma57ad.o mc71ad.o hsl_zd11d.o

hsl_zd11d.o: hsl_zd11d.f90

ma57ad.o: ma57ad.f fakemetis.o mc47ad.o mc64ad.o mc71ad.o

mc21ad.o: mc21ad.f

mc34ad.o: mc34ad.f

mc47ad.o: mc47ad.f mc34ad.o mc59ad.o

mc59ad.o: mc59ad.f

mc64ad.o: mc64ad.f mc21ad.o

mc71ad.o: mc71ad.f

clean:
	rm -f *.o
	rm -f *.mod

distclean: clean
	rm -f ./lib/libhsl.a
	rm -f ./inc/hsl_ma57_double.mod
