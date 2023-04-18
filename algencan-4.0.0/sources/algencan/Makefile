ALGENCAN := algencan.o feasgencan.o gencan.o lss.o memev.o newtkkt.o

all: lib

lib: $(ALGENCAN)
	mkdir -p ./inc/
	mkdir -p ./lib/
	$(AR) rcs ./lib/libalgencan.a $(ALGENCAN)
ifeq ($(shell test -e ./bmalgencan.mod && echo -n yes),yes)
	mv -f bmalgencan.mod ./inc/
endif
ifeq ($(shell test -e ./bmgencan.mod && echo -n yes),yes)
	mv -f bmgencan.mod ./inc/
endif
ifeq ($(shell test -e ./bmfeasalgencan.mod && echo -n yes),yes)
	mv -f bmfeasgencan.mod ./inc/
endif

%.o: %.f90
	$(FC) $(FFLAGS) -I$(HSLSRC)/inc -I./inc -c -o $@ $<

algencan.o: algencan.f90 lss.o newtkkt.o gencan.o feasgencan.o

feasgencan.o: feasgencan.f90 lss.o memev.o gencan.o

gencan.o: gencan.f90 lss.o

newtkkt.o: newtkkt.f90 lss.o memev.o

clean:
	rm -f *.o *.mod

distclean: clean
	rm -f ./lib/libalgencan.a
	rm -f ./inc/bmalgencan.mod ./inc/bmgencan.mod ./inc/bmfeasgencan.mod
