default: libminpack.a

.phony: clean test

FFLAGS=-O3

OBJS = dpmpar.o enorm.o qrfac.o qrsolv.o fdjac2.o fdjac1.o lmpar.o qform.o dogleg.o r1updt.o r1mpyq.o \
       lmdif1.o lmdif.o lmder.o lmder1.o \
       hybrd.o hybrj.o hybrd1.o hybrj1.o \
       chkder.o

libminpack.a:  $(OBJS)
	ar r $@ $(OBJS); ranlib $@

check: libminpack.a
	cd tests ; $(MAKE) check

clean:
	rm -f *.o libminpack.a
	cd tests ; $(MAKE) clean
