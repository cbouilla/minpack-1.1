CFLAGS= -O3 -Wall -Wextra -g -Wno-builtin-declaration-mismatch -I.
LDFLAGS= -g
default: libminpack.a

.phony: clean test

CFLAGS += $(shell pkg-config --cflags blas lapack)
LDLIBS += $(shell pkg-config --libs blas lapack)

OBJS = enorm.o qrfac.o qrsolv.o fdjac2.o fdjac1.o lmpar.o qform.o dogleg.o r1updt.o r1mpyq.o \
       lmdif1.o lmdif.o lmder.o lmder1.o \
       hybrd.o hybrj.o hybrd1.o hybrj1.o

$(OBJS): minpack.h

libminpack.a:  $(OBJS)
	ar r $@ $(OBJS); ranlib $@

check: libminpack.a
	cd tests ; $(MAKE)

clean:
	rm -f *.o libminpack.a