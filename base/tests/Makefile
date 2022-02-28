CFLAGS = -Wall -Wextra
LDLIBS = -lm

CFLAGS += $(shell pkg-config --cflags blas lapack)
LDLIBS += $(shell pkg-config --libs blas lapack)

all: test bench

bench: speed_lmder1 speed_lmdif1 speed_lmdif1_cross speed_lmder1_cross speed_hybrj1
	./speed_hybrj1
	./speed_lmdif1
	./speed_lmdif1_cross
	./speed_lmder1
	./speed_lmder1_cross


.phony: test clean
	
test_lmdif1: test_lmdif1.o ls.o ../libminpack.a
test_lmder1: test_lmder1.o ls.o ../libminpack.a
test_hybrd1: test_hybrd1.o eq.o ../libminpack.a
test_hybrj1: test_hybrj1.o eq.o ../libminpack.a
test_chkder: test_chkder.o eq.o ../libminpack.a

test: test_lmdif1 test_lmder1 test_hybrd1 test_hybrj1 test_chkder
	prove -f ./test_lmdif1 ./test_lmder1 ./test_hybrd1 ./test_hybrj1 ./test_chkder

speed_hybrj1: speed_hybrj1.o eq.o             ../libminpack.a
speed_lmdif1: speed_lmdif1.o ls.o             ../libminpack.a
speed_lmder1: speed_lmder1.o ls.o             ../libminpack.a
speed_lmdif1_cross: speed_lmdif1_cross.o eq.o ../libminpack.a
speed_lmder1_cross: speed_lmder1_cross.o eq.o ../libminpack.a

clean:
	rm -f *.o test_lmdif1 test_lmder1 test_hybrd1 test_hybrj1 
	rm -f speed_lmdif1 speed_lmder1 speed_lmder1_cross speed_lmdif1_cross speed_hybrj1