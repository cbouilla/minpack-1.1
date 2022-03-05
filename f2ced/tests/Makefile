CFLAGS = -Wall -Wextra -g -I..
LDLIBS = -lf2c -lm

all: check

.phony: check clean

test_hybrd1.o test_hybrj1.o test_chkder.o eq.o: eq.h
test_lmdif1.o test_lmder1.o ls.o: ls.h

test_lmdif1: test_lmdif1.o ls.o ../libminpack.a
test_lmder1: test_lmder1.o ls.o ../libminpack.a
test_hybrd1: test_hybrd1.o eq.o ../libminpack.a
test_hybrj1: test_hybrj1.o eq.o ../libminpack.a
test_chkder: test_chkder.o eq.o ../libminpack.a

check: test_lmdif1 test_lmder1 test_hybrd1 test_hybrj1 test_chkder
	prove -f ./test_lmdif1 ./test_lmder1 ./test_hybrd1 ./test_hybrj1 ./test_chkder

clean:
	rm -f *.o test_lmdif1 test_lmder1 test_hybrd1 test_hybrj1 test_chkder
