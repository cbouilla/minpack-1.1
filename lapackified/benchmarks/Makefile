CFLAGS = -Wall -Wextra -I..
LDLIBS = -lm

CFLAGS += $(shell pkg-config --cflags blas lapack)
LDLIBS += $(shell pkg-config --libs blas lapack)

all: bench

bench: speed_lmder1 speed_lmdif1 speed_lmdif1_cross speed_lmder1_cross speed_hybrj1
	./speed_hybrj1
	./speed_lmdif1
	./speed_lmdif1_cross
	./speed_lmder1
	./speed_lmder1_cross

.phony: bench clean

speed_hybrj1: speed_hybrj1.o eq.o             ../libminpack.a
speed_lmdif1: speed_lmdif1.o ls.o             ../libminpack.a
speed_lmder1: speed_lmder1.o ls.o             ../libminpack.a
speed_lmdif1_cross: speed_lmdif1_cross.o eq.o ../libminpack.a
speed_lmder1_cross: speed_lmder1_cross.o eq.o ../libminpack.a

clean:
	rm -f *.o
	rm -f speed_lmdif1 speed_lmder1 speed_lmder1_cross speed_lmdif1_cross speed_hybrj1
