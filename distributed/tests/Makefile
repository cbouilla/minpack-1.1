CFLAGS += -Wall -Wextra -I.. -g
LDLIBS += -lm 

CFLAGS += $(shell pkg-config --cflags blas lapacke scalapack-openmpi)
LDLIBS += $(shell pkg-config --libs blas lapacke scalapack-openmpi)

#LDFLAGS += -L.
#LDLIBS += -llapack -lrefblas

all: check

.phony: check clean
	
test_plmdif1.o ls.o: ls.h

test_plmdif1: test_plmdif1.o ls.o ../libminpack.a 

check: test_plmdif1 
	prove -f --exec mpiexec ./test_plmdif1

clean:
	rm -f *.o test_plmdif1
