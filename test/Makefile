FC=gfortran
CC=gcc

INC=-I../include
CCFLAGS=-std=c99 -fPIC -g -Wall ${INC}
FCFLAGS=-fPIC -g -Wall -fbacktrace ${INC}
CCLINK=-lgfortran -lm

test: main.c quasif.o quasic.o mief.o miec.o arraycontractions.o
	${CC} ${CCFLAGS} ${CCLINK} -o $@ $^
	./test
	rm test *.o

miec.o: ../src/mie.c arraycontractions.o
	${CC} ${CCFLAGS} -c $< -o $@

quasic.o: ../src/quasi.c arraycontractions.o
	${CC} ${CCFLAGS} -c $< -o $@

arraycontractions.o: ../src/arraycontractions.c
	${CC} ${CCFLAGS} -c $< -o $@

mief.o: mie.f90
	${FC} ${FCFLAGS} -c $< -o $@

quasif.o: quasi.f90
	${FC} ${FCFLAGS} -c $< -o $@

#  Make sure we don't remove this by accident if interrupted at the wrong time.
.PRECIOUS: Makefile
