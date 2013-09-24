include ../../Makefile.inc


gsl: eval.o solve_cubic.o 
	ar cr ../../lib/libgsl.a eval.o solve_cubic.o


clean:
	rm -f *.o

mingw_clean:
	del *.o
