
all: R

R:
	R CMD SHLIB -o ../../../lib/R_hierarchical.so -c util.cpp dag.c cube.cpp betabinomial.cpp Rinterface.cpp

clean:
	rm -f *.o ../../../lib/R_hierarchical.so
