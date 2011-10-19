
all: c_funcs

c_funcs: clean
	R CMD SHLIB ../C/util.c C/EM.cpp C/R_interface.cpp -o c_funcs.so

clean:
	rm -f ../C/*util.o C/*.o *.so

