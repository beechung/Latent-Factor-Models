
all: c_funcs

c_funcs: clean
	R CMD SHLIB C/util.c C/MCEM_EStep.c C/MCEM_EStep2.c C/MC_predict.c C/MC_predict2.c -o c_funcs.so

clean:
	rm -f C/*.o *.so

