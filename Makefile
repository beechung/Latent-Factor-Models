
all: lib c_funcs arslogistic

pagerank: src/C/pagerank.c src/C/util.h
	R CMD SHLIB src/C/pagerank.c -o lib/pagerank.so

c_funcs: src/C/util.h src/C/factor_model_util.h src/C/pagerank.c src/C/util.c src/C/factor_model_util.c src/C/factor_model_util2.cpp src/C/hierarchical.c src/C/hierarchical.h src/C/factor_model_multicontext.c
	R CMD SHLIB src/C/util.c src/C/factor_model_util.c src/C/pagerank.c src/C/hierarchical.c src/C/factor_model_multicontext.c src/C/factor_model_util2.cpp -o lib/c_funcs.so

arslogistic: src/arslogistic/utilR.h src/arslogistic/arms.c src/arslogistic/arsspline.c src/RLFM-ars-logistic/C/util.c src/RLFM-ars-logistic/C/MCEM_EStep_logistic.c src/RLFM-ars-logistic/C/logistic.c
	R CMD SHLIB src/arslogistic/arms.c src/arslogistic/arsspline.c src/RLFM-ars-logistic/C/util.c src/RLFM-ars-logistic/C/MCEM_EStep_logistic.c src/RLFM-ars-logistic/C/logistic.c -o lib/arslogistic.so

clean:
	rm -f src/C/*.o src/C/*.so src/arslogistic/*.o src/RLFM-ars-logistic/C/*.o lib/*.so

lib:
	mkdir lib
