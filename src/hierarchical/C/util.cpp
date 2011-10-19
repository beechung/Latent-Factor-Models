#include "util.h"
#include <math.h>


void error(const char *fmt, ...)
{
	char buf[BUFSIZ];
	va_list ap;
	va_start(ap,fmt);
	vsprintf(buf,fmt,ap);
	va_end(ap);
	fprintf(stderr, buf);
    exit(-1);
}

void error2(const char *filename, int lineno, const char *fmt, ...){
	char buf[BUFSIZ];
	va_list ap;
	va_start(ap,fmt);
	vsprintf(buf,fmt,ap);
	va_end(ap);
	fprintf(stderr, "ERROR in file %s at line %d: ", filename, lineno);
	fprintf(stderr, buf);
	fprintf(stderr, "\n");
    exit(-1);
}

void info(FILE* fp, const int verbose, const char* fmt, ...){
	char buf[BUFSIZ];
	va_list ap;
	va_start(ap,fmt);
	vsprintf(buf,fmt,ap);
	va_end(ap);
	int n = (verbose-1)*3;
	if(n<0) n=0;
	std::string str = indent(buf, n);
	fprintf(fp, str.c_str());
}

CReverseIndex::CReverseIndex(const int* indexArray, const int arrayLength, const int nIndexValues){

    revIndex = Malloc(int*, nIndexValues);
    num      = Malloc(int, nIndexValues);
    index_space = Malloc(int, arrayLength);
    size = nIndexValues;
    
    for(int i=0; i<nIndexValues; i++) num[i] = 0;
    
    for(int k=0; k<arrayLength; k++){
        int i = indexArray[k];
#ifdef SANITY_CHECK
        CHK_C_INDEX(i, nIndexValues);
#endif
        num[i]++;
    }
    
    revIndex[0] = index_space;
    for(int i=1; i<nIndexValues; i++){
        revIndex[i] = &(revIndex[i-1][num[i-1]]);
    }
    
    int *current = Malloc(int, nIndexValues);
    for(int i=0; i<nIndexValues; i++) current[i] = 0;

    for(int k=0; k<arrayLength; k++){
        int i = indexArray[k];
        revIndex[i][current[i]] = k;
        current[i]++;
    }
#ifdef SANITY_CHECK
    for(int i=0; i<nIndexValues; i++){
        if(current[i] != num[i]) STOP_HERE("sanity check error");
        for(int j=0; j<num[i]; j++){
            if(i != indexArray[revIndex[i][j]]) STOP_HERE("sanity check error");
        }
    }
#endif
    free(current);
}

CReverseIndex::~CReverseIndex(){
    free(revIndex);
    free(num);
    free(index_space);
}

void tokenize(const std::string& str, std::vector<std::string>& tokens, const std::string& delimiters)
{
    // Skip delimiters at beginning.
    std::string::size_type lastPos = str.find_first_not_of(delimiters, 0);
    // Find first "non-delimiter".
    std::string::size_type pos     = str.find_first_of(delimiters, lastPos);

    tokens.clear();
    
    while (std::string::npos != pos || std::string::npos != lastPos)
    {
        // Found a token, add it to the vector.
        tokens.push_back(str.substr(lastPos, pos - lastPos));
        // Skip delimiters.  Note the "not_of"
        lastPos = str.find_first_not_of(delimiters, pos);
        // Find next "non-delimiter"
        pos = str.find_first_of(delimiters, lastPos);
    }
}

void print_matrix(double **matrix, const int nrow, const int ncol, FILE* fp, const char* sep, const char* format){
    for(int i=0; i<nrow; i++){
        for(int j=0; j<ncol; j++){
            if(j>0) fprintf(fp, sep);
            fprintf(fp, format, matrix[i][j]);
        }
        fprintf(fp, "\n");
    }
}

void print_row_vector(double *vector, const int num, FILE* fp, const char* sep, const char* format){
    for(int i=0; i<num; i++){
        if(i>0) fprintf(fp, sep);
        fprintf(fp, format, vector[i]);
    }
}

void print_col_vector(double *vector, const int num, FILE* fp, const char* sep, const char* format){
    for(int i=0; i<num; i++){
        fprintf(fp, format, vector[i]);
        fprintf(fp, "\n");
    }
}


void print_matrix(int **matrix, const int nrow, const int ncol, FILE* fp, const char* sep, const char* format){
    for(int i=0; i<nrow; i++){
        for(int j=0; j<ncol; j++){
            if(j>0) fprintf(fp, sep);
            fprintf(fp, format, matrix[i][j]);
        }
        fprintf(fp, "\n");
    }
}

void print_row_vector(int *vector, const int num, FILE* fp, const char* sep, const char* format){
    for(int i=0; i<num; i++){
        if(i>0) fprintf(fp, sep);
        fprintf(fp, format, vector[i]);
    }
}

void print_col_vector(int *vector, const int num, FILE* fp, const char* sep, const char* format){
    for(int i=0; i<num; i++){
        fprintf(fp, format, vector[i]);
        fprintf(fp, "\n");
    }
}



std::string indent(const char* input, int n){
	int len = strlen(input);
	std::string out = "";
	char *buf = (char*) malloc((len+1) * sizeof(char));
	char *space = (char*) malloc((n+1) * sizeof(char));
	for(int i=0; i<n; i++) space[i] = ' ';
	space[n] = '\0';
	int k=0;
	for(int i=0; input[i] != '\0'; i++){
		buf[k] = input[i];
		k++;
		if(input[i] == '\n'){
			buf[k] = '\0';
			if(k > 1) out += space;
			out += buf;
			k=0;
		}
	}
	buf[k] = '\0';
	if(k > 1 || (k==1 && buf[0] != '\n')) out += space;
	out += buf;
	free(buf);
	free(space);
	return out;
}

double* get_doubleArray(const char* list, int* num, const char* sep){
    std::string temp(list);
    std::vector<std::string> tokens;
    tokenize(temp, tokens, sep);
    (*num) = tokens.size();
    double *output = new double[*num];
    for(unsigned int i=0; i<tokens.size(); i++){
    	int status = parse_double(tokens[i].c_str(), &(output[i]));
    	if(status != 0) error("ERROR when parsing '%s': It is expected to be a number!\n", tokens[i].c_str());
    }
    return output;
}

int toLinearIndex(const int *dim, const int *dimSize, int ndim){
	int index = dim[0];
	for(int i=1; i<ndim; i++){
		index = (index*dimSize[i])+dim[i];
	}
	return index;
}

void fillInDimIndex(int index, int *dim, const int *dimSize, int ndim){
	int temp = index;
	for(int i=ndim-1; i>=0; i--){
		dim[i] = temp % dimSize[i];
		temp /= dimSize[i];
	}
}
