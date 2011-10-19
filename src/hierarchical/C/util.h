#ifndef UTIL_H
#define UTIL_H

#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <string.h>
#include <algorithm>
#include <vector>
#include <ctime>
#include <ctype.h>
#include <float.h>
#include <errno.h>
#include <limits.h>
#include <sys/stat.h>
#include <sys/types.h>

#define CHK_C_INDEX(index, num) if(index < 0 || index >= num) error("index out of bound: index=%d, bound=%d (file: %s, line: %d)", index, num, __FILE__, __LINE__)
#define CHK_C_INDEX_2D(row_i, col_j, nrow, ncol) if(row_i < 0 || row_i >= nrow || col_j < 0 || col_j >= ncol) error("index out of bound: i=%d, j=%d nrow=%d ncol=%d (file: %s, line: %d)", row_i, col_j, nrow, ncol, __FILE__, __LINE__)

#define STOP_HERE(msg) error("Error at %s on line %d: %s\n", __FILE__, __LINE__, msg)
#define STOP2(msg,x) error2(__FILE__, __LINE__, msg, x)
#define STOP3(msg,x1,x2) error2(__FILE__, __LINE__, msg, x1, x2)
#define STOP4(msg,x1,x2,x3) error2(__FILE__, __LINE__, msg, x1, x2, x3)

#define CHK_SAME_NUMBER_(x, y, rel_tol, abs_tol) if((x != 0 || y != 0) && (fabs(x-y) / fmax(fabs(x), fabs(y)) > rel_tol) && fabs(x-y) > abs_tol) error("Error: The two number should be the same: %16g vs %16g (file: %s, line: %d)\n", x, y, __FILE__, __LINE__)
#define CHK_SAME_NUMBER(x, y) CHK_SAME_NUMBER_(x, y, 1e-8, 1e-12)
#define WARN_SAME_NUMBER(x, y, prec) if((x != 0 || y != 0) && (fabs(x-y) / fmax(fabs(x), fabs(y)) > prec)) printf("WARN: The two number should be the same: %16g vs %16g (file: %s, line: %d)\n", x, y, __FILE__, __LINE__)

#define C_MAT(i,j,nrow) (((j)*(nrow))+(i))
#define C_3DA(i,j,k,nrow,ncol) ((((k)*(ncol))+(j))*(nrow) + (i))
#define C_4DA(i,j,k,m,dim1,dim2,dim3) (((((m)*(dim3)+(k))*(dim2))+(j))*(dim1) + (i))

#define Malloc(type,n) (type *)malloc((n)*sizeof(type))
#define SQR(x) ((x)*(x))
#define MAX(x,y) ((x) > (y) ? (x) : (y))
#define MIN(x,y) ((x) < (y) ? (x) : (y))
#define LOGIT(x) (log((x)/(1-(x))))
#define INV_LOGIT(x) (1/(1+exp(-(x))))

// Index conversion from Matrix (m x n) to Vector ((m*n) x 1)
#define INDEX_M2V(i,j,nrow) (((j)*(nrow))+(i))

#define FORCE_RANGE(value, min, max) (value > max ? max : (value < min ? min : value))

void error(const char *fmt, ...);

// Consider index value i
//   i = indexArray[ revIndex[i][0:(num[i]-1)] ]
// E.g., consider user i
//      y[ revIndex[i][0:(num[i]-1)] ]
//      are the observations of user i
class CReverseIndex {
private:
    int* index_space;
public:
    int** revIndex;
    int* num;
    int size;
    
    CReverseIndex(const int* indexArray, const int arrayLength, const int nIndexValues);
    ~CReverseIndex();
};

class MyTimer {
private:
    time_t  t_begin;
    clock_t c_begin;
public:
    MyTimer(void){ reset(); }
    void reset(void){ t_begin = time(NULL); c_begin = clock(); }
    int wallTime(void){ return (int)(time(NULL) - t_begin); }
    double cpuTime(void){ return ((double)(clock() - c_begin)) / CLOCKS_PER_SEC; }
};

std::string indent(const char* input, int n);
void info(FILE* fp, const int verboseLevel, const char* fmt, ...);

void error2(const char *filename, int lineno, const char *fmt, ...);

void tokenize(const std::string& str, std::vector<std::string>& tokens, const std::string& delimiters = " ");

void print_matrix(double **matrix, const int nrow, const int ncol, FILE* fp=stdout, const char* sep="\t", const char* format="%.16g");
void print_row_vector(double *vector, const int num, FILE* fp=stdout, const char* sep="\t", const char* format="%.16g");
void print_col_vector(double *vector, const int num, FILE* fp=stdout, const char* sep="\t", const char* format="%.16g");

void print_matrix(int **matrix, const int nrow, const int ncol, FILE* fp=stdout, const char* sep="\t", const char* format="%d");
void print_row_vector(int *vector, const int num, FILE* fp=stdout, const char* sep="\t", const char* format="%d");
void print_col_vector(int *vector, const int num, FILE* fp=stdout, const char* sep="\t", const char* format="%d");

inline void swap(int& a, int& b){
	int c = b; b = a; a = c;
}
inline void swap(double& a, double& b){
	double c = b; b = a; a = c;
}

inline void save_matrix(double **matrix, const int nrow, const int ncol, const char* file, const char* sep="\t", const char* format="%.16g"){
	FILE *fp = fopen(file, "w");
	if(fp == NULL) error("ERROR: Cannot open file %s\n", file);
	print_matrix(matrix, nrow, ncol, fp, sep, format);
	fclose(fp);
}
inline void save_col_vector(double *vector, const int num, const char* file, const char* sep="\t", const char* format="%.16g"){
	FILE *fp = fopen(file, "w");
	if(fp == NULL) error("ERROR: Cannot open file %s\n", file);
	print_col_vector(vector, num, fp, sep, format);
	fclose(fp);
}

inline int parse_double(const char *str, double* ans){
	char *endptr;
    errno = 0;
    (*ans) = strtod(str,&endptr);
    if(errno != 0) return errno;
    if(endptr == str || (*endptr != '\0' && !isspace(*endptr))) return -1;
    return 0;
}

inline int parse_int(const char *str, int* ans, int base=10){
	char *endptr;
    errno = 0;
    (*ans) = (int)strtol(str,&endptr, base);
    if(errno != 0) return errno;
    if(endptr == str || (*endptr != '\0' && !isspace(*endptr))) return -1;
    return 0;
}

inline bool file_exist(const char *file){
    struct stat statbuf;
    if(stat(file, &statbuf) == 0) return true;
	return false;
}

inline void check_file_overwrite(const char *file, bool allow_overwrite){
    struct stat statbuf;
    if(stat(file, &statbuf) == 0){
        if(allow_overwrite) fprintf(stderr, "\nWARNING: File/directory %s exists and will be overwritten!!\n\n", file);
        else error("\nERROR: %s aready exists and you chose not to overwrite it!!\n\n", file);
    }
}

inline void normalizeToSumUpToOne(double *vector, const int length){
    double sum = 0;
    for(int k=0; k<length; k++) sum += vector[k];
    for(int k=0; k<length; k++) vector[k] /= sum;
}

inline int which_max(double *array, int num){
	if(array == NULL || num == 0) error("ERROR: The input array is empty (source: %s@%d)\n", __FILE__, __LINE__);
	int best=0; double max = array[0];
	for(int i=1; i<num ;i++) if(array[i] > max){
		max = array[i];
		best = i;
	}
	return best;
}
inline int which_min(double *array, int num){
	if(array == NULL || num == 0) error("ERROR: The input array is empty (source: %s@%d)\n", __FILE__, __LINE__);
	int best=0; double min = array[0];
	for(int i=1; i<num ;i++) if(array[i] < min){
		min = array[i];
		best = i;
	}
	return best;
}

double* get_doubleArray(const char* list, int* num, const char* sep=",");


/** Convert dimenional index to linear index.
 *  dim[i]: the value of the i-th dimension
 *  dimSize[i]: the number of values in the i-th dimension
 */
int toLinearIndex(const int *dim, const int *dimSize, int ndim);

/** Convert linear index to dimensional index
 */
void fillInDimIndex(int index, int *dim, const int *dimSize, int ndim);

/**
 * A wrapper for R array
 */
template <typename ValueType>
class R_MultiDimArray {
private:
	ValueType *data; // data[0, ..., length-1]
	int *dimSize;    // dimSize[i] is the size of the ith dimension
	int nDim;
public:
	int length;      // length = dimSize[0] x ... x dimSize[nDim-1]

	R_MultiDimArray(ValueType *array, const int *dimSize, int nDim){
		if(nDim <= 0) STOP2("nDim = %d", nDim);
		data=array; this->nDim = nDim;
		this->dimSize = new int[nDim]; for(int i=0; i<nDim; i++) this->dimSize[i] = dimSize[i];
		length = 1; for(int i=0; i<nDim; i++) length *= dimSize[i];
	}
	~R_MultiDimArray(void){ delete[] dimSize; }

	inline unsigned int index(const int *entry) const {
		int index = 0;
		for(int i=nDim-1; i>=0; i--){
			if(entry[i] < 0 || entry[i] >= dimSize[i]) STOP4("Index out of bound: i=%d, entry[i]=%d, dimSize[i]=%d", i, entry[i], dimSize[i]);
			index = (index*dimSize[i])+entry[i];
		}
		return index;
	}

	inline ValueType get(const int *entry) const {
		return data[index(entry)];
	}
	inline void set(const int *entry, ValueType value){
		data[index(entry)] = value;
	}
	inline ValueType getBy1DIndex(int i) const {
		return data[i];
	}
	inline void setBy1DIndex(int i, ValueType value){
		data[i] = value;
	}
};

template <typename ValueType>
class R_2DArray {
private:
	ValueType *data;
	int nrow, ncol;
public:
	R_2DArray(ValueType *array, int nrow, int ncol){
		data = array; this->nrow = nrow; this->ncol = ncol;
	}
	inline ValueType get(int i, int j) const {
		if(i < 0 || i >= nrow || j < 0 || j >= ncol) STOP3("Index out of bound: i=%d, j=%d", i, j);
		return data[C_MAT(i,j,nrow)];
	}
	inline void set(int i, int j, ValueType value){
		if(i < 0 || i >= nrow || j < 0 || j >= ncol) STOP3("Index out of bound: i=%d, j=%d", i, j);
		data[C_MAT(i,j,nrow)] = value;
	}
	inline int nRow() const { return nrow; }
	inline int nCol() const { return ncol; }
	inline void getRow(int i, ValueType row[]) const {
		for(int j=0; j<ncol; j++) row[j] = data[C_MAT(i,j,nrow)];
	}
	inline void getCol(int j, ValueType col[]) const {
		for(int i=0; i<nrow; i++) col[i] = data[C_MAT(i,j,nrow)];
	}
};

#endif
