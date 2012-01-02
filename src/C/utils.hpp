/*
	Copyright (c) 2011, Yahoo! Inc.  All rights reserved.
	Copyrights licensed under the New BSD License. See the accompanying LICENSE file for terms.

    Author: Bee-Chung Chen
 */

/*
 * utils.hpp
 *
 */

#ifndef UTILS_HPP_
#define UTILS_HPP_

// Add the following line to the .cpp file that includes this file to
// disable assert()'s  (before #include "utils.hpp")
//   #define NDEBUG
//

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <R.h>

extern "C"{
	#include "util.h"
}


class Matrix_ColumnMajor{
private:
	double* data;
	int _nrow;
	int _ncol;
	bool resizable;
	int size; // number of numbers allocated by this object

public:

	Matrix_ColumnMajor(void){
		_nrow = 0; _ncol = 0; resizable = true;
		size = 0; data = NULL;
	}
	// This constructor allocate memory space
	// The content of the newly created matrix is undefined
	Matrix_ColumnMajor(const int nrow, const int ncol){
		_nrow = nrow; _ncol = ncol; resizable = true;
		size = nrow*ncol;
		assert(size >= 0);
		if(size > 0) data = new double[size];
		else         data = NULL;
	}
	~Matrix_ColumnMajor(void){
		assert(size >= 0);
		assert(size == 0 || resizable);
		if(size > 0) delete[] data;
	}

	// Wrap the input data (no data copy)
	// NOTE: The resulting matrix is not resizable
	inline void wrap(double* data, const int nrow, const int ncol){
		assert(size >= 0);
		assert(size == 0 || resizable);
		if(size > 0) delete[] this->data;
		this->data = data; _nrow = nrow; _ncol = ncol;
		resizable = false;
		size = 0;
	}
	// The content of the newly resized matrix is undefined
	inline void resize(const int nrow, const int ncol){
		if(nrow == _nrow && ncol == _ncol) return;
		if(!resizable) STOP("Trying to resize a non-resizable matrix");
		if(nrow*ncol > size){
			if(size > 0) delete[] data;
			else         assert(data == NULL);
			size = nrow*ncol;
			assert(size > 0);
			data = new double[size];
		}
		_nrow = nrow; _ncol = ncol;
	}
	inline void setToZero(void){
		for(int k=0; k<_nrow*_ncol; k++) data[k] = 0;
	}

	inline double* getData(void) const {
		return data;
	}

	// A(i,j) gets the value in row i and column j of matrix A
	// NOTE: i and j start from 0 (not 1)
	inline double& operator()(const int i, const int j) const {
		assert(i >= 0 && i < _nrow && j >=0 && j< _ncol);
		return data[C_MAT(i,j,_nrow)];
	}

	// copy the content of the other matrix to this one
	Matrix_ColumnMajor& operator= (const Matrix_ColumnMajor& other){
		if(resizable){
			resize(other._nrow, other._ncol);
		}else{
			if(_nrow != other._nrow) STOP2("Matrix assignment #rows mismatch: %d vs. %d", _nrow, other._nrow);
			if(_ncol != other._ncol) STOP2("Matrix assignment #columns mismatch: %d vs. %d", _ncol, other._ncol);
		}
		memcpy(data, other.data, _nrow*_ncol*sizeof(double));
		return *this;
	}

	inline int nrow(void) const {return _nrow;}
	inline int ncol(void) const {return _ncol;}
	inline int length(void) const {return _nrow*_ncol;}
	inline bool isResizable() const {return resizable;}

	// invert the symmetric matrix (option=1 -> Cholesky)
	// Note: check_sym=0 -> do not check symmetry; otherwise, check for symmetry
	inline void sym_invert(int check_sym=0, int option=1){
		if(_nrow != _ncol) STOP2("Trying to invert a non-square matrix: %d x %d", _nrow, _ncol);
		sym_inv_byCholesky(data, &_nrow, &check_sym);
	}

	// Make a conceptually symmetric matrix exactly symmetric
	inline void symmetrize(int debug=0, double abs_tol=1e-10, double rel_tol=1e-8){
		if(_nrow != _ncol) STOP2("Trying to symmetrize a non-square matrix: %d x %d", _nrow, _ncol);
		double max_value = 0;
		if(debug){
			for(int m=0; m<_nrow*_ncol; m++) if(fabs(data[m]) > max_value) max_value = fabs(data[m]);
		}
		for(int i=0; i<_nrow; i++) for(int j=0; j<i; j++){
			double A_ij = data[C_MAT(i,j,_nrow)], A_ji = data[C_MAT(j,i,_nrow)];
			if(debug && fabs(A_ij - A_ji) > abs_tol && fabs(A_ij - A_ji) > rel_tol*max_value){
				char space[] = "   ";
				if(debug >= 10) print(space);
				error("A symmetric matrix is not numerically symmetric: A[%d,%d] = %.12f, A[%d,%d] = %.12f, diff=%e (source file: %s)", i, j, A_ij, j, i, A_ji, A_ij-A_ji, __FILE__);
			}
			double avg = (A_ij + A_ji)/2;
			data[C_MAT(i,j,_nrow)] = avg;
			data[C_MAT(j,i,_nrow)] = avg;
		}
	}

	inline void elementwiseInverse(void){
		for(int m=0; m<_nrow*_ncol; m++) data[m] = 1.0/data[m];
	}

	inline void negate(void){
		for(int m=0; m<_nrow*_ncol; m++) data[m] = -data[m];
	}

	// this = rbind(A, B)
	inline void rbind(const Matrix_ColumnMajor& A, const Matrix_ColumnMajor& B){
		assert(A._ncol == B._ncol);
		assert(data == NULL || (A.data != data && B.data != data));
		resize(A._nrow + B._nrow, A._ncol);
		for(int i=0; i<A._nrow; i++){
			for(int j=0; j<A._ncol; j++)
				data[C_MAT(i,j,_nrow)] = A.data[C_MAT(i,j,A._nrow)];
		}
		for(int i=0; i<B._nrow; i++){
			for(int j=0; j<B._ncol; j++)
				data[C_MAT(i+A._nrow,j,_nrow)] = B.data[C_MAT(i,j,B._nrow)];
		}
	}

	// this = diagonal(d) %*% B
	// dim(d) is either 1x1 or nrow(B) x 1
	inline void product_1stDiagonal(const Matrix_ColumnMajor& d, const Matrix_ColumnMajor& B){
		assert(d._ncol == 1);
		assert(data == NULL || (d.data != data && B.data != data));
		resize(B._nrow, B._ncol);
		if(d.length() == 1){
			for(int m=0; m<_nrow*_ncol; m++) data[m] = d.data[0] * B.data[m];
		}else if(d.length() == B._nrow){
			for(int i=0; i<_nrow; i++){
				for(int j=0; j<_ncol; j++)
					data[C_MAT(i,j,_nrow)] = d.data[i] * B.data[C_MAT(i,j,B._nrow)];
			}
		}else STOP2("length(d) = %d and nrow(B) = %d", d.length(), B._nrow);
	}
	// this = B %*% diagonal(d)
	// dim(d) is either 1x1 or ncol(B) x 1
	inline void product_2ndDiagonal(const Matrix_ColumnMajor& B, const Matrix_ColumnMajor& d){
		assert(d._ncol == 1);
		assert(data == NULL || (d.data != data && B.data != data));
		resize(B._nrow, B._ncol);
		if(d.length() == 1){
			for(int m=0; m<_nrow*_ncol; m++) data[m] = d.data[0] * B.data[m];
		}else if(d.length() == B._ncol){
			for(int i=0; i<_nrow; i++){
				for(int j=0; j<_ncol; j++)
					data[C_MAT(i,j,_nrow)] = d.data[j] * B.data[C_MAT(i,j,B._nrow)];
			}
		}else STOP2("length(d) = %d and ncol(B) = %d", d.length(), B._ncol);
	}
	// this = A %*% B
	// NOTE: This matrix cannot be either A or B
	inline void product(const Matrix_ColumnMajor& A, const Matrix_ColumnMajor& B, bool allowScalar=false){
		if(allowScalar && A.length() == 1){
			product_1stDiagonal(A, B); return;
		}
		if(allowScalar && B.length() == 1){
			product_2ndDiagonal(A, B); return;
		}
		assert(A._ncol == B._nrow);
		assert(data == NULL || (A.data != data && B.data != data));
		resize(A._nrow, B._ncol);
		for(int i=0; i<_nrow; i++){
			for(int j=0; j<_ncol; j++){
				data[C_MAT(i,j,_nrow)] = 0;
				for(int k=0; k<A._ncol; k++)
					data[C_MAT(i,j,_nrow)] += A.data[C_MAT(i,k,A._nrow)] * B.data[C_MAT(k,j,B._nrow)];
			}
		}
	}


	// this = t(A)
	inline void transpose(const Matrix_ColumnMajor& A){
		assert(data == NULL || A.data != data);
		resize(A._ncol, A._nrow);
		for(int i=0; i<_nrow; i++){
			for(int j=0; j<_ncol; j++)
				data[C_MAT(i,j,_nrow)] = A.data[C_MAT(j,i,A._nrow)];
		}
	}

	// this = this + diagonal(v)
	// length = length(v);  length = 1 to add the same number to all diagonal entries
	inline void addDiagonal(const double* v, int length){
		assert(_nrow == _ncol);
		assert(length == 1 || length == _nrow);
		if(length == 1){
			for(int i=0; i<_nrow; i++) data[C_MAT(i,i,_nrow)] += v[0];
		}else{
			for(int i=0; i<_nrow; i++) data[C_MAT(i,i,_nrow)] += v[i];
		}
	}
	// v must be a column vector
	inline void addDiagonal(const Matrix_ColumnMajor& v){
		assert(_nrow == _ncol);
		assert(v._ncol == 1 && (v._nrow == _nrow || v._nrow == 1));
		if(v.length() == 1){
			for(int i=0; i<_nrow; i++) data[C_MAT(i,i,_nrow)] += v.data[0];
		}else{
			for(int i=0; i<_nrow; i++) data[C_MAT(i,i,_nrow)] += v.data[i];
		}
	}
	// this = this + v I
	inline void addDiagonal(const double v){
		assert(_nrow == _ncol);
		for(int i=0; i<_nrow; i++) data[C_MAT(i,i,_nrow)] += v;
	}

	// v must be a column vector
	inline void subtractDiagonal(const Matrix_ColumnMajor& v){
		assert(_nrow == _ncol);
		assert(v._ncol == 1 && (v._nrow == _nrow || v._nrow == 1));
		if(v.length() == 1){
			for(int i=0; i<_nrow; i++) data[C_MAT(i,i,_nrow)] -= v.data[0];
		}else{
			for(int i=0; i<_nrow; i++) data[C_MAT(i,i,_nrow)] -= v.data[i];
		}
	}
	// this = this + A
	inline void add(const Matrix_ColumnMajor& A, bool allowScalar=false){
		if(allowScalar && A.length() == 1){
			addDiagonal(A); return;
		}
		assert(_nrow == A._nrow && _ncol == A._ncol);
		for(int m=0; m<_nrow*_ncol; m++) data[m] += A.data[m];
	}
	// this = this - A
	inline void subtract(const Matrix_ColumnMajor& A, bool allowScalar=false){
		if(allowScalar && A.length() == 1){
			subtractDiagonal(A); return;
		}
		assert(_nrow == A._nrow && _ncol == A._ncol);
		for(int m=0; m<_nrow*_ncol; m++) data[m] -= A.data[m];
	}

	// this = v * this (v is a scalar)
	inline void scale(double v){
		for(int k=0; k<_nrow*_ncol; k++) data[k] *= v;
	}

	/**
	 * print the matrix for debugging
	 */
	void print(char* prefix, FILE* fp=stdout) const {
		for(int i=0; i<_nrow; i++){
			fprintf(fp, "%s", prefix);
			for(int j=0; j<_ncol; j++)
				fprintf(fp," %11.7f", data[C_MAT(i,j,_nrow)]);
			fprintf(fp, "\n");
		}
	}
};


class EigenMatrix {
private:
	double *value;  // eigen value:  n x 1
	double *vector; // eigen vector: n x n
	double *workspace;
	int n;
	int workspace_size;
public:
	EigenMatrix(const int n /* for n x n matrix */){
		this->n = n;
		value  = new double[n];
		vector = new double[n*n];
	    // prepare work space
		workspace_size = sym_eigen2_workspace_size(&n);
		workspace = new double[workspace_size];
	}
	~EigenMatrix(void){
		delete[] value;  delete[] vector;  delete[] workspace;
	}
	// this = eigen_decomposition(A)
	inline void decompose(const Matrix_ColumnMajor& A, int debug=0){
		if(A.nrow() != n) STOP2("nrow(A) = %d, n = %d", A.nrow(), n);
		if(A.ncol() != n) STOP2("ncol(A) = %d, n = %d", A.ncol(), n);
		sym_eigen2(A.getData(), &n, value, vector, workspace, &workspace_size, &debug);
	}
	inline void invert(void){
		for(int i=0; i<n; i++) value[i] = 1/value[i];
	}
	inline int nrow(void) const { return n; }
	inline int ncol(void) const { return n; }
	inline double* getValue(void)  const { return value; }
	inline double* getVector(void) const { return vector; }

	// copy the content of an EigenMatrix to matrix A
	void outputTo(Matrix_ColumnMajor& A){
		if(A.isResizable()){
			A.resize(n, n);
		}else{
			if(n != A.nrow()) STOP2("Matrix assignment #rows mismatch: %d vs. %d", n, A.nrow());
			if(n != A.ncol()) STOP2("Matrix assignment #columns mismatch: %d vs. %d", n, A.ncol());
		}
		double *data = A.getData();
		for(int i=0; i<n*n; i++) data[i] = 0;
		for(int k=0; k<n; k++){
			for(int i=0; i<n; i++){
				for(int j=0; j<n; j++)
					data[C_MAT(i,j,n)] += value[k] * vector[C_MAT(i,k,n)] * vector[C_MAT(j,k,n)];
			}
		}
	}

};


/**
 * Conceptually: list(M1, M2, ...)
 * Input:
 *  data = a vector that concat all matrices
 *  dim  = a vector of the following form
 *         #matrices, start(M1), nrow(M1), ncol(M1), start(M2), nrow(M2), ncol(M2), ...
 * Note:
 *  (1) For sanity check, a random number is appended to the end of both data and dim.
 *      The two numbers should match.
 *  (2) Input starting indices are R indices (starting from 1)
 */
class ListOfMatrices {
private:
	double* data;
	int nMatrices;
	Matrix_ColumnMajor* matrix; // matrix[k] is the kth matrix
public:
	ListOfMatrices(double* data, const int* dim){
		this->data = data;
		nMatrices = dim[0];
		matrix = new Matrix_ColumnMajor[nMatrices];
		double sig1 = dim[1+3*nMatrices];
		int size = 0;
		for(int k=0; k<nMatrices; k++){
			int start = dim[k*3+1]-1;
			int nrow  = dim[k*3+2];
			int ncol  = dim[k*3+3];
			matrix[k].wrap(&(data[start]), nrow, ncol);
			if(start != size) STOP("Corrupted input dim");
			size += nrow*ncol;
		}
		double sig2 = data[size];
		if(sig1 != sig2) STOP("Sanity check failed");
	}
	~ListOfMatrices(void){
		delete[] matrix;
	}
	/**
	 * length of the list (number of matrices)
	 */
	inline int length(void) const { return nMatrices; }
	/**
	 * get matrix
	 */
	inline Matrix_ColumnMajor& get(int k) const {
		return matrix[k];
	}
	/**
	 * print the list of matrices for debugging
	 */
	void print(FILE* fp=stdout) const {
		for(int k=0; k<nMatrices; k++){
			fprintf(fp,"[[%d]]\n", k);
			char space[] = "  ";
			matrix[k].print(space, fp);
		}
	}
};

/**
 * Reverse index
 *   Consider index value i
 *      i = indexArray[ revIndex[i][0:(num[i]-1)] ] - 1
 *      NOTE: ... - 1 because of R to C conversion
 *   E.g., consider user i
 *      y[ revIndex[i][0:(num[i]-1)] ]
 *      are the observations of user i
 */
class ReverseIndex {
private:
    int* index_space;
public:
    int** revIndex;
    int* num;
    int size;

    /**
     * The indices in the input indexArray start from 1 (not 0)
     */
    ReverseIndex(const int* indexArray, const int arrayLength, const int nIndexValues){

        revIndex = new int*[nIndexValues];
        num      = new int[nIndexValues];
        index_space = new int[arrayLength];
        size = nIndexValues;

        for(int i=0; i<nIndexValues; i++) num[i] = 0;

        for(int k=0; k<arrayLength; k++){
            int i = indexArray[k] - 1;
            CHK_C_INDEX(i, nIndexValues);
            num[i]++;
        }

        revIndex[0] = index_space;
        for(int i=1; i<nIndexValues; i++){
            revIndex[i] = &(revIndex[i-1][num[i-1]]);
        }

        int *current = new int[nIndexValues];
        for(int i=0; i<nIndexValues; i++) current[i] = 0;

        for(int k=0; k<arrayLength; k++){
            int i = indexArray[k] - 1;
            revIndex[i][current[i]] = k;
            current[i]++;
        }
        delete[] current;
    }
    ~ReverseIndex(void){
        delete[] revIndex;
        delete[] num;
        delete[] index_space;
    }
};

#endif /* UTILS_HPP_ */
