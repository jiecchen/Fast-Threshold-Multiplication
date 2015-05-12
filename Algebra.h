#ifndef __ALGEBRA_H__
#define __ALGEBRA_H__
#include <cstdlib>
#include <iostream>

typedef struct {
  int row, col, val;
} Element;

// class Matrix_CSC {
// public:
//   // arr has to be sorted
//   Matrix_CSC(Element *arr, int _nnz, int _m, int _n);
//   Matrix_CSC operator * (const Matrix_CSC &m);
//   friend std::ostream& operator << (std::ostream &os, const Matrix_CSC &m);  
//   ~Matrix_CSC() {
//     delete[] values;
//     delete[] row_ind;
//     delete[] col_ptr;
//   }

// private:
//   int m, n; // m by n matrix
//   int nnz; // # of non-zero entries
//   int* values;
//   int* row_ind;
//   int* col_ptr;
// };



class Matrix_CSR {
public:
  Matrix_CSR(Element *arr, int _nnz, int _m,  int _n);
  friend std::ostream& operator << (std::ostream &os, const Matrix_CSR &mat);  

  ~Matrix_CSR() {
    delete[] values;
    delete[] col_ind;
    delete[] row_ptr;
  }

private:
  int m, n;
  int nnz;
  int *values;
  int *col_ind;
  int *row_ptr;
};

#endif

















