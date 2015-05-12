#include "Algebra.h"
#include <iomanip>

// Matrix_CSC::Matrix_CSC(Element *arr, int _nnz, int _m,  int _n) {
//   this->values = new int[_nnz];
//   this->row_ind = new int[_nnz];
//   this->col_ptr = new int[_n];
//   this->n = _n;
//   this->m = _m;
//   this->nnz = _nnz;

//   bool marked[_n];
//   std::fill(marked, marked + _n, false);

//   for (int i = 0; i < _nnz; ++i) {
//     int col = arr->col;
//     int row = arr->row;
//     values[i] = arr->val;
//     row_ind[i] = row;
//     if (!marked[col]) {
//       col_ptr[col] = row;
//       marked[col] = true;
//     }
//     arr++;
//   }
// }


// std::ostream& operator << (std::ostream &os, const Matrix_CSC &m) {
//   os << "Hello World!";
//   return os;
// }


Matrix_CSR::Matrix_CSR(Element *arr,  int _nnz, int _m, int _n) {
  this->values = new int[_nnz];
  this->col_ind = new int[_nnz];
  this->row_ptr = new int[_m + 1];
  row_ptr[_m] = _nnz;
  this->n = _n;
  this->m = _m;
  this->nnz = _nnz;

  bool marked[_m];
  std::fill(marked, marked + _n, false);
  
  for (int i = 0; i < _nnz; ++i) {
    int row = arr->row; 
    values[i] = arr->val;
    col_ind[i] = arr->col;
    if (!marked[row]) {
      row_ptr[row] = i;
      marked[row] = true;
    }
    ++arr;
  }
}


std::ostream& operator << (std::ostream &os,  const Matrix_CSR &mat) {
  int t = 0;
  for (int i = 0; i < mat.m; ++i) {
    for (int j = 0; j < mat.n; ++j) {
      if (mat.row_ptr[i + 1] > t && mat.col_ind[t] == j) 
	os << std::setw(5) << mat.values[t++];
      else
	os << std::setw(5) << 0;
    }
    os << std::endl;
  }
  return os;
}










