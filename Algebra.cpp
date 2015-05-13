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


CMatrix_CSR::CMatrix_CSR(Element *arr,  int _nnz, int _m, int _n) {
  this->values = new int[_nnz];
  this->col_ind = new int[_nnz];
  this->row_ptr = new int[_m + 1];
  std::fill(row_ptr, row_ptr + _m, -1);
  row_ptr[_m] = _nnz;

  this->n = _n;
  this->m = _m;
  this->nnz = _nnz;

  bool marked[_m];
  std::fill(marked, marked + _m, false);
  
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


  // // debug
  // std::cerr << "m = " << m << std::endl;

  // for (int i = 0; i < nnz; ++i)
  //   std::cerr << values[i] << " ";
  // std::cerr << std::endl;
  // for (int i = 0; i < nnz; ++i)
  //   std::cerr << col_ind[i] << " ";
  // std::cerr << std::endl;
  // for (int i = 0; i < m; ++i)
  //   std::cerr << row_ptr[i] << " ";
  // std::cerr << std::endl;
}


CVector CMatrix_CSR::operator*(const CVector &vec) {
  if ((int) vec.size() != this->n) 
    throw "Matrix and vector dimensions are not compatible!";
  
  CVector result(this->m);

  int t = 0;
  int row = 1;
  int pre_row = 0;
  while (t < nnz) {
    while (row_ptr[row] < 0) ++row;
    while (t < row_ptr[row]) {
      result[pre_row] += values[t] * vec[col_ind[t]];
      ++t;
    }
    pre_row = row;
    ++row;
  }
      
  

  return result;
}



std::ostream& operator << (std::ostream &os,  const CMatrix_CSR &mat) {
  CMatrix_COO coo = mat.toCOO();
  int t = 0;
  for (int i = 0; i < mat.m; ++i) {
    for (int j = 0; j < mat.n; ++j) 
      if (coo[t].row == i && coo[t].col == j)
	os << std::setw(5) << coo[t++].val;
      else
	os << std::setw(5) << 0;
    if (i + 1 != mat.m)
      os << std::endl;   
  }
  return os;
}


CMatrix_COO CMatrix_CSR::toCOO() const {
  CMatrix_COO coo;
  int t = 0;
  int row = 1;
  int pre_row = 0;
  while (t < nnz) {
    while (row_ptr[row] < 0) ++row;
    while (t < row_ptr[row]) {
      Element e(pre_row, col_ind[t], values[t]);
      coo.push_back(e);
      ++t;
    }
    pre_row = row;
    ++row;
  }
  return coo;
}


std::ostream& operator << (std::ostream & os, CVector &vec) {
  for (CVector_iter it = vec.begin(); it != vec.end(); ++it)
    os << std::setw(5) << *it;
  return os;
}

















