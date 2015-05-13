#include "Algebra.h"
#include <iomanip>
#include <algorithm>


CMatrix_CSR::CMatrix_CSR(Element *arr,  int _nnz, int _m, int _n) {
  this->col_val = new VectorElement[_nnz];
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
    
    col_val[i].val = arr->val;
    col_val[i].ind = arr->col;
    

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
      result[pre_row] += col_val[t].val * vec[col_val[t].ind];
      ++t;
    }
    pre_row = row;
    ++row;
  }
  return result;
}

CVector CMatrix_CSR::sumRows(Index_iter ibegin, Index_iter iend) {
  //TODO: check if ind legal
  CVector res(this->n);
  for (auto it = ibegin; it != iend; ++it) {
    if (row_ptr[*it] < 0)  continue;
    int i = *it + 1;
    while (row_ptr[i] < 0) ++i;
    int t = row_ptr[*it];
    while (t < row_ptr[i]) {
      res[col_val[t].ind] += col_val[t].val; // values[t];
      ++t;
    }
  }
  return res;
}

CVector CMatrix_CSR::sumRows(int *ibegin, int *iend) {
  //TODO: check if ind legal
  CVector res(this->n);
  for (auto it = ibegin; it != iend; ++it) {
    if (row_ptr[*it] < 0)  continue;
    int i = *it + 1;
    while (row_ptr[i] < 0) ++i;
    int t = row_ptr[*it];
    while (t < row_ptr[i]) {
      res[col_val[t].ind] += col_val[t].val;
      ++t;
    }
  }
  return res;
}






std::ostream& operator << (std::ostream &os,  const CMatrix_CSR &mat) {
  CMatrix_COO coo = mat.toCOO();
  int t = 0;
  for (int i = 0; i < mat.m; ++i) {
    for (int j = 0; j < mat.n; ++j) 
      if (coo[t].row == i && coo[t].col == j)
	os << std::setw(9) << coo[t++].val;
      else
	os << std::setw(9) << 0;
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
      Element e(pre_row, col_val[t].ind, col_val[t].val);
      coo.push_back(e);
      ++t;
    }
    pre_row = row;
    ++row;
  }
  return coo;
}






void CMatrix_COO::sortByRowColumn() {
  std::sort(data.begin(), data.end(), 
	    [](const Element &a, const Element &b) { return a.row < b.row || (a.row == b.row && a.col < b.col); }
	    );
}

void CMatrix_COO::sortByColumnRow() {
  std::sort(data.begin(), data.end(),
	    [](const Element &a, const Element &b) { return a.col < b.col || (a.col == b.col && a.row < b.row); }
	    );
}




std::ostream& operator << (std::ostream & os, CVector &vec) {
  for (auto it = vec.begin(); it != vec.end(); ++it)
    os << std::setw(9) << *it;
  return os;
}


std::ostream& operator << (std::ostream & os, CMatrix_COO &coo) {
  for (auto it = coo.data.begin(); it != coo.data.end(); ++it)
    os << it->row << ", " << it->col << ", " << it->val << std::endl;
  return os;
}




CMatrix_COO thresh_mult_naive(CMatrix_CSR &A, CMatrix_COO &B, double thresh) {
  CMatrix_COO C = B;
  //TODO: check
  CMatrix_COO coo(A.get_m(), B.get_n());

  C.sortByColumnRow();

  int t = 0;
  int col = C[t].col;
  while (t < C.size()) {
    //TODO: will it be initialized as 0?
    CVector vec(A.get_m());
    std::fill(vec.begin(), vec.end(), 0);
    // construct a CVector
    while (C[t].col == col) {
      vec[C[t].row] = C[t].val;
      ++t;
    }
    std::cout << "col = " << col <<" vec: " << vec << std::endl;
    // now I have vec
    CVector res = A * vec;
    std::cout << "col = " << col << " res: " << res << std::endl;
    for (unsigned int i = 0; i < res.size(); ++i)
      if (res[i] > thresh) 
	coo.push_back(Element(i, col, res[i]));

    if (t < C.size())
      col = C[t].col;
  }

  return coo;
}



















