#include "Algebra.h"
#include <iomanip>
#include <algorithm>
#include <ctime>

CMatrix_CSR::CMatrix_CSR(Element *arr_s, Element *arr_e , int _m, int _n) {
  int _nnz = arr_e - arr_s;
  this->col_val = new VectorElement[_nnz];
  this->row_ptr = new int[_m + 1];
  std::fill(row_ptr, row_ptr + _m, -1);
  row_ptr[_m] = _nnz;

  this->n = _n;
  this->m = _m;
  this->nnz = _nnz;

  bool marked[_m];
  std::fill(marked, marked + _m, false);
  
  Element *arr = arr_s;
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
}







CVector operator*(const CMatrix_CSR &mat,  CVector &vec) {
  if ((int) vec.size() != mat.n) 
    throw "Matrix and vector dimensions are not compatible!";
  
  CVector result(mat.m);

  int t = 0;
  int row = 1;
  int pre_row = 0;
  while (t < mat.nnz) {
    while (mat.row_ptr[row] < 0) ++row;
    while (t < mat.row_ptr[row]) {
      result[pre_row] += mat.col_val[t].val * vec[mat.col_val[t].ind];
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

std::ostream& operator << (std::ostream & os, CDenseMatrix &mat) {
  for (int i = 0; i < mat.get_m(); ++i) {
    int *data = mat.get_row(i);
    for (int j = 0; j < mat.get_n(); ++j) 
      os << std::setw(9) << data[j];
    os << std::endl;
  }
  return os;
}






// complexity \prop  nnz(A) * B.n

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
    //    std::cout << "col = " << col <<" vec: " << vec << std::endl;
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



// complexity \prop h * nnz(B)
CDenseMatrix operator *(CMatrix_COO &A, CMatrix_COO &B) {
  int *data = new int[A.get_m() * B.get_n()]();
  int v[A.get_n()];
  A.sortByRowColumn();
  //  B.sortByColumnRow();
  int t = 0;
  int row = A[t].row;
  while (t < A.size()) {
    // construct the row of A
    std::fill(v, v + A.get_n(), 0);
    v[A[t].col] = A[t].val;
    ++t;
    while (t < A.size() && A[t-1].row == A[t].row) {
      v[A[t].col] = A[t].val;
      ++t;
    }
    CVector cv = CVector(v, v + A.get_n());
    //    std::cout << "row = " << row << " v = " << cv  << std::endl;

    // calculate v' * B
    for (int i = 0; i < B.size(); i++)
      data[row * B.get_n() + B[i].col] += B[i].val * v[B[i].row];
    

    if (t < A.size())
      row = A[t].row;
        
  }
  return CDenseMatrix(data, A.get_m(), B.get_n());
}


// same as the operator *, but return a int*
int* coo_mult(CMatrix_COO &A, CMatrix_COO &B) {
  int *data = new int[A.get_m() * B.get_n()]();
  int v[A.get_n()];
  A.sortByRowColumn();
  //  B.sortByColumnRow();
  int t = 0;
  int row = A[t].row;
  while (t < A.size()) {
    // construct the row of A
    std::fill(v, v + A.get_n(), 0);
    v[A[t].col] = A[t].val;
    ++t;
    while (t < A.size() && A[t-1].row == A[t].row) {
      v[A[t].col] = A[t].val;
      ++t;
    }
    CVector cv = CVector(v, v + A.get_n());
    //    std::cout << "row = " << row << " v = " << cv  << std::endl;

    // calculate v' * B
    for (int i = 0; i < B.size(); i++)
      data[row * B.get_n() + B[i].col] += B[i].val * v[B[i].row];
    

    if (t < A.size())
      row = A[t].row;
        
  }
  return data;
}











void createCountMin(CountMinSketch &sk, double eps, int _u, CMatrix_COO &coo) {
  std::srand(time(NULL));
  sk.w = std::ceil(1. / eps);
  sk.u = _u;
  int h = sk.w * sk.u;
  CMatrix_COO* cm =  new CMatrix_COO(h, coo.get_m());
  std::cout << "w = " << w << " u = " << u << " n = " << coo.get_m() << std::endl; 
  for (int i = 0; i < coo.get_m(); ++i) 
    for (int j = 0; j < u; ++j) 
      cm->push_back(Element(std::rand() % w + j * w, i, 1));
  // get the sketch of the coo
  //  std::cout << "I  am here" << std::endl;
  sk.hash = cm;
  sk.cm = coo_mult(cm, coo);
}





CDenseMatrix::CDenseMatrix(int m, int n): m(m), n(n) {
    data = new int[m * n](); 
}

CDenseMatrix::CDenseMatrix(int *_data, int m, int n): m(m), n(n) {
  this->data = _data;
} 





