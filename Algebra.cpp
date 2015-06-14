#include "Algebra.h"
#include <iomanip>
#include <algorithm>
#include <ctime>
#include <map>





///////////////////////////////////////////////////////////////////////////////
////////////////////////// CMatrix_CSC ////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
// has to be sorted by ColumnRow
void CMatrix_CSC::init(int *_val, int *_row, int *_col, 
			 int _nnz,  int _m,  int _n) {
  this->m = _m;
  this->n = _n;
  this->nnz = _nnz;
  this->val = new int[_nnz];
  this->row = new int[_nnz];
  this->col_ptr = new int[_n + 1];
  std::copy(_val, _val + _nnz, val);
  std::copy(_row, _row + _nnz, row);
  std::fill(col_ptr, col_ptr + _n + 1, _nnz);

  for (int i = 0; i < _nnz; i++)
    if (col_ptr[_col[i]] == _nnz)
      col_ptr[_col[i]] = i;

  for (int i = _n - 1; i >= 0; i--) 
    if (col_ptr[i] == _nnz)
      col_ptr[i] = col_ptr[i + 1];
}

CMatrix_CSC::CMatrix_CSC(CVector_iter _val, CVector_iter _row, CVector_iter _col, 
			 int _nnz,  int _m,  int _n): m(_m), n(_n), nnz(_nnz) {
  this->val = new int[_nnz];
  this->row = new int[_nnz];
  this->col_ptr = new int[_n + 1];
  std::copy(_val, _val + _nnz, val);
  std::copy(_row, _row + _nnz, row);
  std::fill(col_ptr, col_ptr + _n + 1, _nnz);

  for (int i = 0; i < _nnz; i++)
    if (col_ptr[_col[i]] == _nnz)
      col_ptr[_col[i]] = i;

  for (int i = _n - 1; i >= 0; i--) 
    if (col_ptr[i] == _nnz)
      col_ptr[i] = col_ptr[i + 1];
}






// // v_s - v_e keeps the row numbers of non-zero entries
// void thresh_mult(SparseVector &result, CMatrix_CSC &csc, int *v_s, int *v_e, double thresh) {
//   std::map<int, int> res;
//   while (v_s != v_e) {
//     for (int i = csc.col_ptr[*v_s]; i < csc.col_ptr[(*v_s) + 1]; ++i) 
//       res[csc.row[i]]++;
//     ++v_s;
//   }
//   for (auto it = res.begin(); it != res.end(); ++it)
//     if (it->second > thresh)
//       result.push_back(VectorElement(it->first, it->second));
// }



//TODO:
//  + bug: when A * B is zero matrix
CMatrix_CSC operator *(const CMatrix_CSC &A, const CMatrix_CSC &B) {
  CVector val, row, col;
  for (int i = 0; i < B.n; ++i) {
    std::map<int, int> res;
    for (int j = B.col_ptr[i]; j < B.col_ptr[i + 1]; ++j) 
      for (int t = A.col_ptr[B.row[j]]; t < A.col_ptr[B.row[j] + 1]; ++t)
	res[A.row[t]] += A.val[t] * B.val[j];
    
    for (auto it = res.begin(); it != res.end(); ++it) {
      val.push_back(it->second);
      row.push_back(it->first);
      col.push_back(i);
    }
  }
  return CMatrix_CSC(val.begin(), row.begin(), col.begin(), val.size(), A.m, B.n);
}

CMatrix_COO toCoo(const CMatrix_CSC &mat) {
  CMatrix_COO coo;
  for (int i = 0; i < mat.n; i++)
    for (int j = mat.col_ptr[i]; j < mat.col_ptr[i + 1]; ++j)
      coo.push_back(Element(mat.row[j], i, mat.val[j]));
  coo.set_mn(mat.m, mat.n);
  return coo;
}

///////////////////////////////////////////////////////////////////////////////
////////////////////////// CMatrix_COO ////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

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


  // Transpose
CMatrix_CSC CMatrix_CSC::T() const {
  CMatrix_COO&& coo = toCoo(*this);
  CMatrix_COO&& cooT = coo.T();
  return CMatrix_CSC(cooT);
}




std::ostream& operator << (std::ostream & os, CMatrix_COO coo) {
  //  int A[coo.get_m()][coo.get_n()];
  //  std::fill(A[0], A[0] + coo.get_m() * coo.get_n(), 0);
  coo.sortByRowColumn();
  int t = 0;
  for (int i = 0; i < coo.get_m(); ++i) {
    for (int j = 0; j < coo.get_n(); ++j)
      if (t < coo.size() && coo[t].row == i && coo[t].col == j) {
	os << std::setw(9) << coo[t].val;
	++t;
      }
      else 
	os << std::setw(9) << 0;
    os << std::endl;
  }

  // for (auto it = coo.data.begin(); it != coo.data.end(); ++it)
  //   os << it->row << ", " << it->col << ", " << it->val << std::endl;
  return os;
}









// // sum rows of an coo matrix
// void sumRows_Coo(int *result, CMatrix_COO &P) {
//   std::fill(result, result + P.get_n(), 0);
//   for (int i = 0; i < P.size(); ++i)
//     result[P[i].col] += P[i].val;
// }

// int inner_prod(int *v, SparseVector &sv) {
//   int res = 0;
//   for (auto it = sv.begin(); it != sv.end(); ++it)
//     res += v[it->ind] * it->val;
//   return res;
// }



std::ostream& operator << (std::ostream &os,  const SparseVector &vec) {
  for (const VectorElement &elem : vec)
    os << elem.ind << ", " << elem.val << std::endl;
  return os;
}


std::ostream& operator << (std::ostream &os,  const VectorElement& e) {
  os << e.ind << ", " << e.val;
  return os;
}



bool operator < (const Element &lh, const Element &rh) {
  return lh.row < rh.row || (lh.row == rh.row && lh.col < rh.col);
}

bool operator ==(const Element &lh, const Element &rh) {
  return lh.row == rh.row && lh.col == rh.col && lh.val == rh.val;
}






// given to vector, return their inner product
int inner_product(const CMatrix_CSC& a, const CMatrix_CSC& b, double speedup_thresh) {
  int pa = 0;
  int pb = 0;
  int res = 0;
  while (pa < a.col_ptr[1] && pb < b.col_ptr[1]) {
    if (a.row[pa] == b.row[pb]) {
      res += a.val[pa++] * b.val[pb++];
      if (res >= speedup_thresh)
	return res;
    }
    else if (a.row[pa] < b.row[pb]) {
      ++pa;
    }
    else
      ++pb;
  }
  return res;
}








