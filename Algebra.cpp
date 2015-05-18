#include "Algebra.h"
#include <iomanip>
#include <algorithm>
#include <ctime>




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




std::ostream& operator << (std::ostream & os, CMatrix_COO &coo) {
  int A[coo.get_m()][coo.get_n()];
  std::fill(A[0], A[0] + coo.get_m() * coo.get_n(), 0);
  for (auto it = coo.data.begin(); it != coo.data.end(); ++it)
    A[it->row][it->col] = it->val;
  for (int i = 0; i < coo.get_m(); ++i) {
    for (int j = 0; j < coo.get_n(); ++j)
      os << std::setw(6) << A[i][j];
    os << std::endl;
  }
  //  os << it->row << ", " << it->col << ", " << it->val << std::endl;
  return os;
}




// product of two coo matrix, return a int*
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



// r = P * vec
// result = r[r > thresh]
// note:  P has to be sorted by RowColumn
void thresholdMult(SparseVector &result, CMatrix_COO &P, SparseVector &vec, double thresh) {
  int v[vec.size()];
  std::fill(v, v + vec.size(), 0);
  for (const VectorElement &it : vec)
    v[it.ind] = it.val;
  int t = 0; 
  while (t < P.size()) {
    int current_row = P[t].row;
    int prod = 0;
    prod += P[t].val * v[P[t].col];
    ++t;

    while (t < P.size() && P[t].row == P[t-1].row) {
      prod += P[t].val * v[P[t].col];
      ++t;
    }

    if (prod > thresh) 
      result.push_back(VectorElement(current_row, prod));    
  }
}





std::ostream& operator << (std::ostream &os,  const SparseVector &vec) {
  for (const VectorElement &elem : vec)
    os << elem.ind << ", " << elem.val << std::endl;
  return os;
}

















