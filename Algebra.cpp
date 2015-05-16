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
  for (auto it = coo.data.begin(); it != coo.data.end(); ++it)
    os << it->row << ", " << it->col << ", " << it->val << std::endl;
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


/////////////////////////////////////////////////////////////////////////
/////////////////////// Sketch Related //////////////////////////////////
/////////////////////////////////////////////////////////////////////////



// sketching a coo matrix
// keep result in MatrixSketch
void sketchMatrixCOO(MatrixSketch &sk, double eps, int u, CMatrix_COO &B) {
  std::srand(time(NULL));
  sk.w = std::ceil(1. / eps);
  sk.u = u;
  int h = sk.w * u;
  sk.n = B.get_n();
  // create count min sketch
  CMatrix_COO* cm = new CMatrix_COO(h, B.get_m());
  for (int i = 0; i < B.get_m(); ++i) 
    for (int j = 0; j < u; ++j) 
      cm->push_back(Element(std::rand() % sk.w + j * sk.w, i, 1));

  sk.hash = cm;
  sk.data = coo_mult(*cm, B);
  sk.hash->sortByColumnRow();
}


int* sketchVector(MatrixSketch &sk, SparseVector_iter it_s, SparseVector_iter it_e) {
  int h = sk.w * sk.u;
  // init as 0
  int *cm = new int[h]();
  for (int i = 0; i < h; i++)
    for (auto it = it_s; it != it_e; ++it)
      cm[i] += sk.data[i * sk.n + it->ind] * it->val;
  return cm;
}


int recover(int *cm, const MatrixSketch &sk, int coor) {
  int m = INFINITY;
  CMatrix_COO &hash = *(sk.hash);
  for (int i = 0; i < sk.u; ++i) {
    int r = hash[coor * sk.u + i].row;
    m = std::min(cm[r], m);
  }
  return m;
}
















