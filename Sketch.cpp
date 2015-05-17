#include "Algebra.h"
#include "Sketch.h"
#include <cmath>

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

