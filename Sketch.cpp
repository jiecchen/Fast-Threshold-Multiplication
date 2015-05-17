#include "Algebra.h"
#include "Sketch.h"
#include <cmath>

//Principle:
//   + each function only do one thing
//   + keep simple, keep stupid



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


// let oldCoo = [R_0; R_1; R_2; R_3; ...]
// newCoo = [R_0 + R_1; R_2 + R_3; ...]
// oldCoo has to be sorted by ColumnRow
// newCoo is empty
void mergeNeighbor(CMatrix_COO &newCoo, CMatrix_COO &oldCoo) {
  newCoo.set_mn((oldCoo.get_m() + 1) / 2, oldCoo.get_n());
  if (oldCoo.size() <= 0)
    return;
  newCoo.push_back(Element(oldCoo[0].row / 2, oldCoo[0].col, oldCoo[0].val));
  for (int i = 1; i < oldCoo.size(); ++i)
    if (oldCoo[i].col == newCoo.back().col &&
	oldCoo[i].row / 2 == newCoo.back().row) {
      newCoo.back().val += oldCoo[i].val;
    }
    else
      newCoo.push_back(Element(oldCoo[i].row / 2, oldCoo[i].col, oldCoo[i].val));
}

// create O(log m) matrix sketches, form as
// dyadic sketch
void dyadicSketch(CDyadicSketch &dyadic,  double eps, int u, CMatrix_COO &B) {
  dyadic.clear();
  CMatrix_COO coo[MAX_LOGN];
  coo[0] = B;
  coo[0].sortByColumnRow();
  int t = 0;
  //TODO:
  //  + bugs, if coo[t].size() == 0
  std::cout << "coo[" << t << "]: 0 \n" << coo[0] << std::endl;
  while (coo[t].get_m() > 1) {
    mergeNeighbor(coo[t + 1], coo[t]);
    t++;
    std::cout << "coo[" << t << "]: " << coo[t].get_m() << "\n" << coo[t] << std::endl;
  }
  while (t >= 0) {
    MatrixSketch *sk = new MatrixSketch();
    sketchMatrixCOO(*sk, eps, u, coo[t--]);
    dyadic.push_back(sk);
  }
}
















