#ifndef __SKETCH_H__
#define  __SKETCH_H__
#include "Algebra.h"


/////////////////////////////////////////////////////////////////////////
/////////////////////// Sketch Related //////////////////////////////////
/////////////////////////////////////////////////////////////////////////
class MatrixSketch {
public:
  ~MatrixSketch() {
    delete[] data;
    delete hash;
  };
  int w, u; // h = w * u
  int n;
  CMatrix_COO *hash;
  int *data; // will be h * n
};

// MatrixSketch * Vector -> CountMinSketch
// sketching a coo matrix
// keep result in MatrixSketch
void sketchMatrixCOO(MatrixSketch &sk, double eps, int u, CMatrix_COO &B);
// return int* as count min sketch of the vector specfied by it_s -> it_e
// remember to release the memory
int* sketchVector(MatrixSketch &sk, SparseVector_iter it_s, SparseVector_iter it_e);
// recover the coordinate given sketch and hash
int recover(int *cm, const MatrixSketch &sk, int coor);


void dyadicSketch(double eps, int u, CMatrix_COO &B);




#endif

















