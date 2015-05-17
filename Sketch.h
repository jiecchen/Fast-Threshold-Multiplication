#ifndef __SKETCH_H__
#define  __SKETCH_H__
#include "Algebra.h"
#include <vector>
const int MAX_LOGN = 20;

/////////////////////////////////////////////////////////////////////////
/////////////////////// Sketch Related //////////////////////////////////
/////////////////////////////////////////////////////////////////////////
class MatrixSketch {
public:
 /* MatrixSketch(int w, int u, int n, CMatrix_COO* hash, int* data): */
 /*  w(w), u(u), n(n), hash(hash), data(data) {} */
  ~MatrixSketch() {
    delete[] data;
    delete hash;
  };
  int w, u; // h = w * u
  int n;
  CMatrix_COO *hash;
  int *data; // will be h * n
};



// sketch to keep the dyadic structure
typedef std::vector<MatrixSketch> DyadicSketch;
typedef std::vector<MatrixSketch>::iterator DyadicSketch_iter;



// MatrixSketch * Vector -> CountMinSketch
// sketching a coo matrix
// keep result in MatrixSketch
void sketchMatrixCOO(MatrixSketch &sk, double eps, int u, CMatrix_COO &B);
// return int* as count min sketch of the vector specfied by it_s -> it_e
// remember to release the memory
int* sketchVector(MatrixSketch &sk, SparseVector_iter it_s, SparseVector_iter it_e);
// recover the coordinate given sketch and hash
int recover(int *cm, const MatrixSketch &sk, int coor);


// let oldCoo = [R_0; R_1; R_2; R_3; ...]
// newCoo = [R_0 + R_1; R_2 + R_3; ...]
// oldCoo has to be sorted by ColumnRow
// newCoo is empty
void mergeNeighbor(CMatrix_COO &newCoo, CMatrix_COO &oldCoo);




#endif


















