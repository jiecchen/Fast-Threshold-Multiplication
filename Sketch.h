#ifndef __SKETCH_H__
#define  __SKETCH_H__
#include "Algebra.h"
#include <vector>
const int MAX_LOGN = 15;

/////////////////////////////////////////////////////////////////////////
/////////////////////// Sketch Related //////////////////////////////////
/////////////////////////////////////////////////////////////////////////
class MatrixSketch {
public:
 /* MatrixSketch(int w, int u, int n, CMatrix_COO* hash, int* data): */
 /*  w(w), u(u), n(n), hash(hash), data(data) {} */

  //TODO:
  //  + resolve the memory leak issue
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
typedef std::vector<MatrixSketch *> DyadicSketch;
typedef std::vector<MatrixSketch *>::iterator DyadicSketch_iter;

class CDyadicSketch {
public:
  ~CDyadicSketch() {
    for (auto it = sk.begin(); it != sk.end(); ++it)
      delete *it;
  }
  void push_back(MatrixSketch *p) {
    sk.push_back(p);
  }
  int size() { return sk.size(); }
  void clear() { sk.clear(); }
  MatrixSketch* operator[] (int i) { return sk[i]; }

private:
  DyadicSketch sk;
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


// let oldCoo = [R_0; R_1; R_2; R_3; ...]
// newCoo = [R_0 + R_1; R_2 + R_3; ...]
// oldCoo has to be sorted by ColumnRow
// newCoo is empty
void mergeNeighbor(CMatrix_COO &newCoo, CMatrix_COO &oldCoo);
CMatrix_CSC mergeNeighbor(const CMatrix_CSC &oldCsc);


CMatrix_CSC FastThreshMult(const CMatrix_CSC &P, const CMatrix_CSC &Q, const CMatrix_CSC& W,
			    double theta, double rho);
// create a dyadic structure for B
// keep in dyadic
void dyadicSketch(CDyadicSketch &dyadic,  double eps, int u, CMatrix_COO &B);

// given threshhold, recover all
// entries in dyadic(vec_s, vec_e)
// approximately recover entries in PQ that >= thresh
// keep in result
void thresholdRecover(SparseVector &result, CDyadicSketch &dyadic, 
		      SparseVector_iter vec_s, SparseVector_iter vec_e, double thresh);


// P, Q are boolean matrices
// theta is a threshold > 0
// rho \in (0, 1) to control the accuracy
CMatrix_COO atLeastMult(CMatrix_COO &P, CMatrix_COO &Q, double theta, double rho);





// double calcL1Norm(const CMatrix_CSC& P, const CMatrix_CSC& Q, const CMatrix_CSC& W);
CMatrix_CSC createCountMin(int w, int mu, int n);




#endif


















