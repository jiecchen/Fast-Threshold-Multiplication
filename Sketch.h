#ifndef __SKETCH_H__
#define  __SKETCH_H__
#include "Algebra.h"
#include <vector>

const int STEP_SIZE = 2; // how many neighbors to be merged
const VAL_TYPE _INFINITY = 1 << 30;
const int MAX_LOGN = 25;
const int _MU = 3;


/////////////////////////////////////////////////////////////////////////
/////////////////////// Sketch Related //////////////////////////////////
/////////////////////////////////////////////////////////////////////////




CMatrix_CSC FastThreshMult_filter(const CMatrix_CSC &P, const CMatrix_CSC &Q, 
				  double theta, double rho, int w=20);


CMatrix_CSC FastThreshMult(const CMatrix_CSC &P, const CMatrix_CSC &Q,
			   double theta, double rho, int w=0);

CMatrix_CSC FastThreshMult_Simple(const CMatrix_CSC &P, const CMatrix_CSC &Q, 
				  double theta, double rho, int w=0);

// new version of FastThreshMult
CMatrix_CSC FastThreshMult_new(const CMatrix_CSC &P, const CMatrix_CSC &Q, 
			       VAL_TYPE theta, int w);

//merge rows of matrix
CMatrix_CSC mergeNeighbor(const CMatrix_CSC &oldCsc, int step_size=STEP_SIZE);

// create count min with dim w * mu * n
CMatrix_CSC createCountMin(int w, int mu, int n);

// given count-min and sketches, recover large (> theta) entries for i_th column 
// only return k values if have
CMatrix_COO recover_column(const std::vector<CMatrix_CSC*>& CMs, 
			   const std::vector<CMatrix_CSC>& sk, 
			   int i, VAL_TYPE theta, int k = 1 << 30);

// recover the entry rw
VAL_TYPE recover(VAL_TYPE sk[], const CMatrix_CSC &cm, int rw);

// convert csc[i] to int[], keep in arr[]
void slicing(VAL_TYPE arr[], const CMatrix_CSC& csc, int i);


#endif


















