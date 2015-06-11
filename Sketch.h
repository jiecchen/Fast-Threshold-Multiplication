#ifndef __SKETCH_H__
#define  __SKETCH_H__
#include "Algebra.h"
#include <vector>
/////////////////////////////////////////////////////////////////////////
/////////////////////// Sketch Related //////////////////////////////////
/////////////////////////////////////////////////////////////////////////


// let oldCoo = [R_0; R_1; R_2; R_3; ...]
// newCoo = [R_0 + R_1; R_2 + R_3; ...]
// oldCoo has to be sorted by ColumnRow
// newCoo is empty
//void mergeNeighbor(CMatrix_COO &newCoo, CMatrix_COO &oldCoo);
//CMatrix_CSC mergeNeighbor(const CMatrix_CSC &oldCsc);


CMatrix_CSC FastThreshMult(const CMatrix_CSC &P, const CMatrix_CSC &Q,
			   int w, double theta, double rho);

CMatrix_CSC FastThreshMult_Simple(const CMatrix_CSC &P, const CMatrix_CSC &Q, 
				  int w, double theta, double rho);


/* // P, Q are boolean matrices */
/* // theta is a threshold > 0 */
/* // rho \in (0, 1) to control the accuracy */
/* CMatrix_COO atLeastMult(CMatrix_COO &P, CMatrix_COO &Q, double theta, double rho); */





/* // double calcL1Norm(const CMatrix_CSC& P, const CMatrix_CSC& Q, const CMatrix_CSC& W); */
/* CMatrix_CSC createCountMin(int w, int mu, int n); */




#endif


















