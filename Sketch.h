#ifndef __SKETCH_H__
#define  __SKETCH_H__
#include "Algebra.h"
#include <vector>
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
			       double theta, double rho, int w);


// recover the entry rw
VAL_TYPE recover(VAL_TYPE sk[], const CMatrix_CSC &cm, int rw);

// convert csc[i] to int[], keep in arr[]
void slicing(VAL_TYPE arr[], const CMatrix_CSC& csc, int i);


#endif


















