#include "Algebra.h"
#include "Sketch.h"
#include <numeric>
//#include <iomanip>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include "utils.h"

const int _INFINITY = 100000000;

//Principle:
//   + each function only do one thing
//   + keep simple, keep stupid



/////////////////////////////////////////////////////////////////////////
/////////////////////// Sketch Related //////////////////////////////////
/////////////////////////////////////////////////////////////////////////
// recover the entry rw
int recover(int sk[], const CMatrix_CSC &cm, int rw) {
  int m = _INFINITY;
  for (int i = 0; i < mu; ++i) {
    int r = cm.row[rw * mu + i];
    m = std::min(sk[r], m);
  }
  //  std::cerr << "m = " << m << std::endl;
  return m;
}


// let oldCsc = [R_0; R_1; R_2; R_3; ...]
// newCsc = [R_0 + R_1; R_2 + R_3; ...]
CMatrix_CSC mergeNeighbor(const CMatrix_CSC &oldCsc) {
  if (oldCsc.nnz <= 0)
    return CMatrix_CSC();
  CVector val, row, col;
  for (int i = 0; i < oldCsc.n; ++i) { // col number

    for (int j = oldCsc.col_ptr[i]; j < oldCsc.col_ptr[i + 1]; ++j) {
      if (j == oldCsc.col_ptr[i]) {
	val.push_back(oldCsc.val[j]);
	row.push_back(oldCsc.row[j] >> 1);
	col.push_back(i);
	continue;
      }
      
      if ((oldCsc.row[j] >> 1) == row.back()) 
	val.back() += oldCsc.val[j];
      else {
	val.push_back(oldCsc.val[j]);
	row.push_back(oldCsc.row[j] >> 1);
	col.push_back(i);
      }
    }// for j    
  }// for i
  
  return CMatrix_CSC(val.begin(), row.begin(), col.begin(), val.size(), 
		     (oldCsc.m + 1) >> 1, oldCsc.n); 
}



CVector calcL1Norm(const CMatrix_CSC& P, const CMatrix_CSC& Q, const CMatrix_CSC& W) {
  int val[P.m];
  std::fill(val, val + P.m, 1);
  int row[P.m];
  std::fill(row, row + P.m, 0);
  int col[P.m];
  std::iota(col, col + P.m, 0);
  CMatrix_CSC S(val, row, col, P.m, 1, P.m);
  CMatrix_COO&& res = toCoo(((S * P) * Q) * W);

  
  int l1[W.n];
  std::fill(l1, l1 + W.n, 0);

  for (int i = 0; i < res.size(); ++i)
    l1[res[i].col] += res[i].val;

  return CVector(l1, l1 + W.n);
}

// create count min with dim w * mu * n
CMatrix_CSC createCountMin(int w, int mu, int n) {
  std::srand(time(NULL));
  int val[mu * n];
  std::fill(val, val + mu * n, 1);
  int row[mu * n];
  int col[mu * n];
  int t = 0;
  for (int i = 0; i < n; ++i) 
    for (int j = 0; j < mu; ++j) {
      row[t] = std::rand() % w + j * w;
      col[t++] = i;
    }
      
  return CMatrix_CSC(val, row, col, mu * n, w * mu, n);
}


// recover entries > theta in P * Q * W
CMatrix_CSC FastThreshMult(const CMatrix_CSC &P, const CMatrix_CSC &Q, const CMatrix_CSC &W, 
			   int w, double theta, double rho) {

  CMatrix_CSC* Ps[MAX_LOGN];
  CTimer timer;
  int t = 0;
  Ps[0] = new CMatrix_CSC(P);
  //  std::cerr << "Ps[0] = \n" << toCoo(*Ps[0]) << std::endl;
  while (Ps[t]->m > 1) {
    Ps[t + 1] = new CMatrix_CSC(mergeNeighbor(*Ps[t]));
    ++t;
  }
  ++t;

  // calc ||P * Q  * W||_1
  //  double R1 = calcL1Norm(P, Q, W);
  // int w = std::ceil(0.1 * R1 / theta / rho);
  std::cerr << "w = " << w << " mu = " << mu << std::endl;

  // create sketches
  CMatrix_CSC* CMs[MAX_LOGN];
  CMatrix_CSC sk[MAX_LOGN];

  

  timer.start();
  for (int i = 0; i < t; ++i) {
    // create count min sketch
    CMs[i] = new CMatrix_CSC(createCountMin(w, mu, Ps[i]->m));
    sk[i] = ((*CMs[i] * *Ps[i]) * Q) * W;
  }
  timer.stop("Create sketches  ");

  
  timer.start();
  CVector&& R1 = calcL1Norm(P, Q, W);
  timer.stop("calcL1Norm(P, Q, W) ");
  
  timer.start();
  // recover heavy entries
  CVector val, row, col;
  int skCol[w * mu]; // keep extracted column vector
  for (int i = 0; i < W.n; ++i) { // for column

    if (2 * R1[i] > 5 * (2 + rho * (w - 2)) * theta) { // use exact algorithm
      CMatrix_CSC&& res = P * (Q * W[i]);
      for (int j = 0; j < res.nnz; ++j)
	if (res.val[j] > theta) {
	  val.push_back(res.val[j]);
	  row.push_back(res.row[j]);
	  col.push_back(i);
	}
      //      std::cerr << "Column " << i << " uses exact algorithm" << std::endl;
      continue;
    }
    
    //    std::cerr << "Column " << i << " uses count min" << std::endl;
    // use count min
    CVector ind[MAX_LOGN];
    ind[t - 1].push_back(0);
    for (int j = t - 1; j >= 0; --j) {
      // convert csc col to int[]
      std::fill(skCol, skCol + w * mu, 0);
      for (int k = sk[j].col_ptr[i]; k < sk[j].col_ptr[i + 1]; ++k)
	skCol[sk[j].row[k]] += sk[j].val[k];
      // now do recover
      for (auto it = ind[j].begin(); it != ind[j].end(); ++it) {
	if (*it >= Ps[j]->m)
	  break;
	int v = recover(skCol, *CMs[j], *it);
	if (v > theta) {
	  if (j > 0) { // has not reached level 0  
	    ind[j - 1].push_back((*it) * 2);
	    ind[j - 1].push_back((*it) * 2 + 1);
	  }
	  else { // reach level 0, keep the result
	    val.push_back(v);
	    row.push_back(*it);
	    col.push_back(i);

	  } 
	} // if
      } // for (auto
    } // for (int j ..      
  }// for (int i
  timer.stop("Recover heavy coordinates ");


  // release memory
  while (--t >= 0) {
    delete Ps[t];
  }

  return CMatrix_CSC(val.begin(), row.begin(), col.begin(), val.size(), P.m, W.n);
}


