#include "Algebra.h"
#include "Sketch.h"
#include <numeric>
//#include <iomanip>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <map>
#include "utils.h"


const int _INFINITY = 100000000;
const int MAX_LOGN = 25;
const int mu = 2;

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



CVector calcL1Norm(const CMatrix_CSC& P, const CMatrix_CSC& Q) {
  int val[P.m];
  std::fill(val, val + P.m, 1);
  int row[P.m];
  std::fill(row, row + P.m, 0);
  int col[P.m];
  std::iota(col, col + P.m, 0);
  CMatrix_CSC S(val, row, col, P.m, 1, P.m);
  CMatrix_COO&& res = toCoo((S * P) * Q);

  
  int l1[Q.n];
  std::fill(l1, l1 + Q.n, 0);

  for (int i = 0; i < res.size(); ++i)
    l1[res[i].col] += res[i].val;

  return CVector(l1, l1 + Q.n);
}

// create count min with dim w * mu * n
CMatrix_CSC createCountMin(int w, int mu, int n) {
  std::srand(time(NULL));
  int *val = new int[mu * n];
  std::fill(val, val + mu * n, 1);
  int *row = new int[mu * n];
  int *col = new int[mu * n];
  int t = 0;
  for (int i = 0; i < n; ++i) 
    for (int j = 0; j < mu; ++j) {
      row[t] = std::rand() % w + j * w;
      col[t++] = i;
    }
      
  CMatrix_CSC res(val, row, col, mu * n, w * mu, n, true);
  delete[] val;
  delete[] row;
  delete[] col;
  return res;
}


// recover entries > theta in P * Q * W
CMatrix_CSC FastThreshMult(const CMatrix_CSC &P, const CMatrix_CSC &Q, 
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
  //  std::cerr << "w = " << w << " mu = " << mu << std::endl;

  // create sketches
  CMatrix_CSC* CMs[MAX_LOGN];
  CMatrix_CSC sk[MAX_LOGN];


  timer.start();
  for (int i = 0; i < t; ++i) {
    // create count min sketch
    CMs[i] = new CMatrix_CSC(createCountMin(w, mu, Ps[i]->m));
    sk[i] = ((*CMs[i] * *Ps[i]) * Q);// * W;
  }
  timer.stop("Create sketches  ");


  
  timer.start();
  CVector&& R1 = calcL1Norm(P, Q);
  timer.stop("calcL1Norm(P, Q) ");


  int debug_ct_naive = 0;
  int debug_ct_algor = 0;

  
  timer.start();
  ///////////////////////////////////////////
  // recover heavy entries
  ///////////////////////////////////////////
  CVector val, row, col;
  int skCol[w * mu]; // keep extracted column vector

  for (int i = 0; i < Q.n; ++i) { // for column

    if (R1[i] > (0.6 + rho) * w * theta) { // use exact algorithm
      debug_ct_naive++;

      std::map<int, int> res;
      
      for (int j = Q.col_ptr[i]; j < Q.col_ptr[i + 1]; ++j) 
      	for (int k = P.col_ptr[Q.row[j]]; k < P.col_ptr[Q.row[j] + 1]; ++k)
      	  res[P.row[k]] += P.val[k] * Q.val[j];

      for (auto it = res.begin(); it != res.end(); ++it) 
      	if (it->second > theta) {
      	val.push_back(it->second);
      	row.push_back(it->first);
      	col.push_back(i);
      }

      // CMatrix_CSC&& res = P * Q[i];
      // for (int j = 0; j < res.col_ptr[1]; ++j) 
      // 	if (res.val[j] > theta) {
      // 	  val.push_back(res.val[j]);
      // 	  row.push_back(res.row[j]);
      // 	  col.push_back(i);
      // 	}
      continue;
    }

    debug_ct_algor++;

    CVector ind[MAX_LOGN]; // dyadic structure
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
  
  std::cerr << debug_ct_naive << " columns use Naive, " << debug_ct_algor 
	    << " columns use Sketch." << std::endl;

  return CMatrix_CSC(val.begin(), row.begin(), col.begin(), val.size(), P.m, Q.n);
}


