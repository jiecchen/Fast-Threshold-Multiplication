#include "Algebra.h"
#include "Sketch.h"
#include <numeric>
#include <limits>
//#include <iomanip>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <map>
#include "utils.h"


const int _INFINITY = std::numeric_limits<int>::infinity();
const int MAX_LOGN = 25;
const int mu = 6;
const double alpha = 1.;
const int STEP_SIZE = 500; // how many neighbors to be merged


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
	row.push_back(oldCsc.row[j] / STEP_SIZE);
	col.push_back(i);
	continue;
      }
      
      if ((oldCsc.row[j] / STEP_SIZE) == row.back()) 
	val.back() += oldCsc.val[j];
      else {
	val.push_back(oldCsc.val[j]);
	row.push_back(oldCsc.row[j] / STEP_SIZE);
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



// return w as # of buckets to optimize the running time
int optimalBucketsSize(const CVector& R1, const CMatrix_CSC& Q, double rho, double theta) {
  SparseVector R_1; // to keep R1
  SparseVector Q_0; // to keep ||Q||_1
  for (unsigned int i = 0; i < R1.size(); ++i)
    R_1.push_back(VectorElement(i, R1[i]));
  std::sort(
	    R_1.begin(), R_1.end(), 
	    [](const VectorElement& a, const VectorElement& b) {
	      return a.val > b.val;
	    }
	    );


  int tot = 0;
  for (auto it = R_1.begin(); it != R_1.end(); ++it) {
    tot += it->val;
    it->val = tot;
    Q_0.push_back(VectorElement(it->ind, Q.get_nnz(it->ind)));
  }
  
  tot = 0;
  for (auto it = Q_0.end() - 1; it + 1 != Q_0.begin(); --it) {
    tot += it->val;
    it->val = tot;
  }
  Q_0.push_back(VectorElement(-1, 0));



  double min_cost = std::numeric_limits<double>::infinity();
  double w = 0;

  for (unsigned int t = 0; t < R_1.size(); ++t) {
    double _w = 4. * R1[R_1[t + 1].ind] / ((rho + alpha) * theta);
    double new_cost = 25 * R_1[t].val  +  mu * _w * Q_0[t + 1].val;
    // std::cerr << R_1[t].val << "  " << R1[R_1[t + 1].ind] << "  " << _w << "  " <<  Q_0[t + 1].val << "  new_cost = " << new_cost << std::endl;
    if (_w < 5)
      _w = 5;
    if (new_cost < min_cost) {
      min_cost = new_cost;
      w = _w;
      //      std::cerr << "w <- " << w << " cost <- " << min_cost  << " t = " << t << std::endl;
    }
  }
  
  std::cerr << "Optimal w = " << (int) w << std::endl;
  return (int) w;
}


// ugly, put all thing together to get speedup
// recover entries > theta in P * Q * W
CMatrix_CSC FastThreshMult(const CMatrix_CSC &P, const CMatrix_CSC &Q, 
			   int w, double theta, double rho) {


  CMatrix_CSC* Ps[MAX_LOGN]; // keep merged matrices
  CTimer timer;
  int t = 0; // # of levels in dyadic structure
  Ps[0] = new CMatrix_CSC(P);




  while (Ps[t]->m > 1) {
    Ps[t + 1] = new CMatrix_CSC(mergeNeighbor(*Ps[t]));
    ++t;
  }
  ++t;



  // calc L1 norm of P * Q
  CVector&& R1 = calcL1Norm(P, Q);

  int _w = optimalBucketsSize(R1, Q, rho, theta);
  w = _w;

  // for debug purpose
  int debug_ct_naive = 0;
  int debug_ct_algor = 0;

  
  ///////////////////////////////////////////
  // recover heavy entries
  ///////////////////////////////////////////
  CVector val, row, col;
  int skCol[w * mu]; // keep extracted column vector

  std::vector<int> use_exact;
  std::vector<int> use_sketch;
  for (int i = 0; i < Q.n; ++i) // for column
    if (R1[i] > (alpha + rho) / 4. * w * theta) { // use exact algorithm
      debug_ct_naive++;
      use_exact.push_back(i);
    }
    else { // use sketch
      debug_ct_algor++;
      use_sketch.push_back(i);
    }

  std::cerr << debug_ct_naive << " columns use Naive, " << debug_ct_algor 
	    << " columns use Sketch." << std::endl;
  
 
  //////////////////////////////////////////////////////
  ///////////////// Using exact algo ///////////////////
  //////////////////////////////////////////////////////
  timer.start("Use exact algo");
  // slicing
  timer.start("exact matrix multilication");
  CMatrix_CSC&& exact_part = P * (Q[use_exact]);
  timer.stop();
  // for the exact calculation part, append entries > theta to result vector
  timer.start("pick up heavy hitters");
  int ptr = 0;
  for (auto it = use_exact.begin(); it != use_exact.end(); ++it) {
    for (int i = exact_part.col_ptr[ptr]; i < exact_part.col_ptr[ptr + 1]; ++i) 
      if (exact_part.val[i] > theta) {
	val.push_back(exact_part.val[i]);
	row.push_back(exact_part.row[i]);
	col.push_back(*it);
      }
    ++ptr;
  }
  timer.stop();
  timer.stop();



  /////////////////////////////////////////////////////
  /////////////  Using sketch /////////////////////////
  /////////////////////////////////////////////////////
  

  timer.start("Use sketch");
  // create sketch 
  CMatrix_CSC&& slice = Q[use_sketch];
  CMatrix_CSC sk[MAX_LOGN];
  CMatrix_CSC* CMs[MAX_LOGN];


  timer.start("sketch matrix");
  timer.start("Create count-min & sketch P");
  for (int i = 0; i < t; ++i) {
    // create count min sketch
    CMs[i] = new CMatrix_CSC(createCountMin(w, mu, Ps[i]->m));
    sk[i] = (*CMs[i] * *Ps[i]); // * W;
  }
  timer.stop();
  timer.start("Sketch Q");
  for (int i = 0; i < t; ++i)
    sk[i] = sk[i] * slice;
  timer.stop();
  timer.stop();


  timer.start("Do recovering");
  // dealing with use_sketch part
  ptr = 0;
  for (auto i = use_sketch.begin(); i != use_sketch.end(); ++i) {
    CVector ind[MAX_LOGN]; // dyadic structure
    ind[t - 1].push_back(0);
    for (int j = t - 1; j >= 0; --j) {

      // convert csc col to int[]
      std::fill(skCol, skCol + w * mu, 0);
      for (int k = sk[j].col_ptr[ptr]; k < sk[j].col_ptr[ptr + 1]; ++k)
	skCol[sk[j].row[k]] += sk[j].val[k];
    

      // now do recover
      for (auto it = ind[j].begin(); it != ind[j].end(); ++it) {
	if (*it >= Ps[j]->m)
	  break;
	int v = recover(skCol, *CMs[j], *it);
	if (v > theta) {
	  if (j > 0) { // has not reached level 0  
	    for (int j0 = 0; j0 < STEP_SIZE; ++j0)
	      ind[j - 1].push_back((*it) * STEP_SIZE + j0);
	    //ind[j - 1].push_back((*it) * 2);
	    //ind[j - 1].push_back((*it) * 2 + 1);
	  }
	  else { // reach level 0, keep the result
	    val.push_back(v);
	    row.push_back(*it);
	    col.push_back(*i);
	  } 
	} // if
      } // for (auto
    } // for (int j ..      
    ptr++;
  }
  timer.stop();


  // release memory
  while (--t >= 0) {
    delete Ps[t];
  }

  timer.stop();
  
  return CMatrix_CSC(val.begin(), row.begin(), col.begin(), val.size(), P.m, Q.n);
}




CMatrix_CSC FastThreshMult_Simple(const CMatrix_CSC &P, const CMatrix_CSC &Q, 
			   int w, double theta, double rho) {

  // calc L1 norm of P * Q
  CVector&& R1 = calcL1Norm(P, Q);

  CTimer timer;
  // for debug purpose
  int debug_ct_naive = 0;
  int debug_ct_algor = 0;

  
  ///////////////////////////////////////////
  // recover heavy entries
  ///////////////////////////////////////////
  CVector val, row, col;
  int skCol[w * mu]; // keep extracted column vector

  std::vector<int> use_exact;
  std::vector<int> use_sketch;
  for (int i = 0; i < Q.n; ++i) // for column
    if (R1[i] > (alpha + rho) / 4. * w * theta) { // use exact algorithm
      debug_ct_naive++;
      use_exact.push_back(i);
    }
    else { // use sketch
      debug_ct_algor++;
      use_sketch.push_back(i);
    }

  std::cerr << debug_ct_naive << " columns use Naive, " << debug_ct_algor 
	    << " columns use Sketch." << std::endl;

  
 
  //////////////////////////////////////////////////////
  ///////////////// Using exact algo ///////////////////
  //////////////////////////////////////////////////////
  timer.start("Use exact algo");
  // slicing
  CMatrix_CSC&& exact_part = P * (Q[use_exact]);
  // for the exact calculation part, append entries > theta to result vector
  int ptr = 0;
  for (auto it = use_exact.begin(); it != use_exact.end(); ++it) {
    for (int i = exact_part.col_ptr[ptr]; i < exact_part.col_ptr[ptr + 1]; ++i) 
      if (exact_part.val[i] > theta) {
	val.push_back(exact_part.val[i]);
	row.push_back(exact_part.row[i]);
	col.push_back(*it);
      }
    ++ptr;
  }
  timer.stop();



  /////////////////////////////////////////////////////
  /////////////  Using sketch /////////////////////////
  /////////////////////////////////////////////////////
  

  timer.start("Create sketches");
  // create sketch 
  CMatrix_CSC&& slice = Q[use_sketch];
  timer.stop();

  timer.start("Skech the matrix");
  CMatrix_CSC CM(createCountMin(w, mu, P.m)); 
  CMatrix_CSC sk((CM * P) * slice);
  timer.stop();


  timer.start("Use sketch to recover");
  // dealing with use_sketch part
  ptr = 0;
  for (auto i = use_sketch.begin(); i != use_sketch.end(); ++i) {
    // convert csc col to int[]
    std::fill(skCol, skCol + w * mu, 0);
    for (int k = sk.col_ptr[ptr]; k < sk.col_ptr[ptr + 1]; ++k)
      skCol[sk.row[k]] += sk.val[k];
    

    // now do recover
    for (auto it = 0; it < P.m; ++it) {
      int v = recover(skCol, CM, it);
	if (v > theta) {
	  val.push_back(v);
	  row.push_back(it);
	  col.push_back(*i);
	} 
    } // for (auto
    ptr++;
  } // for (int i
    
  timer.stop();
  
  return CMatrix_CSC(val.begin(), row.begin(), col.begin(), val.size(), P.m, Q.n);
}
