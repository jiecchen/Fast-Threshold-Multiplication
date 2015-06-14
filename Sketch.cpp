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


const int _INFINITY = 1 << 30;
const int MAX_LOGN = 25;
const int mu = 3;
const double alpha = 0.5;
int STEP_SIZE = 50; // how many neighbors to be merged


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



// convert csc[i] to int[], keep in arr[]
void slicing(int arr[], const CMatrix_CSC& csc, int i) {
  // convert csc col to int[]
  std::fill(arr, arr + csc.m, 0);
  for (int k = csc.col_ptr[i]; k < csc.col_ptr[i + 1]; ++k)
    arr[csc.row[k]] = csc.val[k];
}







// let oldCsc = [R_0; R_1; R_2; R_3; ...]
// newCsc = [R_0 + R_1; R_2 + R_3; ...]
CMatrix_CSC mergeNeighbor(const CMatrix_CSC &oldCsc, int step_size=STEP_SIZE) {
  if (oldCsc.nnz <= 0)
    return CMatrix_CSC();
  CVector val, row, col;
  for (int i = 0; i < oldCsc.n; ++i) { // col number

    for (int j = oldCsc.col_ptr[i]; j < oldCsc.col_ptr[i + 1]; ++j) {
      if (j == oldCsc.col_ptr[i]) {
	val.push_back(oldCsc.val[j]);
	row.push_back(oldCsc.row[j] / step_size);
	col.push_back(i);
	continue;
      }
      
      if ((oldCsc.row[j] / step_size) == row.back()) 
	val.back() += oldCsc.val[j];
      else {
	val.push_back(oldCsc.val[j]);
	row.push_back(oldCsc.row[j] / step_size);
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



  double min_cost = 1e+15;
  double w = 0;

  for (unsigned int t = 0; t < R_1.size(); ++t) {
    double _w = 4. * R1[R_1[t + 1].ind] / ((rho + alpha) * theta);
    double new_cost =  2.9 * R_1[t].val  +  mu * _w * Q_0[t + 1].val;
    // std::cerr << R_1[t].val << "  " << R1[R_1[t + 1].ind] << "  " << _w << "  " <<  Q_0[t + 1].val << "  new_cost = " << new_cost << std::endl;
    if (_w < 1)
      _w = 1;
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
			   double theta, double rho, int w) {

  // set up the step_size to merge neighbor
  int step_size = (int) std::min(std::sqrt(P.n) * 20 + 1, P.n / 2.);

  CMatrix_CSC* Ps[MAX_LOGN]; // keep merged matrices
  CTimer timer;
  int t = 0; // # of levels in dyadic structure
  Ps[0] = new CMatrix_CSC(P);




  while (Ps[t]->m > 1) {
    Ps[t + 1] = new CMatrix_CSC(mergeNeighbor(*Ps[t], step_size));
    ++t;
  }
  ++t;



  // calc L1 norm of P * Q
  CVector&& R1 = calcL1Norm(P, Q);
  if (w < 1)
    w = optimalBucketsSize(R1, Q, rho, theta);
  
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
  for (int i = 0; i < t; ++i) {
    // create count min sketch
    CMs[i] = new CMatrix_CSC(createCountMin(w, mu, Ps[i]->m));
    sk[i] = (*CMs[i] * *Ps[i]) * slice;// * W;
  }
  timer.stop();


  timer.start("Do recovering");
  // dealing with use_sketch part
  ptr = 0;
  for (auto i = use_sketch.begin(); i != use_sketch.end(); ++i) {
    CVector ind[MAX_LOGN]; // dyadic structure
    ind[t - 1].push_back(0);
    for (int j = t - 1; j >= 0; --j) {

      // convert csc col to int[]
      slicing(skCol, sk[j], ptr);
      // std::fill(skCol, skCol + w * mu, 0);
      // for (int k = sk[j].col_ptr[ptr]; k < sk[j].col_ptr[ptr + 1]; ++k)
      // 	skCol[sk[j].row[k]] += sk[j].val[k];
    

      // now do recover
      for (auto it = ind[j].begin(); it != ind[j].end(); ++it) {
	if (*it >= Ps[j]->m)
	  break;
	int v = recover(skCol, *CMs[j], *it);
	if (v > theta) {
	  if (j > 0) { // has not reached level 0  
	    for (int j0 = 0; j0 < step_size; ++j0)
	      ind[j - 1].push_back((*it) * step_size + j0);
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









CMatrix_CSC FastThreshMult_filter(const CMatrix_CSC &P, const CMatrix_CSC &Q, 
			   double theta, double rho, int w) {

  // set up the step_size to merge neighbor
  int step_size = STEP_SIZE;
  double threshold = (alpha + rho) / 4. * w * theta;

  CMatrix_CSC* Ps[MAX_LOGN]; // keep merged matrices
  CTimer timer;
  int t = 0; // # of levels in dyadic structure
  Ps[0] = new CMatrix_CSC(P);





  while (Ps[t]->m > 1) {
    Ps[t + 1] = new CMatrix_CSC(mergeNeighbor(*Ps[t], step_size));
    ++t;
  }
  ++t;



  // // calc L1 norm of P * Q
  // CVector&& R1 = calcL1Norm(P, Q);
  //  double threshold = (alpha + rho) / 4. * w * theta;  


  
  ///////////////////////////////////////////
  // recover heavy entries
  ///////////////////////////////////////////
  CVector val, row, col;
  int skCol[w * mu]; // keep extracted column vector
  CVector R1 = calcL1Norm(P, Q);
  
 

  


  /////////////////////////////////////////////////////
  /////////////  Using sketch /////////////////////////
  /////////////////////////////////////////////////////
  

  timer.start("Use sketch");
  // create sketch 
  CMatrix_CSC sk[MAX_LOGN];
  CMatrix_CSC* CMs[MAX_LOGN];


  timer.start("sketch matrix");
  for (int i = 0; i < t; ++i) {
    // create count min sketch
    CMs[i] = new CMatrix_CSC(createCountMin(w, mu, Ps[i]->m));
    sk[i] = (*CMs[i] * *Ps[i]) * Q;
  }
  timer.stop();


  timer.start("Do recovering");
  // dealing with use_sketch part

  for (int i = 0; i < Q.n; ++i) {
    
    CMatrix_CSC && Qi = Q[i];

    CVector ind[MAX_LOGN]; // dyadic structure
    ind[t - 1].push_back(0);
    for (int j = t - 1; j >= 0; --j) {

      // convert csc col to int[]
      slicing(skCol, sk[j], i);


      // now do recover
      for (auto it = ind[j].begin(); it != ind[j].end(); ++it) {
	if (*it >= Ps[j]->m)
	  break;
	int v = recover(skCol, *CMs[j], *it);

	if (v >= theta) {
	  if (j > 0) { // has not reached level 0  
	    for (int j0 = 0; j0 < step_size; ++j0)
	      ind[j - 1].push_back((*it) * step_size + j0);
	  }
	  else { // reach level 0, do verification 

	    if (v > theta && R1[i] > threshold) {
	      v = inner_product(Q[*it], Qi, theta);
	    }
	    if (v >= theta) {
	      val.push_back(v);
	      row.push_back(*it);
	      col.push_back(i);
	    }
	  } 
	} // if
      } // for (auto
    } // for (int j ..      

  }
  timer.stop();


  // release memory
  while (--t >= 0) {
    delete Ps[t];
    delete CMs[t];
  }

  timer.stop();
  
  return CMatrix_CSC(val.begin(), row.begin(), col.begin(), val.size(), P.m, Q.n);
}








CMatrix_CSC FastThreshMult_Simple(const CMatrix_CSC &P, const CMatrix_CSC &Q, 
				  double theta, double rho, int w) {

  // calc L1 norm of P * Q
  CVector&& R1 = calcL1Norm(P, Q);  
  double threshold = (alpha + rho)/ 4. * w * theta;
  // if (w < 1)
  //   w = optimalBucketsSize(R1, Q, rho, theta);

  CTimer timer;


  
  ///////////////////////////////////////////
  // recover heavy entries
  ///////////////////////////////////////////
  CVector val, row, col;
  int skCol[w * mu]; // keep extracted column vector


  /////////////////////////////////////////////////////
  /////////////  Using sketch /////////////////////////
  /////////////////////////////////////////////////////
  

  timer.start("Create sketches");
  // create sketch 
  timer.stop();

  timer.start("Skech the matrix");
  CMatrix_CSC CM(createCountMin(w, mu, P.m)); 
  CMatrix_CSC sk((CM * P) * Q);
  timer.stop();


  timer.start("Use sketch to recover");
  // dealing with use_sketch part

  for (int i = 0; i < Q.n; ++i) {

    // convert csc col to int[]
    slicing(skCol, sk, i);

    
    // now do recover
    for (int r = 0; r < P.m; ++r) {
      int v = recover(skCol, CM, r);
      if (v > theta && R1[i] > threshold) {
	v = inner_product(Q[r], Q[i], theta);
      }
      if (v > theta) {
	val.push_back(v);
	row.push_back(r);
	col.push_back(i);
      } 
    } // for (auto
  } // for (int i
    
  timer.stop();
  
  return CMatrix_CSC(val.begin(), row.begin(), col.begin(), val.size(), P.m, Q.n);
}
