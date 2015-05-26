#include "Algebra.h"
#include "Sketch.h"
#include <numeric>
//#include <iomanip>
#include <cmath>
#include <cstdlib>
#include <ctime>

//Principle:
//   + each function only do one thing
//   + keep simple, keep stupid



/////////////////////////////////////////////////////////////////////////
/////////////////////// Sketch Related //////////////////////////////////
/////////////////////////////////////////////////////////////////////////



// sketching a coo matrix
// keep result in MatrixSketch
void sketchMatrixCOO(MatrixSketch &sk, double eps, int u, CMatrix_COO &B) {
  std::srand(time(NULL));
  sk.w = std::ceil(1. / eps);
  sk.u = u;
  int h = sk.w * u;
  sk.n = B.get_n();
  // create count min sketch
  CMatrix_COO* cm = new CMatrix_COO(h, B.get_m());
  for (int i = 0; i < B.get_m(); ++i) 
    for (int j = 0; j < u; ++j) 
      cm->push_back(Element(std::rand() % sk.w + j * sk.w, i, 1));

  sk.hash = cm;
  sk.data = coo_mult(*cm, B);

  sk.hash->sortByColumnRow();
}


int* sketchVector(MatrixSketch &sk, SparseVector_iter it_s, SparseVector_iter it_e) {
  int h = sk.w * sk.u;
  // init as 0
  int *cm = new int[h]();
  for (int i = 0; i < h; i++)
    for (auto it = it_s; it != it_e; ++it)
      cm[i] += sk.data[i * sk.n + it->ind] * it->val;
  return cm;
}


int recover(int *cm, const MatrixSketch &sk, int coor) {
  int m = INFINITY;
  CMatrix_COO &hash = *(sk.hash);
  if (coor * sk.u >= hash.size())
    return 0;
  for (int i = 0; i < sk.u; ++i) {
    int r = hash[coor * sk.u + i].row;
    m = std::min(cm[r], m);
  }
  return m;
}

// recover the entry rw
int recover(int sk[], const CMatrix_CSC &cm, int rw) {
  int m = INFINITY;
  for (int i = 0; i < mu; ++i) {
    int r = cm.row[rw * mu + i];
    m = std::min(sk[r], m);
  }
  //  std::cerr << "m = " << m << std::endl;
  return m;
}

// let oldCoo = [R_0; R_1; R_2; R_3; ...]
// newCoo = [R_0 + R_1; R_2 + R_3; ...]
// oldCoo has to be sorted by ColumnRow
// newCoo is empty
void mergeNeighbor(CMatrix_COO &newCoo, CMatrix_COO &oldCoo) {
  newCoo.set_mn((oldCoo.get_m() + 1) / 2, oldCoo.get_n());
  if (oldCoo.size() <= 0)
    return;
  newCoo.push_back(Element(oldCoo[0].row / 2, oldCoo[0].col, oldCoo[0].val));
  for (int i = 1; i < oldCoo.size(); ++i)
    if (oldCoo[i].col == newCoo.back().col &&
	oldCoo[i].row / 2 == newCoo.back().row) {
      newCoo.back().val += oldCoo[i].val;
    }
    else
      newCoo.push_back(Element(oldCoo[i].row / 2, oldCoo[i].col, oldCoo[i].val));
}



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



double calcL1Norm(const CMatrix_CSC& P, const CMatrix_CSC& Q, const CMatrix_CSC& W) {
  int val[P.m];
  std::fill(val, val + P.m, 1);
  int row[P.m];
  std::fill(row, row + P.m, 0);
  int col[P.m];
  std::iota(col, col + P.m, 0);
  CMatrix_CSC S(val, row, col, P.m, 1, P.m);
  CMatrix_CSC&& res = ((S * P) * Q) * W;
  int r = 10;
  for (int i = 0; i < W.n; ++i)
    if (i >= res.nnz) {
      break;
    }
    else
      r = std::max(r, res.val[i]);
  std::cout << "r = " << r << std::endl;
  return r;
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
			    double theta, double rho) {

  CMatrix_CSC* Ps[MAX_LOGN];
  
  int t = 0;
  Ps[0] = new CMatrix_CSC(P);
  //  std::cerr << "Ps[0] = \n" << toCoo(*Ps[0]) << std::endl;
  while (Ps[t]->m > 1) {
    Ps[t + 1] = new CMatrix_CSC(mergeNeighbor(*Ps[t]));
    ++t;
  }
  ++t;

  // calc ||P * Q  * W||_1
  double R1 = calcL1Norm(P, Q, W);
  int w = std::ceil(R1 / theta / rho);
  std::cerr << "w = " << w << " mu = " << mu << std::endl;

  // create sketches
  CMatrix_CSC* CMs[MAX_LOGN];
  CMatrix_CSC sk[MAX_LOGN];
  for (int i = 0; i < t; ++i) {
    // create count min sketch
    CMs[i] = new CMatrix_CSC(createCountMin(w, mu, Ps[i]->m));
    sk[i] = ((*CMs[i] * *Ps[i]) * Q) * W;

  }
  

  
  // recover heavy entries
  CVector val, row, col;
  int skCol[w * mu]; // keep extracted column vector
  for (int i = 0; i < W.n; ++i) { // for column
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


  // release memory
  while (--t >= 0) {
    delete Ps[t];
  }

  return CMatrix_CSC(val.begin(), row.begin(), col.begin(), val.size(), P.m, W.n);
}

// // recover entries > theta in P * Q * W
// CMatrix_CSC FastThreshMult(CMatrix_CSC &P, CMatrix_CSC &Q, CMatrix_CSC &W, double theta) {
// }



// create O(log m) matrix sketches, form as
// dyadic sketch
void dyadicSketch(CDyadicSketch &dyadic,  double eps, int u, CMatrix_COO &B) {
  dyadic.clear();
  CMatrix_COO coo[MAX_LOGN];
  coo[0] = B;
  coo[0].sortByColumnRow();
  int t = 0;
  //TODO:
  //  + bugs, if coo[t].size() == 0
  while (coo[t].get_m() > 1) {
    mergeNeighbor(coo[t + 1], coo[t]);
    t++;
  }
  while (t >= 0) {
    MatrixSketch *sk = new MatrixSketch();
    sketchMatrixCOO(*sk, eps, u, coo[t--]);
    dyadic.push_back(sk);
  }
}




// given threshhold, recover all
// entries in dyadic(vec_s, vec_e)
// approximately recover entries in PQ that >= thresh
// keep in result
void thresholdRecover(SparseVector &result, CDyadicSketch &dyadic, 
		      SparseVector_iter vec_s, SparseVector_iter vec_e, double thresh) {
  std::vector<int> ind[MAX_LOGN];
  ind[0].push_back(0);
  for (int i = 0; i < dyadic.size(); ++i) {
    int *cm = sketchVector(*dyadic[i], vec_s, vec_e);
    for (auto it = ind[i].begin(); it != ind[i].end(); ++it) {
      int v = recover(cm, *dyadic[i], *it); 
      if (v > thresh) {
	if (i + 1 < dyadic.size()) {
	  ind[i + 1].push_back((*it) * 2);
	  ind[i + 1].push_back((*it) * 2 + 1);
	}
	else
	  result.push_back(VectorElement(*it, v));	
      }
    }
    delete[] cm; 
  }
 
}







// P, Q are boolean matrices
// theta is a threshold > 0
// rho \in (0, 1) to control the accuracy
CMatrix_COO atLeastMult(CMatrix_COO &P, CMatrix_COO &Q, double theta, double rho) {
  CMatrix_COO thresh_R(P.get_m(), Q.get_n());

  P.sortByRowColumn();
  int sumP[P.get_n()];
  sumRows_Coo(sumP, P);



  //assume R = PQ
  //now calculate ||R||_1
  int normR = 0;
  for (int i = 0; i < Q.size(); ++i)
    normR += sumP[Q[i].row] * Q[i].val;
  std::cerr << "||R|| = " << normR << std::endl;

  double W = std::sqrt(normR * theta * rho);
  double eps = rho * theta / W;
  std::cerr << "eps = " << eps << " W = " << W << std::endl;

  CDyadicSketch dyadic;

  std::cerr << "Creating dyadic structure ..." << std::endl;
  dyadicSketch(dyadic, eps, 20, P);
  std::cerr << "Creating dyadic structure DONE" << std::endl;
  P.sortByRowColumn();
  Q.sortByColumnRow();
  int t = 0;
  while (t < Q.size()) {
    SparseVector tmp;
    int current_col = Q[t].col;

    // extract one column from coo format
    // put the result to tmp
    tmp.push_back(VectorElement(Q[t].row, Q[t].val));
    ++t;
    while (t < Q.size() && Q[t].col == Q[t - 1].col) {
      tmp.push_back(VectorElement(Q[t].row, Q[t].val));
      ++t;
    }
    
    
    //DEBUG
    //  W = 0;

    // calculate inner product sumP * tmp
    double prod = inner_prod(sumP, tmp); 
    if (prod > W) { // use exact algorithm
      //TODO
      std::cerr << "use exact algorithm" << std::endl;
      SparseVector result;
      thresholdMult(result, P, tmp, theta);
      for (auto it = result.begin(); it != result.end(); ++it)
	thresh_R.push_back(Element(it->ind, current_col, it->val));
    }
    else { // use dyadic structure to recover
      std::cerr << "use dyadic structure" << std::endl;
      SparseVector result;
      thresholdRecover(result, dyadic, tmp.begin(), tmp.end(), theta);
      for (auto it = result.begin(); it != result.end(); ++it)
	thresh_R.push_back(Element(it->ind, current_col, it->val));
    }

  }

  return thresh_R;
  
}










