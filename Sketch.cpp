#include "Algebra.h"
#include "Sketch.h"
#include <cmath>

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










