#include <iostream>
#include "Algebra.h"
#include "Sketch.h"
#include <cmath>
#include <numeric>
#include <iomanip>

// sum rows of an coo matrix
void sumRows_Coo(int *result, CMatrix_COO &P) {
  for (int i = 0; i < P.size(); ++i)
    result[P[i].col] += P[i].val;
}

int inner_prod(int *v, SparseVector &sv) {
  int res = 0;
  for (auto it = sv.begin(); it != sv.end(); ++it)
    res += v[it->ind] * it->val;
  return res;
}



// r = P * vec
// result = r[r > thresh]
// note:  P has to be sorted by RowColumn
void thresholdMult(SparseVector &result, CMatrix_COO &P, SparseVector &vec, double thresh) {
  int v[vec.size()];
  std::fill(v, v + vec.size(), 0);
  for (const VectorElement &it : vec)
    v[it.ind] = it.val;
  int t = 0; 
  while (t < P.size()) {
    int current_row = P[t].row;
    int prod = 0;
    prod += P[t].val * v[P[t].col];
    ++t;

    while (P[t].row == P[t-1].row) {
      prod += P[t].val * v[P[t].col];
      ++t;
    }

    if (prod > thresh) 
      result.push_back(VectorElement(current_row, prod));    
  }
}


// P, Q are binary matrices
// theta is a threshold > 0
// rho \in (0, 1) to control the accuracy
CMatrix_COO atLeastMult(CMatrix_COO &P, CMatrix_COO &Q, double theta, double rho) {
  CMatrix_COO thresh_R;
  int sumP[P.get_n()];
  sumRows_Coo(sumP, P);
  Q.sortByColumnRow();
  //assume R = PQ
  //now calculate ||R||_1
  int normR = 0;
  for (int i = 0; i < Q.size(); ++i)
    normR += sumP[Q[i].row] * Q[i].val;
  double eps = std::sqrt(rho * theta / (normR + 0.0));
  double K = 1. / eps;
  

  CDyadicSketch dyadic;

  dyadicSketch(dyadic, eps, 10, P);

  int t = 0;
  //TODO:
  //  + bug, when Q.size() == 0
  while (t < Q.size()) {
    SparseVector tmp;
    int current_col = Q[t].col;
    tmp.push_back(VectorElement(Q[t].row, Q[t].val));
    ++t;
    while (t < Q.size() && Q[t].col == Q[t - 1].col) {
      tmp.push_back(VectorElement(Q[t].row, Q[t].val));
      ++t;
    }
    
    double prod = inner_prod(sumP, tmp); // inner product
    P.sortByRowColumn();
    if (prod > K) { // use exact algorithm
      //TODO
      SparseVector result;
      thresholdMult(result, P, tmp, theta);
      for (auto it = result.begin(); it != result.end(); ++it)
	thresh_R.push_back(Element(it->ind, current_col, it->val));
    }
    else { // use dyadic structure to recover
      SparseVector result;
      thresholdRecover(result, dyadic, tmp.begin(), tmp.end(), theta);
      for (auto it = result.begin(); it != result.end(); ++it)
	thresh_R.push_back(Element(it->ind, current_col, it->val));
    }

  }

  return thresh_R;
  
}






int main() {
  // read P
  int M, N;
  std::cin >> M >> N;
  int row, col, val;
  std::vector<Element> arr;
  while (std::cin >> row >> col >> val) {
    arr.push_back(Element(row, col, val));
  }   

  CMatrix_COO B(arr.begin(), arr.end(), M, N);
  B.sortByColumnRow();
  std::cout << "B = \n" << B << std::endl;

  MatrixSketch sk;
  sketchMatrixCOO(sk, 0.05, 5, B);
  
  SparseVector vec;
  vec.push_back(VectorElement(0, 500));
  vec.push_back(VectorElement(1, 100));
  vec.push_back(VectorElement(2, 10));
  vec.push_back(VectorElement(3, 55));

  

  // try to construct a dyadic structure
  // and do recovering
  // CDyadicSketch dyadic;
  // dyadicSketch(dyadic, 0.1, 10, B);
  // int thresh = 201;
  // SparseVector result;
  // thresholdRecover(result, dyadic, vec.begin(), vec.end(), thresh);

  
  // for (unsigned int i = 0; i < result.size(); ++i)
  //   std::cout << result[i].ind << ", " << result[i].val << std::endl;
  // std::cout << std::endl;



  std::cout << B << std::endl;
  std::cout << vec << std::endl;
  SparseVector mulRes;
  B.sortByRowColumn();
  thresholdMult(mulRes, B, vec, 1);
  std::cout << mulRes << std::endl;

  
  // int *cm = sketchVector(sk, vec.begin(), vec.end());
  // for (int i  = 0; i < M; i++) 
  //   std::cout << recover(cm, sk, i) << std::endl;

  return 0;
}
















