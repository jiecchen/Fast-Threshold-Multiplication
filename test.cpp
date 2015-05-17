#include <iostream>
#include "Algebra.h"
#include "Sketch.h"
#include <numeric>
#include <iomanip>
const int M = 10;
const int N = 4;


int main() {
  int A[M][N];
  std::fill(A[0], A[0] + M * N, 0);
  A[0][0] = 21;
  A[1][2] = 2;
  A[1][3] = 12;
  A[2][2] = 12;
  A[3][0] = 2;
  A[3][2] = 2;
  A[7][1] = 100;
  Element arr[M * N];
  int nnz = 0;
  for (int i = 0; i < M; ++i)
    for (int j = 0; j < N; ++j) 
      if (A[i][j] != 0) {
	arr[nnz].row = i;
	arr[nnz].col = j;
	arr[nnz].val = A[i][j];
	++nnz;
      }
   

  CMatrix_COO B(arr, arr + nnz, M, N);
  B.sortByColumnRow();
  std::cout << "B = \n" << B << std::endl;

  CMatrix_COO newCoo;
  mergeNeighbor(newCoo, B);
  std::cout << "newCoo =\n" << newCoo << std::endl;

  CMatrix_COO newCoo1;
  mergeNeighbor(newCoo1, newCoo);
  std::cout << "newCoo1 =\n" << newCoo1 << std::endl;



  MatrixSketch sk;
  sketchMatrixCOO(sk, 0.05, 5, B);
  
  SparseVector vec;
  vec.push_back(VectorElement(1, 10));
  vec.push_back(VectorElement(2, 100));
  vec.push_back(VectorElement(3, 500));


  // try to construct a dyadic structure
  // and do recovering
  CDyadicSketch dyadic;
  dyadicSketch(dyadic, 0.1, 10, B);
  int thresh = 201;
  SparseVector result;
  thresholdRecover(result, dyadic, vec.begin(), vec.end(), thresh);

  // std::vector<int> result(B.get_m());
  // for (int i = 0; i < dyadic.size(); ++i) {
  //   int *cm = sketchVector(*dyadic[i], vec.begin(), vec.end());
  //   for (auto it = ind[i].begin(); it != ind[i].end(); ++it) {
  //     int v = recover(cm, *dyadic[i], *it); 
  //     if (v > thresh) {
  // 	ind[i + 1].push_back((*it) * 2);
  // 	ind[i + 1].push_back((*it) * 2 + 1);
  // 	if (i + 1 == dyadic.size()) { // last level
  // 	  result[*it] = v;
  // 	}
  //     }
  //   }
  //   delete[] cm; 
  // }
  
  for (unsigned int i = 0; i < result.size(); ++i)
    std::cout << result[i].ind << ", " << result[i].val << std::endl;
  std::cout << std::endl;
  
  // int *cm = sketchVector(sk, vec.begin(), vec.end());
  // for (int i  = 0; i < M; i++) 
  //   std::cout << recover(cm, sk, i) << std::endl;

  return 0;
}
















