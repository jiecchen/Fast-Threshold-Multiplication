#include <iostream>
#include "Algebra.h"
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
  std::cout << "B = \n" << B << std::endl;
  
  MatrixSketch sk;
  sketchMatrixCOO(sk, 0.05, 5, B);
  
  SparseVector vec;
  vec.push_back(VectorElement(1, 10));
  vec.push_back(VectorElement(2, 100));
  vec.push_back(VectorElement(3, 500));


  
  int *cm = sketchVector(sk, vec.begin(), vec.end());
  for (int i  = 0; i < M; i++) 
    std::cout << recover(cm, sk, i) << std::endl;

  return 0;
}
















