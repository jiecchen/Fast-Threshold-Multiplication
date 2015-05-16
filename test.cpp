#include <iostream>
#include "Algebra.h"
#include <numeric>
#include <iomanip>
const int M = 4;
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
   

  CMatrix_COO B(arr, arr + nnz, N, M);
  std::cout << "B = \n" << B << std::endl;
  

  return 0;
}












