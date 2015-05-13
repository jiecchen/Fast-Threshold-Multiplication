#include <iostream>
#include "Algebra.h"
#include <numeric>
const int M = 4;
const int N = 4;


int main() {
  int A[M][N];
  std::fill(A[0], A[0] + M * N, 1);
  A[0][0] = 1;
  A[1][2] = 2;
  A[1][3] = 10;
  A[2][2] = 1000;
  A[3][0] = -1;
  A[3][2] = 100;
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
  //  std::cerr << "nnz = " << nnz << std::endl;
  CMatrix_CSR mat(arr, arr + nnz, M, N);
  std::cout << mat << std::endl << std::endl;
  int ind[2];
  ind[0] = 0;
  ind[1] = 3;
  //  std::iota(ind.begin(), ind.end(), 1);
  CVector vsum = mat.sumRows(ind, ind + 2);
  std::cout << "sum of rows: " << vsum << std::endl;
  
  CVector vec(N);
  std::fill(vec.begin(), vec.end(), 1);
  vec[0] = 100;
  vec[3] = -2;
  std::cout << vec << std::endl;
  CVector res = mat * vec;
  std::cout << res << std::endl;

   

  for (int i = 0; i < nnz; ++i) 
    std::swap(arr[i].row, arr[i].col);
  
  CMatrix_COO B(arr, arr + nnz, N, M);
  B.sortByColumnRow();
  std::cout << B << std::endl;
  CMatrix_COO ans = thresh_mult_naive(mat, B, 10);
  std::cout << ans << std::endl;

  return 0;
}



