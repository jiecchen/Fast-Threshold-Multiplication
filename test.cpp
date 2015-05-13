#include <iostream>
#include "Algebra.h"
const int M = 3;
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
  CMatrix_CSR mat(arr, nnz, M, N);
  std::cout << mat << std::endl << std::endl;


  CVector vec(N);
  std::fill(vec.begin(), vec.end(), 1);
  vec[0] = 100;
  vec[3] = -2;
  std::cout << vec << std::endl;
  CVector res = mat * vec;
  std::cout << res << std::endl;

  CMatrix_COO coo = mat.toCOO();
  coo.sortByRowColumn();
  std::cout << coo << std::endl;
  coo.sortByColumnRow();
  std::cout << coo << std::endl;
  
  


  return 0;
}



