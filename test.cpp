#include <iostream>
#include "Algebra.h"
const int N = 4;

int main() {
  int A[N][N];
  A[1][2] = 1;
  A[0][3] = 2;
  Element arr[2];

  for (int i = 0; i < N; ++i)
    for (int j = 0; j < N; ++j) 
      if (A[i][j] != 0) {
	arr[i].row = i;
	arr[i].col = j;
	arr[i].val = A[i][j];
      }
  Matrix_CSR mat(arr, 2, N, N);
  std::cout << mat << std::endl;
  return 0;
}

