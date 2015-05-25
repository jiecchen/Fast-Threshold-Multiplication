#include <iostream>
#include "Algebra.h"
#include "Sketch.h"
#include <cmath>
#include <numeric>
#include <iomanip>



void test_mergeNeighbor(CMatrix_CSC &P) {
  std::cout << "P = \n" << toCoo(P) << std::endl;
  CMatrix_CSC Q = mergeNeighbor(P);
  std::cout << "mergeNeighbor(P) = \n" << toCoo(Q) << std::endl;
}


void test_csc(CMatrix_CSC &P, CMatrix_CSC &Q) {
  std::cout << "P = \n" << toCoo(P) << std::endl;
  std::cout << "Q = \n" << toCoo(Q) << std::endl;
  CMatrix_CSC mult = P * Q;
  std::cout << "P * Q = \n" << toCoo(mult) << std::endl;
}








int main() {
  const int MAX_NNZ = 1000;
  int Pm, Pn, nnz;
  std::cin >> Pm >> Pn;
  std::cin >> nnz;
  int row[MAX_NNZ];
  int col[MAX_NNZ];
  int val[MAX_NNZ];
  for (int i = 0; i < nnz; ++i) {
    std::cin >> row[i] >> col[i] >> val[i];
  }   
  CMatrix_CSC P(val, row, col, nnz, Pm, Pn);
  
  int Qm, Qn;
  std::cin >> Qm >> Qn;
  std::cin >> nnz;
  for (int i = 0; i < nnz; ++i) {
    std::cin >> row[i] >> col[i] >> val[i];
  }   
  CMatrix_CSC Q(val, row, col, nnz, Qm, Qn);


  test_csc(P, Q);
  test_mergeNeighbor(P);
  return 0;
}
















