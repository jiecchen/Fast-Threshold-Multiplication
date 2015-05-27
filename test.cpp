#include <iostream>
#include "Algebra.h"
#include "Sketch.h"
#include <cmath>
#include <numeric>
#include <iomanip>
#include <fstream>
#include "utils.h"

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



const int MAX_NNZ = 2000000;
int row[MAX_NNZ];
int col[MAX_NNZ];
int val[MAX_NNZ];




int main() {
  int Pm, Pn, nnz;
  CTimer timer;
  timer.start();
  std::cin >> Pm >> Pn;
  std::cin >> nnz;
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

  int Wm, Wn;
  std::cin >> Wm >> Wn;
  std::cin >> nnz;
  for (int i = 0; i < nnz; ++i) {
    std::cin >> row[i] >> col[i] >> val[i];
  }   
  CMatrix_CSC W(val, row, col, nnz, Wm, Wn);
  timer.stop("Read P Q W ");
  //  test_csc(P, Q);
  //  test_mergeNeighbor(P);

  // std::cout << "P * Q = \n" << toCoo(P * Q) << std::endl 
  // 	    << "Q * W = \n" << toCoo(Q * W) << std::endl;
  timer.start();
  // CMatrix_CSC ans = P * Q * W;
  timer.stop("Calculate P * Q * W ");
  //std::cout << "P * Q * W =\n" << toCoo(ans) << std::endl;
  timer.start();
  CMatrix_CSC A = FastThreshMult(P, Q, W, 3, 200000, 0.5);
  
  timer.stop("FastThreshMult ");
  std::cout << toCoo(A) << std::endl;
  return 0;

}
















