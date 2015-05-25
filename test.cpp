#include <iostream>
#include "Algebra.h"
#include "Sketch.h"
#include <cmath>
#include <numeric>
#include <iomanip>






void test_csc() {
  // read matrix P
  int Pm, Pn, nnz;
  std::cin >> Pm >> Pn;
  std::cin >> nnz;
  int row[nnz];
  int col[nnz];
  int val[nnz];
  std::vector<Element> arr;
  for (int i = 0; i < nnz; ++i) {
    std::cin >> row[i] >> col[i] >> val[i];
    arr.push_back(Element(row[i], col[i], val[i]));
  }   
  CMatrix_COO coo(arr.begin(), arr.end(), Pm, Pn);
  CMatrix_CSC csc(val, row, col, nnz, Pm, Pn);
  
  // int v[] = {0, 1, 2, 3, 4, 5};
  
  // SparseVector result;
  // thresh_mult(result, csc, v, v + 6, 0);
  std::cout << "mat = \n" << coo << std::endl;
  //  std::cout << result << std::endl;
    


  // read matrix Q
  int Qm, Qn;
  std::cin >> Qm >> Qn;
  std::cin >> nnz;
  arr.clear();
  for (int i = 0; i < nnz; ++i) {
    std::cin >> row[i] >> col[i] >> val[i];
    arr.push_back(Element(row[i], col[i], val[i]));
  }   
  CMatrix_COO Q(arr.begin(), arr.end(), Qm, Qn);
  CMatrix_CSC new_csc(val, row, col, nnz, Qm, Qn);
  std::cout << "new mat = \n" << Q << std::endl;
  CMatrix_CSC qq = csc * new_csc;
  Q = toCoo(qq);
  std::cout << "mul = \n" << Q << std::endl;

}









int main() {
  test_csc();
  // // read matrix P
  // int Pm, Pn, Qm, Qn, nnz;
  // std::cin >> Pm >> Pn;
  // std::cin >> nnz;
  // int row, col, val;
  // std::vector<Element> arr;
  // for (int i = 0; i < nnz; ++i) {
  //   std::cin >> row >> col >> val;
  //   arr.push_back(Element(row, col, val));
  // }   
  // CMatrix_COO P(arr.begin(), arr.end(), Pm, Pn);
  
  // // read matrix Q
  // std::cin >> Qm >> Qn;
  // std::cin >> nnz;
  // arr.clear();
  // for (int i = 0; i < nnz; ++i) {
  //   std::cin >> row >> col >> val;
  //   arr.push_back(Element(row, col, val));
  // }   
  // CMatrix_COO Q(arr.begin(), arr.end(), Qm, Qn);
  
 

  // std::cout << "P = \n" << P << std::endl;
  // std::cout << "Q = \n" << Q << std::endl;

  // double theta = 2;
  // double rho = 0.1;
  // CMatrix_COO result = atLeastMult(P, Q, theta, rho);
  // std::cout << "PQ = \n" << result << std::endl;
  // result.print();

  

  // // try to construct a dyadic structure
  // // and do recovering
  // // CDyadicSketch dyadic;
  // // dyadicSketch(dyadic, 0.1, 10, B);
  // // int thresh = 201;
  // // SparseVector result;
  // // thresholdRecover(result, dyadic, vec.begin(), vec.end(), thresh);

  
  // // for (unsigned int i = 0; i < result.size(); ++i)
  // //   std::cout << result[i].ind << ", " << result[i].val << std::endl;
  // // std::cout << std::endl;



  // // std::cout << B << std::endl;
  // // std::cout << vec << std::endl;
  // // SparseVector mulRes;
  // // B.sortByRowColumn();
  // // thresholdMult(mulRes, B, vec, 1);
  // // std::cout << mulRes << std::endl;

  
  // // int *cm = sketchVector(sk, vec.begin(), vec.end());
  // // for (int i  = 0; i < M; i++) 
  // //   std::cout << recover(cm, sk, i) << std::endl;

  return 0;
}
















