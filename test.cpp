#include <iostream>
#include "Algebra.h"
#include "Sketch.h"
#include <cmath>
#include <numeric>
#include <iomanip>
#include <fstream>
#include "utils.h"
#include <ctime>


void test_utils_functions() {
  int M = 5;
  int N = 1000;
  CMatrix_COO coo(M, N);
  std::srand(time(NULL));
  for (int i = 0; i < M; ++i)
    for (int j = 0; j < N; ++j)
      if (std::rand() % 3 == 0)
	coo.push_back(Element(i, j, 1));

  CMatrix_CSC P(coo); 
  CMatrix_COO cooT = coo.T();
  CMatrix_CSC Q(cooT);
  //  std::cerr << "P = \n" << coo << std::endl;
  //  std::cerr << "Q = \n" << cooT << std::endl;
  std::cerr << "P * Q =\n" << toCoo(P * Q) << std::endl;
  CMatrix_COO coo_new(M, M);
  for (int i = 0; i < M; ++i)
    for (int j = 0; j < M; ++j)
      coo_new.push_back(Element(i, j, inner_product(Q[i], Q[j])));
  std::cerr << "P * Q =\n" << coo_new << std::endl;
}





int main() {
  test_utils_functions();

  return 0;
}
