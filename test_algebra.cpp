#include "Algebra.h"
#include <cstdlib>
#include <ctime>
#include <iostream>

int main() {
  int M = 5;
  int N = 10;

  CMatrix_COO coo(M, N);
  std::srand(time(NULL));
  for (int i = 0; i < M; ++i)
    for (int j = 0; j < N; ++j)
      if (std::rand() % 3 == 0)
	coo.push_back(Element(i, j, 0.5));

  CMatrix_CSC P(coo); 
  CMatrix_CSC&& Q = P.T();
  std::cerr << "P = \n" << coo << "\nQ = \n" << toCoo(Q) << std::endl;
  
  std::cerr << "P * Q =\n" << toCoo(P * Q) << std::endl;

  return 0;
}
