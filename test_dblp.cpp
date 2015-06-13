
/*
  Run test for dblp data.

  # Undirected graph: ../../data/output/dblp.ungraph.txt
  # DBLP
  # Nodes: 317080 Edges: 1049866
  # FromNodeId    ToNodeId
*/
#include <cstdlib>
#include <iostream>
#include <fstream>
#include "Sketch.h"
#include "Algebra.h"
#include "utils.h"


int main(int argc, char **argv) {
  if (argc != 2) {
    std::cerr << "Check your input!" << std::endl;
    std::exit(1);
  }
  CTimer timer;
  timer.start("Read data");
  CMatrix_COO coo;
  std::ifstream fin(argv[1]);
  int M = 0;
  int N = 0;
  int row, col;
  while (fin >> row >> col) {
    M = std::max(M,row);
    N = std::max(N, col);
    coo.push_back(Element(row, col, 1));
    //    coo.push_back(Element(col, row, 1));
  }
  ++N;
  ++M;
  coo.set_mn(M, N);
  fin.close();
  timer.stop();
  timer.start("Create CSC Matrix ");
  CMatrix_CSC P(coo, false);
  CMatrix_COO cooT(coo.T());
  CMatrix_CSC Q(cooT, false);
  timer.stop();
  std::cerr << "Matrix: " << P.m << "x" << P.n << " nnz = " << P.nnz << std::endl;
  
  double theta= 40;
  double rho = 1.1;

  int w = 300;


  // std::cerr << std::endl;
  // timer.start("Naive P * Q ");
  // CMatrix_COO&& res = toCoo(P * Q);
  // timer.stop();
  // res.print(theta);


  
  
  // std::cerr << std::endl;
  // timer.start("Use FastThreshMult_Simple");
  // CMatrix_COO &&new_res = toCoo(FastThreshMult_Simple(P, Q, theta, rho, w));
  // timer.stop();
  // new_res.print();

 

  std::cerr << std::endl;
  timer.start("Use FastThreshMult_filter");
  CMatrix_COO &&new_res1 = toCoo(FastThreshMult_filter(P, Q, theta, rho, w));
  timer.stop();
  new_res1.print();
 
  return 0;
}
















