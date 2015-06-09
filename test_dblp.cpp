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



const int MAX_N = 425876; // limit for # of nodes
const int MAX_M = 3000000; // limit for # of edges

int main(int argc, char **argv) {
  if (argc != 2) {
    std::cerr << "Check your input!" << std::endl;
    std::exit(1);
  }
  int N = 34580;
  int M = 34540;
  CTimer timer;
  timer.start();
  CMatrix_COO coo(N, M);
  std::ifstream fin(argv[1]);
  int row, col;
  while (fin >> row >> col) {
    coo.push_back(Element(row, col, 1));
    //    coo.push_back(Element(col, row, 1));
  }
  fin.close();
  timer.stop("Read data ");
  timer.start();
  CMatrix_CSC P(coo, false);
  CMatrix_COO cooT(coo.T());
  CMatrix_CSC Q(cooT, false);
  timer.stop("Create CSC Matrix ");
  std::cout << "Matrix: " << P.m << "x" << P.n << " nnz = " << P.nnz << std::endl;
  
  double theta= 300;

  // timer.start();
  // CMatrix_COO&& res = toCoo(P * Q);
  // res.print(theta);
  // timer.stop("Naive P * P ");




  timer.start();
  CMatrix_COO&& new_res = toCoo(FastThreshMult(P, Q, 10, theta, 0.3));
  timer.stop("Use Our Algorithm ");
  new_res.print();
  return 0;
}
















