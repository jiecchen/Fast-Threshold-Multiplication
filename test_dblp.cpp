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
#include "Algebra.h"
#include "utils.h"


const int MAX_N = 425876; // limit for # of nodes
const int MAX_M = 3000000; // limit for # of edges

int main(int argc, char **argv) {
  if (argc != 2) {
    std::cerr << "Check your input!" << std::endl;
    std::exit(1);
  }
  
  CTimer timer;
  timer.start();
  CMatrix_COO coo(MAX_N, MAX_N);
  std::ifstream fin(argv[1]);
  int row, col;
  while (fin >> row >> col) {
    coo.push_back(Element(row, col, 1));
    coo.push_back(Element(col, row, 1));
  }
  fin.close();
  timer.stop("Read data ");
  timer.start();
  CMatrix_CSC P(coo, false);
  timer.stop("Create CSC Matrix ");
  std::cout << "Matrix: " << P.m << "x" << P.n << " nnz = " << P.nnz << std::endl;

  timer.start();
  CMatrix_CSC&& res = P * P[1];
  timer.stop("Matrix-Vector Product ");
  std::cout << "Matrix: " << res.m << "x" << res.n << " nnz = " << res.nnz << std::endl;
  CMatrix_COO res_coo = toCoo(res);
  for (int i = 0; i < res_coo.size(); ++i)
    if (res_coo[i].val > 20)
      std::cout << res_coo[i].row << ", " << res_coo[i].val << std::endl;
  return 0;
}
















