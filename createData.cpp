// create experimental data
#include <iostream>
#include <cstdlib>
#include "Algebra.h"
#include <vector>
#include <random>
#include <ctime>
using namespace std;

int main() {
  
  default_random_engine g(time(NULL));
  discrete_distribution<int> dist {39,18,11,11,2,2,1,1};


  int Pm, Pn, Qn, Wn;
  int nnz;
  cerr << "P_m P_n Q_n Wn = ";
  cin >> Pm >> Pn >> Qn >> Wn;
  cerr << "nnz(P) = ";
  cin >> nnz;
  int P[Pm][Pn];
  vector<Element> nonz;
  fill(P[0], P[0] + Pm * Pn, 0);
  for (int i = 0; i < nnz; ++i)
    P[dist(g) % Pm][dist(g) % Pn] = 1;
  for (int j = 0; j < Pn; ++j)
    for (int i = 0; i < Pm; ++i)
      if (P[i][j] > 0)
	nonz.push_back(Element(i, j, 1));
  cout << Pm << " " << Pn << endl;
  cout << nonz.size() << endl;
  for (const Element &e : nonz)
    cout << e.row << " " << e.col << " " << e.val << endl;

  nonz.clear();
  cerr << "nnz(Q) = ";
  cin >> nnz;
  int Q[Pn][Qn];
  fill(Q[0], Q[0] + Pn * Qn, 0);
  for (int i = 0; i < nnz; ++i)
    Q[dist(g) % Pn][dist(g) % Qn] = 1;
  for (int j = 0; j < Qn; ++j)
    for (int i = 0; i < Pn; ++i)
      if (Q[i][j] > 0)
	nonz.push_back(Element(i, j, 1));
  cout << Pn << " " << Qn << endl;
  cout << nonz.size() << endl;
  for (const Element &e : nonz)
    cout << e.row << " " << e.col << " " << e.val << endl;


  nonz.clear();
  cerr << "nnz(W) = ";
  cin >> nnz;
  int W[Qn][Wn];
  fill(W[0], W[0] +  Qn * Wn, 0);
  for (int i = 0; i < nnz; ++i)
    W[dist(g) % Qn][dist(g) % Wn] = 1;
  for (int j = 0; j < Wn; ++j)
    for (int i = 0; i < Qn; ++i)
      if (W[i][j] > 0)
	nonz.push_back(Element(i, j, 1));
  cout << Qn << " " << Wn << endl;
  cout << nonz.size() << endl;
  for (const Element &e : nonz)
    cout << e.row << " " << e.col << " " << e.val << endl;



  return 0;  
}










