// create experimental data
#include <iostream>
#include <cstdlib>
#include "Algebra.h"
#include <vector>
using namespace std;

int main() {
  int Pm, Pn, Qn;
  int nnz;
  cerr << "P_m P_n Q_n = ";
  cin >> Pm >> Pn >> Qn;
  cerr << "nnz(P) = ";
  cin >> nnz;
  int P[Pm][Pn];
  vector<Element> nonz;
  fill(P[0], P[0] + Pm * Pn, 0);
  for (int i = 0; i < nnz; ++i)
    P[rand() % Pm][rand() % Pn] = 1;
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
    Q[rand() % Pn][rand() % Qn] = 1;
  for (int i = 0; i < Pn; ++i)
    for (int j = 0; j < Qn; ++j)
      if (Q[i][j] > 0)
	nonz.push_back(Element(i, j, 1));
  cout << Pn << " " << Qn << endl;
  cout << nonz.size() << endl;
  for (const Element &e : nonz)
    cout << e.row << " " << e.col << " " << e.val << endl;

  return 0;  
}










