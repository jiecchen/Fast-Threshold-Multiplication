#ifndef __ALGEBRA_H__
#define __ALGEBRA_H__
#include <vector>
#include <iostream>

//TODO:
//  + combine CSC & CSR to create CMatrix type?

int const INFINITY = 10000000;


////////////////////////////////////////////////////////////////////////////////
class Element {
public:
  Element(): row(0), col(0), val(0) {};
  Element(int row, int col, int val): row(row), col(col), val(val) {};
  int row, col, val;
};





////////////////////////////////////////////////////////////////////////////////
// can be used to represent sparse vector
class VectorElement {
public:
  VectorElement(): ind(0), val(0) {};
  VectorElement(int ind, int val): ind(ind), val(val) {};
  int ind, val;
};

typedef std::vector<int> CVector; // dense vector
typedef CVector::iterator CVector_iter;
typedef std::vector<int>::iterator Index_iter;
typedef std::vector<VectorElement> SparseVector;
typedef std::vector<VectorElement>::iterator SparseVector_iter;





///////////////////////////////////////////////////////////////////////////////
////////////////////////// CMatrix_COO ////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

class CMatrix_COO {
public:
  CMatrix_COO(): m(0), n(0) {};
  CMatrix_COO(int m, int n): m(m), n(n) {};
  CMatrix_COO(Element *arr_s, Element *arr_e, int m, int n): m(m), n(n) {
    data.assign(arr_s, arr_e);
  };
  void sortByColumnRow();
  void sortByRowColumn();
  int size() { return data.size(); }
  void push_back(Element e) { 
    data.push_back(e); 
    if (e.row > m)
      m = e.row;
    if (e.col > n)
      n = e.col;
  }
  Element& operator[] (int i) { return data[i]; }
  friend std::ostream& operator << (std::ostream & os, CMatrix_COO &coo);
  
  int get_m() { return m; }
  int get_n() { return n; }
  // return transpose
  CMatrix_COO T() {
    CMatrix_COO cm = *this;
    std::swap(cm.m, cm.n);
    for (auto it = cm.data.begin(); it != cm.data.end(); ++it)
      std::swap(it->row, it->col);
    return cm;
  };
private:
  int m, n;
  std::vector<Element> data;
};

// multiplication for coo format
// return an int* as a result
// complexity O(nnz(B) * A.m)
int* coo_mult(CMatrix_COO &A, CMatrix_COO &B);



/////////////////////////////////////////////////////////////////////////
/////////////////////// Sketch Related //////////////////////////////////
/////////////////////////////////////////////////////////////////////////
class MatrixSketch {
public:
  ~MatrixSketch() {
    delete[] data;
    delete hash;
  };
  int w, u; // h = w * u
  int n;
  CMatrix_COO *hash;
  int *data; // will be h * n
};

// MatrixSketch * Vector -> CountMinSketch
// sketching a coo matrix
// keep result in MatrixSketch
void sketchMatrixCOO(MatrixSketch &sk, double eps, int u, CMatrix_COO &B);
// return int* as count min sketch of the vector specfied by it_s -> it_e
// remember to release the memory
int* sketchVector(MatrixSketch &sk, SparseVector_iter it_s, SparseVector_iter it_e);
// recover the coordinate given sketch and hash
int recover(int *cm, const MatrixSketch &sk, int coor);

#endif


















