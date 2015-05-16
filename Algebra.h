#ifndef __ALGEBRA_H__
#define __ALGEBRA_H__
#include <vector>
#include <iostream>

//TODO:
//  + combine CSC & CSR to create CMatrix type?



////////////////////////////////////////////////////////////////////////////////
class Element {
public:
  Element(): row(0), col(0), val(0) {};
  Element(int row, int col, int val): row(row), col(col), val(val) {};
  int row, col, val;
};



////////////////////////////////////////////////////////////////////////////////
class CDenseMatrix {
public:
  CDenseMatrix(int m, int n);
  CDenseMatrix(int *data, int m, int n);
  ~CDenseMatrix() {
    //TODO: add it back, memory leak otherwise
    delete[] data;
  };
  int get_m() { return m; }
  int get_n() { return n; }
  int* get_row(int i) { return data + i * n; }
private:
  int m, n;
  int* data;
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
//typedef std::vector<Element> CMatrix_COO;






////////////////////////////////////////////////////////////////////////////////
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




////////////////////////////////////////////////////////////////////////////////
// sparse matrix using CSR format
class CMatrix_CSR {
public:
  // arr should be sorted, from left-to-right, top-to-bottom
  CMatrix_CSR(Element *arr_s, Element *arr_e, int _m,  int _n);

  // sparse matrix *  dense vector
  friend CVector operator *(const CMatrix_CSR &mat,  CVector &vec);
  friend CDenseMatrix operator *(const CMatrix_CSR &mat,  CMatrix_COO &coo);
  CMatrix_COO toCOO() const;
  CVector sumRows(Index_iter ibegin, Index_iter iend); // given indices of rows, return the sum
  CVector sumRows(int *s, int *e);
  friend std::ostream& operator << (std::ostream &os, const CMatrix_CSR &mat);  
  // create list of dyatic count sketches
  //  friend DyadicCountMin createDyadicCountMin(double eps, double mu, CMatrix_CSR &P);
  
  ~CMatrix_CSR() {
    delete[] col_val;
    delete[] row_ptr;
  }

  int get_m() { return m; }
  int get_n() { return n; }
private:
  int m, n;
  int nnz;
  VectorElement *col_val;
  int *row_ptr;
};



////////////////////////////////////////////////////////////////////////////////
struct CountMinSketch {
  int recover(int cood);
  ~CountMinSketch() {
    delete hash;
    delete[] cm;
  }
  int w, u;
  CMatrix_COO* hash;
  int *cm;
};




////////////////////////////////////////////////////////////////////////////////
// naive threshhold-multiplication for matrix_csr * matrix_coo 
CMatrix_COO thresh_mult_naive(CMatrix_CSR &A, CMatrix_COO &coo, double thresh);

CountMinSketch createCountMin(CountMinSketch &sk, double eps, int u, CMatrix_COO &coo);


// override << for CVector
std::ostream& operator << (std::ostream & os, CVector &vec);
// override << for CDenseMatrix
std::ostream& operator << (std::ostream & os, CDenseMatrix &mat);


CDenseMatrix operator *(CMatrix_COO &A, CMatrix_COO &B);




#endif


















