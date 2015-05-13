#ifndef __ALGEBRA_H__
#define __ALGEBRA_H__
#include <vector>
#include <iostream>

//TODO:
//  + combine CSC & CSR to create CMatrix type?



class Element {
public:
  Element(): row(0), col(0), val(0) {};
  Element(int row, int col, int val): row(row), col(col), val(val) {};
  int row, col, val;
};


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

class CMatrix_COO {
public:
  CMatrix_COO(): m(0), n(0) {};
  CMatrix_COO(int m, int n): m(m), n(n) {};
  CMatrix_COO(Element *arr, int length, int m, int n): m(m), n(n) {
    data.assign(arr, arr + length);
  }
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
  
private:
  int m, n;
  std::vector<Element> data;
};

// sparse matrix using CSR format
class CMatrix_CSR {
public:
  // arr should be sorted, from left-to-right, top-to-bottom
  CMatrix_CSR(Element *arr, int _nnz, int _m,  int _n);

  // sparse matrix *  dense vector
  CVector operator *(const CVector &vec);
  CMatrix_COO toCOO() const;
  CVector sumRows(Index_iter ibegin, Index_iter iend); // given indices of rows, return the sum
  CVector sumRows(int *s, int *e);
  friend std::ostream& operator << (std::ostream &os, const CMatrix_CSR &mat);  

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

// naive threshhold-multiplication for matrix_csr * matrix_coo 
CMatrix_COO thresh_mult_naive(CMatrix_CSR &A, CMatrix_COO &coo, double thresh);



// override << for CVector
std::ostream& operator << (std::ostream & os, CVector &vec);





#endif

















