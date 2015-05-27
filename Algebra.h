#ifndef __ALGEBRA_H__
#define __ALGEBRA_H__
#include <iterator>
#include <vector>
#include <iostream>

//TODO:
//  + combine CSC & CSR to create CMatrix type?

int const INFINITY = 100000000;


////////////////////////////////////////////////////////////////////////////////
class Element {
public:
  Element(): row(0), col(0), val(0) {};
  Element(int row, int col, int val): row(row), col(col), val(val) {};
  int row, col, val;
};

bool operator < (const Element &lh, const Element &rh);
bool operator == (const Element &lh, const Element &rh);


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
////////////////////////// CMatrix_CSC ////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
// once created, can not change

class CMatrix_CSC {
public:
  // arr should be sorted, from  top-to-bottom, left-to-right
 CMatrix_CSC(): m(0), n(0), nnz(0), val(NULL), row(NULL), col_ptr(NULL) {};
  
  // move constructor
 CMatrix_CSC(CMatrix_CSC&& A): m(A.m), n(A.n), nnz(A.nnz) {
    val = A.val;
    A.val = NULL;
    row = A.row;
    A.row = NULL;
    col_ptr = A.col_ptr;
    A.col_ptr = NULL;
  };
  CMatrix_CSC& operator =(CMatrix_CSC&& A) {
    m = A.m;
    n = A.n;
    nnz = A.nnz;
    delete val;
    val = A.val;
    A.val = NULL;
    delete row;
    row = A.row;
    A.row = NULL;
    delete col_ptr;
    col_ptr = A.col_ptr;
    A.col_ptr = NULL;
    return *this;
  }
  // copy constructor
 CMatrix_CSC(const CMatrix_CSC &A): m(A.m), n(A.n), nnz(A.nnz) {
    val = new int[A.nnz];
    std::copy(A.val, A.val + A.nnz, val);
    row = new int[A.nnz];
    std::copy(A.row, A.row + A.nnz, row);
    col_ptr = new int[A.n + 1];
    std::copy(A.col_ptr, A.col_ptr + A.n + 1, col_ptr);
  };
  CMatrix_CSC(CVector_iter _val, CVector_iter _row, CVector_iter _col, int _nnz, int _m,  int _n);
  CMatrix_CSC(int *_val, int *_row, int *_col, int _nnz, int _m,  int _n);

  CMatrix_CSC operator[](int i) const {// return i_th column
    CVector _val, _row, _col;
    for (int t = col_ptr[i]; t < col_ptr[i + 1]; ++t) {
      _val.push_back(val[t]);
      _row.push_back(row[t]);
      _col.push_back(0);
    }
    return CMatrix_CSC(_val.begin(), _row.begin(), _col.begin(), _val.size(), m, 1);
  };

  ~CMatrix_CSC() {
    delete[] val;
    delete[] row;
    delete[] col_ptr;
  };
  int m, n;
  int nnz;
  int *val;
  int *row;
  int *col_ptr;
};


CMatrix_CSC operator *(const CMatrix_CSC &A, const CMatrix_CSC &B);

void thresh_mult(SparseVector &result, CMatrix_CSC &csc, int *v_s, int *v_e, double thresh);




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
  CMatrix_COO(std::vector<Element>::iterator arr_s, std::vector<Element>::iterator arr_e,
	      int m, int n): m(m), n(n) {
    data.assign(arr_s, arr_e);
  };
  void sortByColumnRow();
  void sortByRowColumn();
  int size() { return data.size(); }
  Element & back() {
    return data.back();
  }
  void push_back(Element e) { 
    data.push_back(e); 
    if (e.row > m)
      m = e.row;
    if (e.col > n)
      n = e.col;
  }
  Element& operator[] (int i) { return data[i]; }
  friend std::ostream& operator << (std::ostream & os, CMatrix_COO coo);
  void print() {
    for (const Element &e: data) 
      std::cout << e.row << " " << e.col << " " << e.val << std::endl;
  }
  void set_mn(int _m, int _n) {
    m = _m;
    n = _n;
  }
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


// r = P * vec
// result = r[r > thresh]
// note:  P has to be sorted by RowColumn
// complexity O(nnz(P))
void thresholdMult(SparseVector &result, CMatrix_COO &P, SparseVector &vec, double thresh);

std::ostream& operator << (std::ostream &os,  const SparseVector &vec); 

// sum rows of an coo matrix
void sumRows_Coo(int *result, CMatrix_COO &P);
// inner product v * sv
int inner_prod(int *v, SparseVector &sv);
CMatrix_COO toCoo(const CMatrix_CSC &mat);

#endif


















