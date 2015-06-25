#ifndef __ALGEBRA_H__
#define __ALGEBRA_H__
#include <algorithm>
#include <iterator>
#include <vector>
#include <iostream>

#define VAL_TYPE double


////////////////////////////////////////////////////////////////////////////////
class Element {
public:
  Element(): row(0), col(0), val(0) {};
  Element(int row, int col, VAL_TYPE val): row(row), col(col), val(val) {};
  int row, col;
  VAL_TYPE val;
};

bool operator < (const Element &lh, const Element &rh);
bool operator == (const Element &lh, const Element &rh);


////////////////////////////////////////////////////////////////////////////////
// can be used to represent sparse vector
class VectorElement {
public:
  VectorElement(): ind(0), val(0) {};
  VectorElement(int ind, VAL_TYPE val): ind(ind), val(val) {};
  int ind; 
  VAL_TYPE val;
};

typedef std::vector<VAL_TYPE> CVector; // dense vector
typedef CVector::iterator CVector_iter;
typedef std::vector<int>::iterator Index_iter;
typedef std::vector<VectorElement> SparseVector;
typedef std::vector<VectorElement>::iterator SparseVector_iter;






///////////////////////////////////////////////////////////////////////////////
////////////////////////// CMatrix_COO ////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

class CMatrix_COO {
public:
  // constructors
  CMatrix_COO(): m(0), n(0) {};
  CMatrix_COO(int m, int n): m(m), n(n) {};
  CMatrix_COO(Element *arr_s, Element *arr_e, int m, int n): m(m), n(n) {
    data.assign(arr_s, arr_e);
  };

  CMatrix_COO(std::vector<Element>::iterator arr_s, std::vector<Element>::iterator arr_e,
	      int m, int n): m(m), n(n) {
    data.assign(arr_s, arr_e);
  };


  // sorting
  void sortByColumnRow();
  void sortByRowColumn();


  int size() { return data.size(); }
  Element & back() {
    return data.back();
  }

  void push_back(Element e);

  Element& operator[] (int i) { return data[i]; }
 
  friend std::ostream& operator << (std::ostream & os, CMatrix_COO coo);

  void print(double theta = 0) const;   // print entries > theta

  void set_mn(int _m, int _n) {
    m = _m;
    n = _n;
  }

  int get_m() { return m; }
  int get_n() { return n; }

  CMatrix_COO T();   // return transpose

private:
  int m, n;
  std::vector<Element> data;
};








///////////////////////////////////////////////////////////////////////////////
////////////////////////// CMatrix_CSC ////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
// once created, can not change

class CMatrix_CSC {
public:
  // arr should be sorted, from  top-to-bottom, left-to-right
 CMatrix_CSC(): m(0), n(0), nnz(0), val(NULL), row(NULL), col_ptr(NULL) {};
  
  // move constructor
  CMatrix_CSC(CMatrix_CSC&& A);

  // copy constructor
  CMatrix_CSC(const CMatrix_CSC &A);

  CMatrix_CSC(CVector_iter _val, Index_iter _row, Index_iter _col, int _nnz, int _m,  int _n, bool sorted=false);

  CMatrix_CSC(VAL_TYPE *_val, int *_row, int *_col, int _nnz, int _m,  int _n, bool sorted=false);


  CMatrix_CSC(CMatrix_COO& coo, bool sorted=false) {
    this->init(coo, sorted);

  };

  void init(CMatrix_COO& coo, bool sorted=false);

  // move assignment
  CMatrix_CSC& operator =(CMatrix_CSC&& A);

  // slicing
  CMatrix_CSC operator[](int i) const;
  CMatrix_CSC operator[](const std::vector<int>& idx) const;

  // Transpose
  CMatrix_CSC T() const;

  

  ~CMatrix_CSC() {
    delete[] val;
    delete[] row;
    delete[] col_ptr;
  };



  // return nnz of a column
  int get_nnz(int ind_col) const {
    return col_ptr[ind_col + 1] - col_ptr[ind_col];
  }


  
  // bad design, but much earsier to use
  int m, n;
  int nnz;
  VAL_TYPE *val;
  int *row;
  int *col_ptr;

private:
  template<typename T1, typename T2>
  void init(T1 _val, T2 _row, T2 _col, int _nnz, int _m,  int _n);

};


CMatrix_CSC operator *(const CMatrix_CSC &A, const CMatrix_CSC &B);



// convert CSC to COO format
CMatrix_COO toCoo(const CMatrix_CSC &mat);


// given to vector, return their inner product
VAL_TYPE inner_product(const CMatrix_CSC& a, const CMatrix_CSC& b, double speedup_thresh=1e+20);



std::ostream& operator << (std::ostream &os,  const VectorElement& e);

#endif


















