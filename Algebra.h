#ifndef __ALGEBRA_H__
#define __ALGEBRA_H__
#include <vector>
#include <cstdlib>
#include <iostream>


// class Matrix_CSC {
// public:
//   // arr has to be sorted
//   Matrix_CSC(Element *arr, int _nnz, int _m, int _n);
//   Matrix_CSC operator * (const Matrix_CSC &m);
//   friend std::ostream& operator << (std::ostream &os, const Matrix_CSC &m);  
//   ~Matrix_CSC() {
//     delete[] values;
//     delete[] row_ind;
//     delete[] col_ptr;
//   }

// private:
//   int m, n; // m by n matrix
//   int nnz; // # of non-zero entries
//   int* values;
//   int* row_ind;
//   int* col_ptr;
// };

class Element {
public:
  Element(): row(0), col(0), val(0) {};
  Element(int row, int col, int val): row(row), col(col), val(val) {};
  int row, col, val;
};

typedef std::vector<int> CVector; // dense vector
typedef CVector::iterator CVector_iter;
typedef std::vector<Element> CMatrix_COO;

// sparse matrix using CSR format
class CMatrix_CSR {
public:
  CMatrix_CSR(Element *arr, int _nnz, int _m,  int _n);

  // sparse matrix *  dense vector
  CVector operator *(const CVector &vec);
  CMatrix_COO toCOO() const;
  friend std::ostream& operator << (std::ostream &os, const CMatrix_CSR &mat);  

  ~CMatrix_CSR() {
    delete[] values;
    delete[] col_ind;
    delete[] row_ptr;
  }

private:
  int m, n;
  int nnz;
  int *values;
  int *col_ind;
  int *row_ptr;
};




// override << for CVector
std::ostream& operator << (std::ostream & os, CVector &vec);




#endif

















