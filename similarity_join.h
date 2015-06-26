#ifndef __SIMILARITY_JOIN__
#define __SIMILARITY_JOIN__
#include "Algebra.h"


CMatrix_COO prefix_sketch_join(const CMatrix_CSC& _M, int n_prefix, int theta, int w=10);

CMatrix_COO prefix_matrix_join(const CMatrix_CSC& _M, int n_prefix, int theta);

#endif



















