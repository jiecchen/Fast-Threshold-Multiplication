#include "sjoin.h"



// build inverted index incrementally
InvertedIndex buildInvertedIndex(const CMatrix_CSC& M, int n_prefix, int theta) {
  if (n_prefix > M.m)
    throw "n_token is too large!";
  
  InvertedIndex result;

  for (int i = 0; i < M.n; ++i) {
    int n_token = M.get_nnz(i) - theta + n_prefix;

    if (n_token < n_prefix) // impossible to have similar pairs, skip this column
      continue;

    for (int j = M.col_ptr[i]; j < M.col_ptr[i + 1]; ++j)
      if (n_token-- > 0)
	result[M.row[j]].push_back(i);
      else 
	break;
  }

  return result;
} 



// each column of M is a data point
CMatrix_COO prefix_sjoin(const CMatrix_CSC& M, int n_prefix, int theta) {
  CMatrix_COO result(M.n, M.n); // to keep the result
  InvertedIndex&& inIdx = buildInvertedIndex(M, n_prefix, theta);

  for (int i = 0; i < M.n; ++i) {
    std::map<int, int> counter;
    int n_token = M.get_nnz(i) - theta + n_prefix;
    if (n_token < n_prefix) // impossible to have pairs, skip
      continue;

    CMatrix_CSC&& Mi = M[i];

    for (int j = M.col_ptr[i]; j < M.col_ptr[i + 1]; ++j) {
      for (auto it = inIdx[M.row[j]].begin(); it != inIdx[M.row[j]].end(); ++it)
	counter[*it]++;
    }
   
    // now choose candidates and verify
    for (auto it = counter.begin(); it != counter.end(); ++it)
      if (it->second >= n_prefix) { // it's a legal candidates
	VAL_TYPE v = inner_product(M[it->first], Mi);
	if (v >= theta) {
	  result.push_back(Element(it->first, i, v));
	}
      }
     
  }

  return result;
}
















