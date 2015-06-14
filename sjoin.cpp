#include "sjoin.h"



void buildInvertedIndex(InvertedIndex& inIdx, const CMatrix_CSC& P, const CMatrix_CSC& q, int n) {
  if (n > P.n)
    throw "n is too large!";

  inIdx.clear();

  for (int i = 0; i < n; ++i) {
    int c = q.row[i]; // column of P should be processed
    for (int j = P.col_ptr[c]; j < P.col_ptr[c + 1]; ++j)
      inIdx[c].push_back(P.row[j]);
  }
} 




// n-prefix - n-prefix scheme
// theta - similarity threshold
CMatrix_COO prefix_sjoin(const CMatrix_CSC& P, const CMatrix_CSC& Q, int n_prefix, int theta) {
  CMatrix_COO res(P.m, Q.n); // to keep the result
  CMatrix_CSC&& PT = P.T();

  for (int i = 0; i < Q.n; ++i) {

    // setup n_token, first n_tokens will be used to build Inverted index
    int n_token = Q.get_nnz(i) - theta + n_prefix;

    if (n_token < n_prefix) // impossible to have similiar pairs
      continue;

    InvertedIndex inIdx;
    buildInvertedIndex(inIdx, P, Q[i], n_token);

    std::map<int, int> counter;
    for (auto it = inIdx.begin(); it != inIdx.end(); ++it) {
      std::vector<int>& vec = it->second;
      for (auto it_vec = vec.begin(); it_vec != vec.end(); ++it_vec)
	counter[*it_vec]++;
    }
    // not choose candidates and verify
    for (auto it = counter.begin(); it != counter.end(); ++it)
      if (it->second >= n_prefix) { // it's a legal candidates
	//TODO, now I assume P = Q^T, try to generalizd it
	int v = inner_product(PT[it->first], PT[i]);
	if (v >= theta) {
	  res.push_back(Element(it->first, i, v));
	}
      }
  }

  return res;
}














