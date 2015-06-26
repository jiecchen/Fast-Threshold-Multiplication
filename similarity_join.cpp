#include "similarity_join.h"
#include "Sketch.h"
#include "Algebra.h"
#include <vector>
#include "utils.h"


CMatrix_CSC remove_tails(const CMatrix_CSC& M, int n_prefix, int theta) {
  CMatrix_COO res(M.m, M.n);
  for (int i = 0; i < M.n; ++i) {
    int n_token = M.get_nnz(i) - theta  + n_prefix;

    if (n_token < n_prefix) // impossible to have similar pairs
      continue;

    // will only keep the first n_token non-zero entries in column i
    for (int j = M.col_ptr[i]; j < M.col_ptr[i + 1]; ++j) 
      if (n_token-- > 0) {
	res.push_back(Element(M.row[j], i, M.val[j]));	
      }
      else
	break;
  }
  //  std::cerr << "after remove tails: size = " << res.size() << std::endl;
  return CMatrix_CSC(res, true);
}


CMatrix_COO prefix_matrix_join(const CMatrix_CSC& _M, int n_prefix, int theta) {
  CTimer timer;
  CMatrix_COO result(_M.n, _M.n); // to keep the result

  // remove tails to speed up calculation
  timer.start("removing tails ...");
  CMatrix_CSC&& M = remove_tails(_M, n_prefix, theta);
  timer.stop();
  
  timer.start("verifying ...");

  int debug_ct = 0;

  CMatrix_CSC&& MT = M.T();
  CMatrix_CSC&& res = MT * M;
  for (int i = 0; i < res.n; ++i) {
    CMatrix_CSC&& _Mi = _M[i];
    for (int j = res.col_ptr[i]; j < res.col_ptr[i + 1]; ++j)
      if (res.val[j] >= n_prefix) {
	debug_ct++;
	VAL_TYPE v = inner_product(_M[res.row[j]], _Mi);
	if (v >= theta)
	  result.push_back(Element(res.row[j], i, v));
      }
  }
  timer.stop();

  std::cerr << "size of candidates: " << debug_ct << std::endl;

  return result;
}






// each column of _M is a data point
CMatrix_COO prefix_sketch_join(const CMatrix_CSC& _M, int n_prefix, int theta, int w) {
  CTimer timer;
  CMatrix_COO result(_M.n, _M.n); // to keep the result

  // remove tails to speed up calculation
  timer.start("removing tails ...");
  CMatrix_CSC&& M = remove_tails(_M, n_prefix, theta);
  timer.stop();
  CMatrix_CSC&& MT = M.T();
  

  /////////////////////////////////////////////////////
  /////////////  Using sketch /////////////////////////
  /////////////////////////////////////////////////////

  // create Ms[..] by merging neighbors
  std::vector<CMatrix_CSC*> Ms; // keep merged matrices
  Ms.push_back(new CMatrix_CSC(MT));


  
  while (Ms.back()->m > 1) {
    Ms.push_back(new CMatrix_CSC(mergeNeighbor(*Ms.back())));
  }


  
  std::vector<CMatrix_CSC> sk; // to sketched matrix
  std::vector<CMatrix_CSC*> CMs; // to keep count-min matrices




  // - create count-min
  // - sketch matrices
  timer.start("create sketch ...");
  for (unsigned int i = 0; i < Ms.size(); ++i) {
    // create count min sketch
    CMs.push_back(new CMatrix_CSC(createCountMin(w, _MU, Ms[i]->m)));
    sk.push_back( (*CMs[i] * *Ms[i]) * M );
  }
  timer.stop();



  timer.start("Do recovery & verification ...");
  // now do recovery and verification
  for (auto i = 0; i < M.n; ++i) {

    // use count-min to do recovery
    CMatrix_COO&& res = recover_column(CMs, sk, i, n_prefix);
    CMatrix_CSC&& _Mi = _M[i];
    CMatrix_CSC&& Mi = M[i];

    // do verification
    for (int j = 0; j < res.size(); ++j) 
      if (inner_product(M[res[j].row], Mi) >= n_prefix) {
	VAL_TYPE v = inner_product(_M[res[j].row], _Mi);
	if (v >= theta) {
	  result.push_back(Element(res[j].row, i, v));
	}
      }
    
  } 
  timer.stop();

  // return result as a CSC matrix
  return result;
}
