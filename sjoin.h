// to implement similarity join algorithms
#include "Algebra.h"


#include <vector>
#include <map>

typedef std::map<int, std::vector<int>> InvertedIndex;


// build inverted index on P, use the first n tokens of vector q
void buildInvertedIndex(InvertedIndex& inIdx, const CMatrix_CSC& P, const CMatrix_CSC& q, int n);


// n-prefix - n-prefix scheme
// theta - similarity threshold
CMatrix_COO prefix_sjoin(const CMatrix_CSC& P, const CMatrix_CSC& Q, int n_prefix, int theta);






