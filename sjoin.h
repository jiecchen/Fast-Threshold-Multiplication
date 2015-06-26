// to implement similarity join algorithms
#include "Algebra.h"


#include <vector>
#include <map>

typedef std::map<int, std::vector<int>> InvertedIndex;


// n-prefix - n-prefix scheme
// theta - similarity threshold
// each column of M is a data point
CMatrix_COO prefix_sjoin(const CMatrix_CSC& M, int n_prefix, int theta);



