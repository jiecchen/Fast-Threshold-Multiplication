Notes
========

## ideas
+ The problem to use CountSketch  is that it is difficult to
  decide when we should use exact algorithm

+ Need to balance running time for "Sketching Matrix" && "Recovery"
  not necessary to merge just two neighbors, instead, we can merge many

+ basically, after we get R1 = [r0, r1, r2, r3, ...], we automatically set the
  value of w, so that running time T = exact + sketch can be minimized.
  group the columns in Q be Qh and Ql, the T will be
  T = ||R*Qh||_1 + beta * mu * w * nnz(Ql)
  beta is a constant to be estimated.
  also need to optimize STEP_SIZE

+ add comparasison to matlab if my implementation of sparse multiplication
  is truely faster than matlab.

+ read some papers for paralell sparse matrix multiplication,
  discuss how to extend our algorithm in that seting.
  read this http://gauss.cs.ucsb.edu/~aydin/spgemm_sisc12.pdf
  for parallel sparse matrix multiplication.
  need new formula to estimate the running time.
  also this http://www.cs.berkeley.edu/~yelick/cs267-sp04/lectures/

+ why not splitting the vector R_i to blocks? Will have some speedup?

+ compare with the nips best paper?


## What I learned from this work
+ implement the algorithm will help you to **think**.
  easy to see the bottleneck.

+ work on the problem **after** we found strong motivation.

+ need sort of theory to support the idea.


