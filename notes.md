Notes
========

## idea
+ The problem to use CountSketch  is that it is difficult to
  decide when we should use exact algorithm

+ Need to balance running time for "Sketching Matrix" && "Recovery"
  not necessary to merge just two neighbors, instead, we can merge many

+ basically, after we get R1 = [r0, r1, r2, r3, ...], we automatically set the
  value of w, so that running time T = exact + sketch can be minimized.
  group the columns in Q be Qh and Ql, the T will be
  T = ||R*Qh||_1 + beta * mu * w * nnz(Ql)
  beta is a constant to be estimated.



