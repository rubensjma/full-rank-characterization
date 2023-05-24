# full-rank-characterization
This code is concerned with the characterization of the rank of convex combinations of given full-rank matrices. 

The test is based on the solution of Quadratic Programming problems. After a finite number of iterations,
two outcomes are possible: (i) no convex combination is rank-deficient within the unit simplex or (ii)
convex coefficients are returned, such that the smallest singular of the convex combination of the matrices is
smaller than an arbitrary user-defined tolerance.

The file examples_manuscript.m is the main routine to run three examples. Requires MATLAB (https://www.mathworks.com/products/matlab.html) including the Optimization Toolbox, and MPT (http://control.ee.ethz.ch/~mpt).

An alternative method for comparison purposes for the second example is presented in the file example2.m, which requires MATLAB (https://www.mathworks.com/products/matlab.html), YALMIP (https://yalmip.github.io/) and the SeDuMi solver (https://github.com/sqlp/sedumi). The reference for the alternative method is "Comments on 'Less conservative conditions for robust LQR-state-derivative controller design: an LMI approach’ and new sufficient LMI conditions for invertibility of a convex combination of matrices." (https://doi.org/10.1080/00207721.2021.2023689), and the example also requires the function remark3.m available in Appendix 2 of the paper.
