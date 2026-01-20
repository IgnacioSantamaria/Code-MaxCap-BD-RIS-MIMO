# Code-MaxCap-BD-RIS-MIMO/paper_ICASSP26
This folder contains code related to the paper "Riemannian optimization on the manifold of unitary and symmetric matrices with application to BD-RIS-assisted systems,"  I. Santamaria, M. Soleymani, E. Jorswieck, J. Gutierrez, and C. Beltran.

# Content of Code Package
The code is implemented in Matlab.
The function implementing Algorithm 1 in the paper is OptimizeBDRIS_MOUs_FullRank.m. The function OptimizeBDRIS_MOU implements a MO algorithm on the manifold of unitary (but not symmetric) matrices. Therefore, we need a final projection step to obtain a symmetric and unitary matrix. For comparison, we also include the code of the SPAWC'24 paper, but it is rather slow for large M so we recommend using the ICASSP'26 version instead. The man script is Script_ICASSP26. 



