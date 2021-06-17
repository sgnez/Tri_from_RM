# Tri from RM
SageMath codes for construction of triorthogonal codes (and computing their distances) from the Reed-Muller codewords.

Please make sure to install SageMath prior to running the code.
https://www.sagemath.org

If you are new to SageMath, or running mac Big Sur or a later macOS, it might be easier to install SageMath through Miniconda or Anaconda. 

The main file is the Jupyter notebook, and codebuilder.sage. 

The C++ code classes.cpp and related library polynomial.h are fast multithread C++ codes that find the affine conjugacy classes of polynomials. 
C++ code is designed to be called from th Jupyter notebooks. See the notebook.
Please do not modify codebuilder.py yourself. It should be automatically constructed from codebuilder.sage. See codebuilder.sage for instruction.
