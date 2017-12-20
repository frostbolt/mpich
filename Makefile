default: matrixXmatrix.cpp matrixXvector.cpp vectorXvector.cpp
		mpicxx matrixXmatrix.cpp -o matrixXmatrix && mpicxx matrixXvector.cpp -o matrixXvector && mpicxx vectorXvector.cpp -o vectorXvector

mxm:
		mpiexec -n 4 ./matrixXmatrix

mxv:
		mpiexec -n 4 ./matrixXvector

vxv:
		mpiexec -n 4 ./vectorXvector