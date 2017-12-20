#include <mpi.h>
#include <cstdlib>
#include <iostream>
#include <ctime>
#include <stdio.h>
#include <iomanip>


void rndVector(double *vector, int vector_size) {
	for (unsigned int i = 0; i < vector_size; i++) {
		//vector[i] = ((double)rand() / (double)(RAND_MAX));
		vector[i] = 1.0f;
	}
}

void printVector(double *vector, int vector_size) {
	for (unsigned int i = 0; i < vector_size; i++) {
		std::cout << "|" << std::setw(5) << vector[i] << "|" << std::endl;
	}
}

void rndMatrixRow(int size, double *matrix) {
	for (unsigned int i = 0; i < size; i++)	{
		for (unsigned int j = 0; j < size; j++)	{
			matrix[i*size+j] = 1.0f;
		}
	}
}


void printMatrix(int size, double *matrix){ 
	for (unsigned int i = 0; i < size; i++)	{
		for (unsigned int j = 0; j < size; j++)	{
			std::cout << "|" << std::setw(3) << matrix[i*size+j] << "|";
		}
		std::cout << "\n";
	}
}

void initMatrixVector(int size, double *matrix, double *vector, int ProcRank){
	rndVector(vector, size);
	rndMatrixRow(size, matrix);		
}

void mulMatrixVector(int size,double *matrix,double *vector,int ProcRank,int ProcSize,double *resultVector){
	//колво строк на каждый проц
	int rowsPerProc = size / ProcSize;
	//отрезки для расчета
	int startingPoint = ProcRank * rowsPerProc;
	int endingPoint = (ProcRank == ProcSize - 1) ? (size- 1) : (startingPoint + rowsPerProc - 1);
	
	for (int i = startingPoint; i <= endingPoint; i++) {
		resultVector[i] = 0.0;
		for (int j = 0; j < size; j++)
			resultVector[i] += matrix[i * size + j] * vector[j];
	}

	int *rcounts = new int[ProcSize]; // массив, содержащий порцию от каждого процесса 
	int *displs = new int[ProcSize]; // массив смещений

	for (int i = 0; i < ProcSize; i++) {
		rcounts[i] = (i == ProcSize - 1) ? size - i * rowsPerProc : rowsPerProc;
		displs[i] = (i > 0) ? displs[i - 1] + rcounts[i - 1] : 0;
	}
	//сборка всеми процессорами группы из кусочков в резалт вектор
	MPI_Allgatherv(MPI_IN_PLACE, 0, MPI_DOUBLE, resultVector, rcounts, displs, MPI_DOUBLE, MPI_COMM_WORLD);
}

void benchmark(int procrank, int procsize) {
	int size = 1000;
	double tStart, tEnd, duration;
	double *vector;
	double *matrix=new double[size*size] ;
	double *resultVector;
	vector = new double[size];
	resultVector = new double[size];

	initMatrixVector(size, matrix, vector,procrank);
	tStart = MPI_Wtime();
	mulMatrixVector(size, matrix, vector, procrank, procsize, resultVector);

	tEnd = MPI_Wtime();

	duration = tEnd - tStart;

	printVector(resultVector, size);
	std::cout << "\nDURATION:" << duration << std::endl;

	delete[] vector;
	delete[] matrix;
	delete[] resultVector;
}

int main(int argc, char *argv[]) {
	int procsize;
	int procrank;

	MPI_Init(&argc, &argv);
	srand(time(NULL));

	MPI_Comm_size(MPI_COMM_WORLD, &procsize);
	MPI_Comm_rank(MPI_COMM_WORLD, &procrank);
 
	benchmark(procrank, procsize);

	MPI_Finalize();
	return 0;
}