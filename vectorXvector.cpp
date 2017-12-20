#include <mpi.h>
#include <cstdlib>
#include <iostream>
#include <ctime>
#include <stdio.h>
#include <iomanip>

void rndVector(double *vector, int vector_size) {
	for (unsigned int i = 0; i < vector_size; i++) {
		//vector[i] = ((double)rand() / (double)(RAND_MAX));
		vector[i] =1.0f;
	}
}

void rndMatrixRow(int size, double *matrix) {
	for (unsigned int i = 0; i < size; i++)	{
		for (unsigned int j = 0; j < size; j++)	{
			matrix[i*size+j] = 1.0f;
		}
	}
}

void printVector(double*vector, int vector_size) {
	for (unsigned int i = 0; i < vector_size; i++) {
		std::cout << "|" << std::setw(5) << vector[i] << "|" << std::endl;
	}
}

void initVectors(double*vector1, double*vector2,int ProcRank,int size) {
	if (ProcRank == 0) {
		rndVector(vector1, size);
		rndVector(vector2, size);
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

//distribute data and reduce on proc
double sMulVectors(double*vector1, double*vector2,int ProcRank,int ProcSize,int size,double localResult,double result) {
	int segment = size / ProcSize; //get yummy portions 
	int length1=segment * ProcRank;
	int length2=segment * (ProcRank+1);

	MPI_Bcast(vector1, size, MPI_DOUBLE, 0, MPI_COMM_WORLD); //send vectors 
	MPI_Bcast(vector2, size, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	if (ProcRank == ProcSize - 1) {
		length2 = size;
	}

	for (int i = length1; i < length2; i++) {
		localResult += vector1[i] * vector2[i]; 
		//std::cout << "proc result:" <<localResult << std::endl;
	}
		
	
	MPI_Reduce(&localResult, &result, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD); //collection local results in one total result
	return result;
}

void benchmark(int procrank,int procsize) {

	int size = 7000; 
	double tStart, tEnd, duration;
	double *vector1=0;
	double *vector2=0;
	vector1 = new double[size]; 
	vector2 = new double[size];
	double localResult = 0.0;
	double totalResult = 0.0;

	initVectors(vector1, vector2, procrank, size);

	tStart = MPI_Wtime();
	totalResult = sMulVectors(vector1, vector2, procrank, procsize, size, localResult, totalResult);
	tEnd = MPI_Wtime();

	duration = tEnd - tStart; 

	if (procrank == 0) {

		std::cout.precision(10);
		std::cout << "RESULT:" << totalResult << std::endl;
		std::cout << "\nDURATION:" << duration << std::endl;
	}

	delete[] vector1;
	delete[] vector2;

}

int main(int argc, char *argv[]) {
	MPI_Init(&argc, &argv);
	srand(time(NULL));
	int procsize;
	int procrank;

	MPI_Comm_size(MPI_COMM_WORLD, &procsize);
	MPI_Comm_rank(MPI_COMM_WORLD, &procrank);

	benchmark(procrank, procsize); //matrix matrix

	MPI_Finalize();
	return 0;
}
