#include <mpi.h>
#include <cstdlib>
#include <iostream>
#include <ctime>
#include <stdio.h>
#include <iomanip>

void matrixMultiply(double *first, double *second, double *result, int size, int startingPoint) {
	for (unsigned int i = 0; i<startingPoint; i++) {
		for (unsigned int j = 0; j<size; j++) {
			for (unsigned int k = 0; k<size; k++) {
				result[i*size + j] += first[i*size + k] * second[k*size + j];
			}
		}
	}
}

void rndMatrixCol(int size, double *matrix) {
	for (unsigned int i = 0; i < size; i++)	{
		for (unsigned int j = 0; j < size; j++)	{
			matrix[j*size + i] = 1.0f;
		}
	}
}

void rndMatrixRow(int size, double *matrix) {
	for (unsigned int i = 0; i < size; i++)	{
		for (unsigned int j = 0; j < size; j++)	{
			matrix[i*size+j] = 1.0f;
		}
	}
}


void resultMatrixInit(int size, double *first, double *second, int ProcRank) {
	if (!ProcRank) {
		rndMatrixCol(size, first);
		rndMatrixRow(size, second);
	}
}

void matrixMul(int size,double *first, double *second, int ProcRank, int ProcSize, double *result) {	
	int currentSegment, startingPoint, endingPoint;
	int sendingSize, add;	
	MPI_Bcast(&size, 1, MPI_INT, 0, MPI_COMM_WORLD); //отправка размера
	currentSegment = (size / ProcSize) + (size%ProcSize > ProcRank);
	startingPoint = size*currentSegment; //верхний предел
	endingPoint = size*size; //нижний предел
	second = new double[endingPoint];

	if (ProcRank) {
		first = new double[startingPoint];
		result = new double[startingPoint];
	} else {
		first = new double[endingPoint];
		result = new double[endingPoint];
		for (unsigned int i = 0; i<endingPoint; i++) {
			first[i] = 1.0;
			second[i] = 1.0;
		}
	}

	MPI_Bcast(second, size * size, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	if (!ProcRank) {
		sendingSize =startingPoint;
		for (unsigned int i = 1; i<ProcSize; i++) {
			add = size * ((size / ProcSize) + (size % ProcSize > ProcRank));
			MPI_Send(first + sendingSize, add, MPI_DOUBLE, i, 31337,MPI_COMM_WORLD);
			sendingSize += add;
		}
	} else { 
		MPI_Recv(first, startingPoint, MPI_DOUBLE, 0, 31337, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
	}

	for (unsigned int i = 0; i<startingPoint; i++) {
		result[i] = 0.0;
	}

	matrixMultiply(first, second, result, size, currentSegment);

	if (!ProcRank) {
		sendingSize = startingPoint;

		for (unsigned int i = 1; i<ProcSize; i++) {
			add = size * ((size / ProcSize) + (size % ProcSize > ProcRank));
			MPI_Recv(result + sendingSize, add, MPI_DOUBLE, i, 42,MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			sendingSize += add;
		}
	} else { 
		MPI_Send(result,startingPoint, MPI_DOUBLE, 0, 42, MPI_COMM_WORLD);
	}

	if (!ProcRank) {
		for (unsigned int i = 0; i<size*size; i++) {
			printf("%5.1lf\t", result[i]);
			if (i % size == 0) printf("\n");
		}
	}	
}

void benchmark(int procrank, int procsize) {
	int size = 500;
	double tStart, tEnd, duration;
	double *matrixA=0;
	double *matrixB=0;
	double *matrixC=0;

	//resultMatrixInit(size, matrixA, matrixB, procrank);

	tStart = MPI_Wtime();
	matrixMul(size, matrixA, matrixB, procrank, procsize, matrixC);
	tEnd = MPI_Wtime();
	duration = tEnd - tStart;
	//print_matrix(size, matrixC);
	std::cout << "\ntime:" << duration << std::endl;
	delete[] matrixA;
	delete[] matrixB;
	delete[] matrixC;	

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
