#include <stdio.h>
#include <stdlib.h>
#include <conio.h>
#include <time.h>
#include <math.h>
#include "SerialAlgorithm.h"

double SerialAlgorithm::execute(int argc, char* argv[])
{
	double* pMatrix;
	double* pVector;
	double* pResult;
	int Size;
	time_t start, finish;
	double duration;

	int N[] = { 10,100, 500, 1000, 1500, 2000,2500,3000 };
	for (int i = 0; i < 8; i++) {
		Size = N[i];
		printf("Serial Gauss algorithm for solving linear systems\n");
		ProcessInitialization(pMatrix, pVector, pResult, Size);

		start = clock();
		SequentialResultCalculation(pMatrix, pVector, pResult, Size);
		finish = clock();
		duration = (finish - start) / double(CLOCKS_PER_SEC);

		printf("\n Time of execution: %f\n", duration);

		ProcessTermination(pMatrix, pVector, pResult);
	}
}

void SerialAlgorithm::ProcessInitialization(double*& pMatrix, double*& pVector, double*& pResult, int& Size)
{
	pMatrix = new double[Size * Size];
	pVector = new double[Size];
	pResult = new double[Size];

	RandomDataInitialization(pMatrix, pVector, Size);
}

void SerialAlgorithm::RandomDataInitialization(double* pMatrix, double* pVector, int Size)
{
	int i, j;
	srand(unsigned(clock()));
	for (i = 0; i < Size; i++) {
		pVector[i] = rand() / double(1000);
		for (j = 0; j < Size; j++) {
			if (j <= i)
				pMatrix[i * Size + j] = rand() / double(1000);
			else
				pMatrix[i * Size + j] = 0;
		}
	}
}

void SerialAlgorithm::SequentialResultCalculation(double* pMatrix, double* pVector, double* pResult, int Size)
{
	pSerialPivotPos = new int[Size];
	pSerialPivotIter = new int[Size];
	for (int i = 0; i < Size; i++) {
		pSerialPivotIter[i] = -1;
	}
	SerialGaussianElimination(pMatrix, pVector, Size);

	SerialBackSubstitution(pMatrix, pVector, pResult, Size);

	delete[] pSerialPivotPos;
	delete[] pSerialPivotIter;
}

void SerialAlgorithm::ProcessTermination(double* pMatrix, double* pVector, double* pResult)
{
	delete[] pMatrix;
	delete[] pVector;
	delete[] pResult;
}

void SerialAlgorithm::SerialGaussianElimination(double* pMatrix, double* pVector, int Size)
{
	int Iter;
	int PivotRow;    
	for (Iter = 0; Iter < Size; Iter++) {
		PivotRow = FindPivotRow(pMatrix, Size, Iter);
		pSerialPivotPos[Iter] = PivotRow;
		pSerialPivotIter[PivotRow] = Iter;
		SerialColumnElimination(pMatrix, pVector, PivotRow, Iter, Size);
	}
}

int SerialAlgorithm::FindPivotRow(double* pMatrix, int Size, int Iter)
{
	int PivotRow = -1;
	int MaxValue = 0;
	int i;

	for (i = 0; i < Size; i++) {
		if ((pSerialPivotIter[i] == -1) &&
			(fabs(pMatrix[i * Size + Iter]) > MaxValue)) {
			PivotRow = i;
			MaxValue = fabs(pMatrix[i * Size + Iter]);
		}
	}
	return PivotRow;
}

void SerialAlgorithm::SerialColumnElimination(double* pMatrix, double* pVector, int Pivot, int Iter, int Size)
{
	double PivotValue, PivotFactor;
	PivotValue = pMatrix[Pivot * Size + Iter];
	for (int i = 0; i < Size; i++) {
		if (pSerialPivotIter[i] == -1) {
			PivotFactor = pMatrix[i * Size + Iter] / PivotValue;
			for (int j = Iter; j < Size; j++) {
				pMatrix[i * Size + j] -= PivotFactor * pMatrix[Pivot * Size + j];
			}
			pVector[i] -= PivotFactor * pVector[Pivot];
		}
	}
}

void SerialAlgorithm::SerialBackSubstitution(double* pMatrix, double* pVector, double* pResult, int Size)
{
	int RowIndex, Row;
	for (int i = Size - 1; i >= 0; i--) {
		RowIndex = pSerialPivotPos[i];
		pResult[i] = pVector[RowIndex] / pMatrix[Size * RowIndex + i];
		for (int j = 0; j < i; j++) {
			Row = pSerialPivotPos[j];
			pVector[j] -= pMatrix[Row * Size + i] * pResult[i];
			pMatrix[Row * Size + i] = 0;
		}
	}
}
