#pragma once
class ParallelGaussAlgorithm
{
public:
	double execute(int argc, char* argv[]);
private:
	int ProcNum;
	int ProcRank;
	int* pParallelPivotPos;
	int* pProcPivotIter;

	int* pProcInd;
	int* pProcNum;

	void ProcessInitialization(double*& pMatrix, double*& pVector,
		double*& pResult, double*& pProcRows, double*& pProcVector,
		double*& pProcResult, int& Size, int& RowNum);

	void DataDistribution(double* pMatrix, double* pProcRows, double* pVector,
		double* pProcVector, int Size, int RowNum);

	void ResultCollection(double* pProcResult, double* pResult);

	void RandomDataInitialization(double* pMatrix, double* pVector, int Size);

	void ParallelResultCalculation(double* pProcRows, double* pProcVector,
		double* pProcResult, int Size, int RowNum);

	void ProcessTermination(double* pMatrix, double* pVector, double* pResult,
		double* pProcRows, double* pProcVector, double* pProcResult);
	
	// Gausian elimination
	void ParallelGaussianElimination(double* pProcRows, double* pProcVector, int Size, int RowNum);

	void ParallelEliminateColumns(double* pProcRows, double* pProcVector, double* pPivotRow, int Size, int RowNum, int Iter);

	void ParallelBackSubstitution(double* pProcRows, double* pProcVector,
		double* pProcResult, int Size, int RowNum);

	void FindBackPivotRow(int RowIndex, int Size, int& IterProcRank, int& IterPivotPos);
};


