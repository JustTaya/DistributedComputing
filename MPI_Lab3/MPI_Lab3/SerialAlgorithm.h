#pragma once
class SerialAlgorithm
{
public:
	double execute(int argc, char* argv[]);
private:
	int* pSerialPivotPos;
	int* pSerialPivotIter;

	void ProcessInitialization(double*& pMatrix, double*& pVector,
		double*& pResult, int& Size);

	void RandomDataInitialization(double* pMatrix, double* pVector, int Size);

	void SequentialResultCalculation(double* pMatrix, double* pVector,
		double* pResult, int Size);

	void ProcessTermination(double* pMatrix, double* pVector, double* pResult);

	// Gaussian Elimination
	void SerialGaussianElimination(double* pMatrix, double* pVector, int Size);

	int FindPivotRow(double* pMatrix, int Size, int Iter);

	void SerialColumnElimination(double* pMatrix, double* pVector, int Pivot,
		int Iter, int Size);

	// Back Substitution
	void SerialBackSubstitution(double* pMatrix, double* pVector,
		double* pResult, int Size);

};

