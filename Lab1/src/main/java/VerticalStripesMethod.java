import mpi.MPI;
import mpi.MPIException;

import java.util.concurrent.TimeUnit;

public class VerticalStripesMethod {

    public static void compute(String[] args, int matrixSize, Integer maxNum) throws MPIException {
        Matrix matrix = new Matrix();
        Vector vector = new Vector();
        Vector resultVector = new Vector();

        Vector processRows = new Vector();
        Vector processResult = new Vector();

        Integer rowNum = null;
    
        MPI.Init(args);

        int rank = MPI.COMM_WORLD.Rank();
        Integer size = MPI.COMM_WORLD.Size();

        processInitialization(matrix, vector, resultVector, processRows, processResult, matrixSize, rowNum, rank, size, maxNum);
        long start = System.nanoTime();
        parallelResultCalculation(processRows, vector, processResult, matrixSize, rowNum);
        resultReplication(processResult, resultVector, matrixSize, rowNum, size, rank);
        long finish = System.nanoTime();

        MPI.Finalize();

//        System.out.println(resultVector);
        System.out.println(TimeUnit.NANOSECONDS.toMicros(finish - start));
    }

    private static void processInitialization(Matrix inputMatrix, Vector inputVector, Vector resultVector,
                                              Vector processRows, Vector processResult, int matrixSize, Integer rowNum,
                                              int rank, Integer size, int maxValue) throws MPIException {
        int restRows;

        if (rank == 0) {
            do {
                if (matrixSize < Runtime.getRuntime().availableProcessors()) {
                    System.out.println("Size of the objects must be greater than number of processes!");
                }
                size = MPI.COMM_WORLD.Size();
            } while (matrixSize < size);
        }

        MPI.COMM_WORLD.Bcast(matrixSize, 0, 1, MPI.INT, 0);

        restRows = matrixSize;
        for (int i = 0; i < rank; i++) {
            restRows = restRows - restRows / (size - rank);
        }
        rowNum = restRows / (size - rank);

        inputMatrix.init(matrixSize);
        inputVector.init(matrixSize);
        resultVector.init(matrixSize);

        processResult.init(matrixSize);
        processRows.init(rowNum * matrixSize);

        if (rank == 0) {
            inputMatrix.fillRandom(maxValue);
            inputVector.fillRandom(maxValue);
        }
    }

    private static void dataDistribution(Matrix inputMatrix, Vector processRows, Vector inputVector, int matrixSize, Integer rowNum, int size, int rank) throws MPIException {
        int restRows = matrixSize;

        MPI.COMM_WORLD.Bcast(inputVector, 0, matrixSize, MPI.DOUBLE, 0);

        int[] sendInd = new int[size];
        int[] sendNum = new int[size];

        rowNum = matrixSize / size;

        sendNum[0] = rowNum * matrixSize;
        sendInd[0] = 0;

        for (int i = 1; i < size; i++) {
            restRows -= rowNum;
            rowNum = restRows / (size - i);
            sendNum[i] = rowNum * matrixSize;
            sendInd[i] = sendInd[i - 1] + sendNum[i - 1];
        }

        MPI.COMM_WORLD.Scatterv(inputMatrix, 0, sendNum, sendInd, MPI.DOUBLE, processRows,
                sendNum[rank], 0, MPI.DOUBLE, 0);
    }

    private static void parallelResultCalculation(Vector processRows, Vector inputVector, Vector processResult, int matrixSize, int rowNum) {
        for (int i = 0; i < rowNum; i++) {
            processResult.set(i, 0);
            for (int j = 0; j < matrixSize; j++) {
                processResult.add(i, processRows.get(i * matrixSize + j) * inputVector.get(j));
            }
        }
    }

    private static void resultReplication(Vector processResult, Vector resultVector, int matrixSize, int rowNum, int size, int rank) throws MPIException {
        int restRows = size;

        int[] receiveInd = new int[size];
        int[] receiveNum = new int[size];

        receiveInd[0] = 0;
        receiveNum[0] = matrixSize / size;

        for (int i = 1; i < size; i++) {
            restRows -= receiveNum[i - 1];
            receiveNum[i] = restRows / (size - i);
            receiveInd[i] = receiveInd[i - 1] + receiveNum[i - 1];
        }

        MPI.COMM_WORLD.Allgatherv(processResult, receiveNum[rank], receiveNum[rank], MPI.DOUBLE, resultVector, receiveNum[rank], receiveNum, receiveInd, MPI.DOUBLE);
    }
}
