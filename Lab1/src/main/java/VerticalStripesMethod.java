import mpi.MPI;

public class VerticalStripesMethod {
    public static void compute(String[] args, int matrixSize, int maxNum) throws mpi.MPIException {
        MPI.Init(args);

        int rank = MPI.COMM_WORLD.Rank();
        int size = MPI.COMM_WORLD.Size();

        if (rank == 0) {

        }

        MPI.Finalize();
    }

    void processInitialization(Matrix inputMatrix, Vector inputVector, Vector outputVector, int size) {


        if (rank == 0) {
            inputMatrix = new Matrix();
        }
    }
}
