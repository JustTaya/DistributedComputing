import mpi.MPIException;

public class Main {
    public static void main(String[] args) {
        try {
            VerticalStripesMethod.compute(args, 20, 20);
        } catch (MPIException e) {
            e.printStackTrace();
        }
    }
}
