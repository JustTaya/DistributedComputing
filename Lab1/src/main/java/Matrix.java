import java.util.Random;

public class Matrix {
    private double[][] matrix;
    
    private int size;
    
    Matrix() {
    
    }
    
    Matrix(int size) {
        this.size = size;
        matrix = new double[size][];
        for (int i = 0; i < size; i++)
            matrix[i] = new double[size];
    }
    
    Matrix(int[] lineMatrix) {
        size = (int) Math.sqrt(lineMatrix.length);
        
        matrix = new double[size][];
        for (int i = 0; i < size; i++) {
            matrix[i] = new double[size];
            System.arraycopy(lineMatrix, i * size, matrix[i], 0, size);
        }
    }
    
    public void init(int size) {
        this.size = size;
        matrix = new double[size][];
        for (int i = 0; i < size; i++)
            matrix[i] = new double[size];
    }
    
    public void copyFromLine(int[] lineMatrix) {
        for (int i = 0; i < size; i++) {
            System.arraycopy(lineMatrix, i * size, matrix[i], 0, size);
        }
    }
    
    public void fillRandom(int maxNumber) {
        Random rand = new Random(System.currentTimeMillis());
        for (int i = 0; i < size; i++)
            for (int j = 0; j < size; j++)
                matrix[i][j] = rand.nextDouble() * maxNumber;
    }
    
    public int[] getMatrix() {
        int[] linedMatrix = new int[size * size];
        for (int i = 0; i < size; i++)
            System.arraycopy(matrix[i], 0, linedMatrix, i * size, size);
        return linedMatrix;
    }
    
    public int getSize() {
        return size;
    }
    
    public double get(int i, int j) {
        return matrix[i][j];
    }
    
    public void set(int i, int j, int value) {
        matrix[i][j] = value;
    }
    
    public void add(int i, int j, int value) {
        matrix[i][j] += value;
    }
}
