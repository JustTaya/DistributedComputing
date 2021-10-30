import java.util.Random;

public class Vector {
    private double[] vector;
    
    private int size;
    
    Vector() {
    
    }
    
    Vector(int size) {
        this.size = size;
        this.vector = new double[size];
    }
    
    Vector(double[] vector) {
        this.size = vector.length;
        this.vector = vector;
    }
    
    public void init(int size) {
        this.size = size;
        this.vector = new double[size];
    }
    
    public void fillRandom(int maxNumber) {
        Random rand = new Random(System.currentTimeMillis());
        for (int i = 0; i < size; i++)
            vector[i] = rand.nextDouble() * maxNumber;
    }
    
    public double[] getVector() {
        return vector;
    }
    
    public int getSize() {
        return size;
    }
    
    public double get(int i) {
        return vector[i];
    }
    
    public void set(int i, int value) {
        vector[i] = value;
    }
    
    public void add(int i, double value) {
        vector[i] += value;
    }
}
