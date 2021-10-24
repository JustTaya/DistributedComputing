import java.util.Random;

public class Vector {
    private final int[] vector;

    private final int size;

    Vector(int size) {
        this.size = size;
        this.vector = new int[size];
    }

    Vector(int[] vector) {
        this.size = (int) Math.sqrt(vector.length);
        this.vector = vector;
    }

    public void fillRandom(int maxNumber) {
        var rand = new Random(System.currentTimeMillis());
        for (int i = 0; i < size; i++)
            vector[i] = rand.nextInt(maxNumber);
    }

    public int[] getVector() {
        return vector;
    }

    public int getSize() {
        return size;
    }

    public int get(int i) {
        return vector[i];
    }

    public void set(int i, int value) {
        vector[i] = value;
    }

    public void add(int i, int value) {
        vector[i] += value;
    }
}
