import java.util.Arrays;
import java.util.BitSet;

public class Population {

    Genotype[] encodedSolutions;

    public Population(int popSize){
        encodedSolutions = new Genotype[popSize];
    }

    public Genotype getSolution(int i) {
        return encodedSolutions[i];
    }

    public int Size() {
        return encodedSolutions.length;
    }

    public void Insert(int index, Genotype newSolution){
        encodedSolutions[index] = newSolution;
    }

    public boolean hasConverged(){
        //TODO
        // Check for convergence
        return false;
    }

    @Override
    public String toString() {
        return "Population{" +
                "encodedSolutions=" + Arrays.toString(encodedSolutions) +
                '}';
    }
}
