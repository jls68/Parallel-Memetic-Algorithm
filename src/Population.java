import java.util.BitSet;

public class Population {

    Genotype[] encodedSolutions;

    public Population(int popSize){
        encodedSolutions = new Genotype[popSize];
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

    public Genotype ExtractBest(){
        Genotype best = encodedSolutions[0];
        //TODO
        // Search for best solution
        return best;
    }
}
