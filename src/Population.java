import java.util.BitSet;

public class Population {

    BitSet[] encodedSolutions;

    public Population(int popSize){
        encodedSolutions = new BitSet[popSize];
    }

    public int Size() {
        return encodedSolutions.length;
    }

    public void Insert(int index, BitSet newSolution){
        encodedSolutions[index] = newSolution;
    }

    public boolean hasConverged(){
        //TODO
        // Check for convergence
        return false;
    }

    public BitSet ExtractBest(){
        BitSet best = encodedSolutions[0];
        //TODO
        // Search for best solution
        return best;
    }
}
