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
        return false;
    }
}
