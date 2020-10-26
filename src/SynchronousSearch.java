import java.util.ArrayList;
import java.util.List;
import java.util.Random;

public class SynchronousSearch extends Thread {

    Random rand;
    Genotype bestSolution;
    int kMax;
    List<Integer> searchArea;

    public SynchronousSearch(Random rand, Genotype currentGeno, int kMax){
        this.rand = rand;
        bestSolution = currentGeno;
        this.kMax = kMax;
        searchArea = new ArrayList<>();
    }

    public void addIndex(int value){
        searchArea.add(value);
    }

    /**
     * Performs local search for the best solution until the termination criteria is reached
     */
    public void run(){
        int bestScore = Search.Evaluate(Search.Growth(bestSolution));
        int k = 1;
        do{
            Genotype newGeno = getBestInNeighborhood(bestSolution, k);
            newGeno = Search.Repair(newGeno);
            int newScore = Search.Evaluate(Search.Growth(newGeno));

            if(newScore < bestScore){
                bestSolution = newGeno;
                bestScore = newScore;
                k = 1;
            }
            else{
                k++;
            }
        } while(k < kMax);
    }

    public Genotype getBestInNeighborhood(Genotype currentGeno, int k) {
        Genotype xBest = null;
        int bestScore = 1000;

        // Create neighbours of solution that have k difference
        for (int i : searchArea) {
            Genotype xNew = GenerateNeighbour(currentGeno, i, k);
            int newScore = Search.Evaluate(Search.Growth(xNew));

            if (xBest == null || newScore < bestScore) {
                xBest = xNew;
                bestScore = newScore;
            }
        }

        return xBest;
    }

    /**
     * Generate a new solution in the neighbourhood of the given solution
     * @param current solution in genotype space
     * @param i the index of the first bit to flip
     * @param k the number of other bits to also flip
     * @return a neighbour solution in genotype space
     */
    private Genotype GenerateNeighbour(Genotype current, int i, int k){
        current = current.clone();
        current.flip(i);
        for(int j = 1; j < k; j++){
            int r = rand.nextInt(current.length());
            current.flip(r);
        }
        return current;
    }

    public Genotype getResult() {
        return bestSolution;
    }
}
