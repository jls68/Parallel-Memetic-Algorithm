import java.util.BitSet;
import java.util.Random;

public class Main {

    static Random rand = new Random();
    static int numberOfNodes;
    static int maxConnection;

    /**
     * Clever method to get starting population
     * @return the starting population
     */
    private static Population GenerateInitialPopulation(int popSize){
        int numberOfUniqueLinks = ( numberOfNodes * (numberOfNodes - 1) ) / 2;

        Population pop = new Population(popSize);

        // Generate and add the genotype encoded solutions to the initial population
        for(int i = 0; i < popSize; i++) {
            BitSet g = new BitSet(numberOfUniqueLinks);

            //TODO
            // Create random solution

            // Preform local Search
            g = LocalSearch(g);

            pop.Insert(i, g);
        }

        return pop;
    }


    private static BitSet LocalSearch(BitSet current){
        do{
            BitSet newSolution = GenerateNeighbour(current);
            if(Evaluate(newSolution) < Evaluate(current)){
                current = newSolution;
            }
        } while(!TerminationCriterion());
        return current;
    }


    private static BitSet GenerateNeighbour(BitSet current){
        int i = rand.nextInt(current.length());
        BitSet neighbour = (BitSet)current.clone();
        neighbour.flip(i);
        return neighbour;
    }

    private static BitSet Repair(BitSet solution) {

        int[] connectedNodes = new int[numberOfNodes];
        int nodeIndex = 0;
        int genotypeIndex = 0;
        // Look through each link making up the nodes
        while(genotypeIndex < solution.length()){
            // Find the range of new links connected to the node at the current nodeIndex
            int range = (numberOfNodes - 1) - nodeIndex;
            BitSet nodeConnections = solution.get(genotypeIndex, genotypeIndex + range);

            // Check each set link for the current node
            for(int i = 0; i < nodeConnections.length(); i++){
                // Set i to the index of the next set bit
                i = nodeConnections.nextSetBit(i);
                // Add one to the current node
                connectedNodes[nodeIndex]++;
                // Add one to the connecting node
                connectedNodes[nodeIndex + 1 + i]++;

                // If either nodes has too many connections then remove the last added
                if(connectedNodes[nodeIndex] > maxConnection || connectedNodes[nodeIndex + 1 + i] > maxConnection) {
                    solution.set(i, false);
                    connectedNodes[nodeIndex]--;
                    connectedNodes[nodeIndex + 1 + i]--;
                }
            }

            //TODO
            // Check that the node is has some connection to all other nodes

            // Increment to look at the next node and set of links
            nodeIndex++;
            genotypeIndex += range;
        }
        return solution;
    }

    private static int Evaluate(BitSet g){
        //TODO
        return 0;
    }

    /**
     * Apply recombination and mutation to get new population
     * @param pop the current population
     * @return the new population
     */
    private static Population GenerateNewPopulation(Population pop){
        Population newpop = new Population(pop.Size());
        //TODO
        return newpop;
    }

    /**
     * Select a subset of newpop for the next pop
     * @param pop the current population
     * @param newpop the new population from recombination, mutation
     * @return a subset for the next population
     */
    private static Population UpdatePopulation(Population pop, Population newpop){
        Population nextpop = new Population(pop.Size());
        //TODO
        return nextpop;
    }

    /**
     * Restart the population
     * @param pop the stagnated population
     * @return a new population
     */
    private static Population RestartPopulation(Population pop){
        Population restartedpop = new Population(pop.Size());
        //TODO
        return restartedpop;
    }

    /**
     * Check if the algorithm needs to stop
     * @return true if the termination criteria is met
     */
    private static boolean TerminationCriterion(){
        //TODO
        return true;
    }

    public static void main(String[] args) {
        //TODO
        int popSize = 2;
        numberOfNodes = 4;
        maxConnection = 2;
        // Initialise starting population
        Population pop = GenerateInitialPopulation(popSize);
        do{
            // Apply recombination, mutation
            Population newpop = GenerateNewPopulation(pop);
            // Select a subset of newpop for the next pop
            pop = UpdatePopulation(pop, newpop);
            // Decide if we need to restart due to stagnation
            if(pop.hasConverged()){
                pop = RestartPopulation(pop);
            }
        } while(!TerminationCriterion());
    }
}
