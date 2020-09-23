import java.util.ArrayList;
import java.util.List;
import java.util.Random;

public class Main {


    static Random rand = new Random();
    static int numberOfNodes;
    static int numberOfUniqueLinks;
    static int maxConnection;
    static int[] linkLengths;
    static double preserve;

    /**
     * Clever method to get starting population
     * @return the starting population
     */
    private static Population GenerateInitialPopulation(int popSize){
        Population pop = new Population(popSize);

        // Generate and add the genotype encoded solutions to the initial population
        for(int i = 0; i < popSize; i++) {
            Genotype s = GenerateRandomConfiguration();

            // Preform local Search
            s = LocalSearch(s);

            // Add the new solution to the population
            pop.Insert(i, s);
        }

        return pop;
    }

    /**
     * Create a random solution configuration
     * @return a new solution in genotype space
     */
    private static Genotype GenerateRandomConfiguration(){
        Genotype s = new Genotype(numberOfUniqueLinks);

        // Create random solution
        for(int j = 0; j < s.length(); j++){
            // Turn each bit on at random chance of 50%
            if(rand.nextBoolean()) {
                s.flip(j);
            }
        }

        // Make the solution feasible
        s = Repair(s);

        return s;
    }

    /**
     * Performs local search for the best solution until the termination criteria is reached
     * @param currentGeno the current solution in genotype space
     * @return the local best solution in genotype space
     */
    private static Genotype LocalSearch(Genotype currentGeno){
        List<Integer> currentPheno = Growth(currentGeno);
        do{
            Genotype newGeno = GenerateNeighbour(currentGeno);
            newGeno = Repair(newGeno);
            List<Integer> newPheno = Growth(newGeno);
            if(Evaluate(newPheno) < Evaluate(currentPheno)){
                currentGeno = newGeno;
                currentPheno = newPheno;
            }
        } while(!TerminationCriterion());
        return currentGeno;
    }

    /**
     * Generate a random solution in the neighbourhood of the given solution
     * @param current solution in genotype space
     * @return a neighbour solution in genotype space
     */
    private static Genotype GenerateNeighbour(Genotype current){
        int i = rand.nextInt(current.length());
        Genotype neighbour = (Genotype)current.clone();
        neighbour.flip(i);
        return neighbour;
    }

    /**
     * Convert the genotype into an Integer list phenotype
     * @param genotype of the solution
     * @return the phenotype of that solution
     */
    private static List<Integer> Growth(Genotype genotype){
        List<Integer> phenotype = new ArrayList<>();

        // Convert genotype into phenotype
        for(int i = 0; i < genotype.length(); i++){
            if(genotype.isSet(i)){
                phenotype.add(linkLengths[i]);
            }
        }

        return phenotype;
    }

    /**
     * Turn infeasible solutions into feasible ones
     * @param solution in genotype space
     * @return genotype of a feasible solution
     */
    private static Genotype Repair(Genotype solution) {

        int[] connectionOnNode = new int[numberOfNodes];
        int nodeIndex = 0;
        int genotypeIndex = 0;
        // Look through each link making up the nodes
        while(genotypeIndex < solution.length()){
            // Find the range of new links connected to the node at the current nodeIndex
            int range = (numberOfNodes - 1) - nodeIndex;
            Genotype nodeConnections = solution.get(genotypeIndex, genotypeIndex + range);

            // Check each set link for the current node
            for(int i = 0; i < nodeConnections.length(); i++){
                // Set i to the index of the next set bit
                i = nodeConnections.nextSetBit(i);
                // if a set bit was found
                if(i < nodeConnections.length()){
                    // Add one to the current node
                    connectionOnNode[nodeIndex]++;
                    // Add one to the connecting node
                    connectionOnNode[nodeIndex + 1 + i]++;

                    // If either nodes has too many connections then remove the last added
                    if (connectionOnNode[nodeIndex] > maxConnection || connectionOnNode[nodeIndex + 1 + i] > maxConnection) {
                        solution.set(i, false);
                        connectionOnNode[nodeIndex]--;
                        connectionOnNode[nodeIndex + 1 + i]--;
                    }
                }
            }

            //TODO
            // Check that the node has some connection to all other nodes

            // Increment to look at the next node and set of links
            nodeIndex++;
            genotypeIndex += range;
        }
        return solution;
    }

    /**
     * Evaluate the solution
     * @param phenotype space of the solution
     * @return the length of all the links in the solution
     */
    private static int Evaluate(List<Integer> phenotype){
        int total = 0;
        for (int i: phenotype) {
            total += i;
        }
        return total;
    }

    /**
     * Apply recombination and mutation to get new population
     * @param pop the current population
     * @return the new population
     */
    private static Population GenerateNewPopulation(Population pop){
        Population newpop = new Population(pop.Size());
        //TODO
        // Recombination
        // mutation
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
        // Select subset of newpop
        // Combine the subset with pop to get the nextpop
        return nextpop;
    }

    /**
     * Restart the population
     * @param pop the stagnated population
     * @return a new population
     */
    private static Population RestartPopulation(Population pop){
        Population newpop = new Population(pop.Size());
        // preserve is the percent of the solution to keep when restarting
        int numberPreserved = (int)(pop.Size() * preserve);
        for (int i = 0; i < numberPreserved; i++){
            Genotype s = pop.ExtractBest();
            newpop.Insert(i, s);
        }
        for(int i = numberPreserved; i < pop.Size(); i++){
            Genotype s = GenerateRandomConfiguration();
            s = LocalSearch(s);
            newpop.Insert(i, s);
        }
        return newpop;
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
        // Read in these parameters
        int popSize = 2;
        numberOfNodes = 4;
        numberOfUniqueLinks = ( numberOfNodes * (numberOfNodes - 1) ) / 2;
        maxConnection = 2;
        linkLengths = new int[numberOfUniqueLinks];
        for (int i = 0; i < numberOfUniqueLinks; i++) {
            linkLengths[i] = i + 1;
        }

        preserve = 0.2; // This is the only optional parameter that can have a default value

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

        //TODO
        // Output result
    }
}
