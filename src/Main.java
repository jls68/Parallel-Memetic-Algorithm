import javax.print.DocFlavor;
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
            if(IsBetterThan(newPheno, currentPheno)){
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

        Node[] nodes = new Node[numberOfNodes];
        for (int i = 0; i < numberOfNodes; i++) {
            nodes[i] = new Node(i);
        }

        int nodeIndex = 0;
        int genotypeIndex = 0;

        // Look through each link making up the nodes
        while(genotypeIndex < solution.length()){
            // Find the range of new links connected to the node at the current nodeIndex
            int range = (numberOfNodes - 1) - nodeIndex;
            Genotype nodeConnections = solution.getSubset(genotypeIndex, genotypeIndex + range);

            // Check each set link for the current node
            for(int i = 0; i < nodeConnections.length(); i++){
                // Set i to the index of the next set bit
                i = nodeConnections.nextSetBit(i);
                // if a set bit was found
                if(i < nodeConnections.length()){
                    // Record the link between the two nodes
                    new Link(nodes[nodeIndex], nodes[nodeIndex + 1 + i]);

                    // If either nodes has too many connections then remove a random connection from that node
                    if (nodes[nodeIndex].exceedLimit(maxConnection)) {
                        nodes[nodeIndex].removeRandomLink(rand);
                    }
                    // Check the other node in case it still has too many connections
                    if (nodes[nodeIndex + 1 + i].exceedLimit(maxConnection)){
                        nodes[nodeIndex + 1 + i].removeRandomLink(rand);
                    }
                }
            }

            // Increment to look at the next node and set of links
            nodeIndex++;
            genotypeIndex += range;
        }

        // Check that each node has some connection to all other nodes
        boolean infeasibleSolution = true;
        while(infeasibleSolution) {
            infeasibleSolution = false;
            Node[] connectedNodes = nodes[0].checkNodesConnected(new Node[numberOfNodes]);
            // Start checking at a random spot to allow random nodes to be repaired
            int startIndex = rand.nextInt(connectedNodes.length);
            // Check second partition
            for (int i = startIndex; i < connectedNodes.length; i++) {
                if (connectedNodes[i] == null) {
                    infeasibleSolution = true;
                    nodes[i].addNewRandomLink(connectedNodes, rand, maxConnection);
                    break; // Break to check if the change has connected all other nodes
                }
            }
            // Check first partition if no node has been found to not be connected
            for (int i = 0; i < startIndex && infeasibleSolution == false; i++) {
                if (connectedNodes[i] == null) {
                    infeasibleSolution = true;
                    nodes[i].addNewRandomLink(connectedNodes, rand, maxConnection);
                    break; // Break to check if the change has connected all other nodes
                }
            }
        }

        // Update solution
        solution = new Genotype(solution.length());
        genotypeIndex = 0;
        for(nodeIndex = 0; nodeIndex < numberOfNodes; nodeIndex++){
            List<Integer> links = nodes[nodeIndex].getNodesLinked();
            for (int l: links) {
                solution.set(genotypeIndex + l, true);
            }
            genotypeIndex += numberOfNodes - (nodeIndex + 1);
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
     * Compare two solutions in phenotype space
     * @param isPhenotype to check if better
     * @param thanPhenotype to compare against
     * @return true if the first solution is better
     */
    private static Boolean IsBetterThan(List<Integer> isPhenotype, List<Integer> thanPhenotype){
        return Evaluate(isPhenotype) < Evaluate(thanPhenotype);
    }

    /**
     * Apply Recombination that has implicit mutation to get a new population
     * @param pop the current population
     * @return the new population
     */
    private static Population GenerateNewPopulation(Population pop, int numParents, double mutationPercent){
        Population newpop = new Population(pop.Size());

        // For each child to be added to the new population
        for (int childIndex = 0; childIndex < pop.Size(); childIndex++) {

            Genotype[] parents = new Genotype[numParents];
            // Get random parents
            for (int i = 0; i < numParents; i++) {
                parents[i] = pop.getSolution(rand.nextInt(pop.Size()));
            }

            int i = 0;
            // Get solution genotype length
            int genoLength = pop.getSolution(0).length();
            // Initialise the child solution
            Genotype newChild = new Genotype(genoLength);

            // Fill the new population, the child, using bits from the parents solution in the current population, pop
            while(i < genoLength){
                // Recombination of each parent in turn
                for(int j = 0; j < numParents && i < genoLength; j++){
                    newChild.set(i, parents[j].getBit(i));
                    i++;
                }
            }

            // Add mutation at random
            if(rand.nextDouble() < mutationPercent){
                newChild = Mutate(newChild);
            }

            // Repair solution
            newChild = Repair(newChild);

            // Add child to the new population
            newpop.Insert(childIndex, newChild);
        }

        return newpop;
    }

    /**
     * Mutate solution
     * @param solution
     * @return a solution with a random bit flipped
     */
    private static Genotype Mutate(Genotype solution){
        int i = rand.nextInt(solution.length());
        solution.flip(i);
        return solution;
    }

    /**
     * Select a subset of newpop for the next pop using the plus strategy
     * @param pop the current population
     * @param newpop the new population from recombination, mutation
     * @return a subset for the next population
     */
    private static Population UpdatePopulationPlus(Population pop, Population newpop){
        Population nextpop = new Population(pop.Size());
        List<Integer>[] currPheno = ConvertPopToPhenotype(pop);
        List<Integer>[] newPheno = ConvertPopToPhenotype(newpop);
        // Select subset of newpop
        // In the plus strategy, we add pop and newpop together and select the best popsize solutions for the final version of nextpop
        int currBest = FindBestSolution(currPheno);
        int newBest = FindBestSolution(newPheno);
        for(int i = 0; i < nextpop.Size(); i++){
            // Add the better solution to the next population
            if(IsBetterThan(currPheno[currBest], newPheno[newBest])){
                nextpop.Insert(i, pop.getSolution(currBest));
                // Get the next best solution from the current population
                currPheno[currBest] = null;
                currBest = FindBestSolution(currPheno);
            }
            else {
                nextpop.Insert(i, newpop.getSolution(newBest));
                // Get the next best solution from the new population
                newPheno[newBest] = null;
                newBest = FindBestSolution(newPheno);
            }
        }

        // Combine the subset with pop to get the nextpop
        return nextpop;
    }

    /**
     * Select a subset of newpop for the next pop using the comma strategy
     * this works when the child population is larger than the population size
     * @param pop the current population
     * @param newpop the new population from recombination, mutation
     * @return a subset for the next population
     */
    private static Population UpdatePopulationComma(Population pop, Population newpop){
        Population nextpop = new Population(pop.Size());
        List<Integer>[] newPheno = ConvertPopToPhenotype(newpop);
        // Select subset of newpop
        // we select the best popsize elements from newpop
        int newBest = FindBestSolution(newPheno);
        for(int i = 0; i < nextpop.Size(); i++) {
            // Add the best solutions to the next population
            nextpop.Insert(i, newpop.getSolution(newBest));
            // Get the next best solution from the new population
            newPheno[newBest] = null;
            newBest = FindBestSolution(newPheno);
        }
        return nextpop;
    }

    /**
     * Convert Population from genotype space to phenotype space
     * @param pop
     * @return array of phenotype solutions
     */
    private static List<Integer>[] ConvertPopToPhenotype(Population pop) {
        List<Integer>[] phenotypes = new List[pop.Size()];
        for (int i = 0; i < pop.Size(); i++){
            phenotypes[i] = Growth(pop.getSolution(i));
        }
        return phenotypes;
    }

    /**
     * Find the best solution
     * @param phenotypes
     * @return the index of solution that evaluates to the smallest score
     */
    private static int FindBestSolution(List<Integer>[] phenotypes){
        int best = 0;
        for (int i = 1; i < phenotypes.length; i++){
            if(phenotypes[i] != null && (phenotypes[best] == null ||Evaluate(phenotypes[i]) < Evaluate(phenotypes[best]))){
                best = i;
            }
        }
        return best;
    }

    /**
     * Restart the population
     * @param pop the stagnated population
     * @return a new population
     */
    private static Population RestartPopulation(Population pop){
        Population newpop = new Population(pop.Size());
        List<Integer>[] currPheno = ConvertPopToPhenotype(pop);
        // preserve is the percent of the solution to keep when restarting
        int numberPreserved = (int)(pop.Size() * preserve);
        for (int i = 0; i < numberPreserved; i++){
            int index = FindBestSolution(currPheno);
            Genotype s = pop.getSolution(index);
            newpop.Insert(i, s);
            currPheno[index] = null;
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
        int popSize = 4;
        int numParents = 2;
        double mutatePercent = 0.2;
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
            Population newpop = GenerateNewPopulation(pop, numParents, mutatePercent);
            // Select a subset of newpop for the next pop
            pop = UpdatePopulationPlus(pop, newpop);
            // Decide if we need to restart due to stagnation
            if(pop.hasConverged()){
                pop = RestartPopulation(pop);
            }
        } while(!TerminationCriterion());

        //TODO
        // Output result
        List<Integer>[] popPheno = ConvertPopToPhenotype(pop);
        int bestSolution = FindBestSolution(popPheno);
        System.out.println("Genotype: " + pop.getSolution(bestSolution).toString());
        System.out.println("Phenotype: " + popPheno[bestSolution].toString());
        System.out.println("Score: " + Evaluate(popPheno[bestSolution]));
    }
}
