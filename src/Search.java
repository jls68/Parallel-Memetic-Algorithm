import java.util.ArrayList;
import java.util.List;
import java.util.Random;
import java.util.concurrent.*;

public class Search extends Thread{

    static Random rand;
    int convergeAmount;
    int searchAmount;
    int popSize;
    int numParents;
    int numChildren;
    static int numberOfNodes;
    int numberOfUniqueLinks;
    static int maxConnection;
    int n_pr;
    static int[] linkLengths;
    double mutatePercent;
    double preservePercent;
    boolean sectionInheritance;
    boolean plusInsteadOfComma;
    Genotype bestSolution;
    double tMax;
    int kMax;
    private long initialRunTime;
    boolean VND;


    Search(Random rand, int popSize, int numParents, int numChildren, int numberOfNodes, int numberOfUniqueLinks, int maxConnection, int VNDn_pr, int[] linkLengths,
           double mutatePercent, double preservePercent, boolean sectionInheritance, boolean plusInsteadOfComma, double tMax, int kMax, boolean VND){
        convergeAmount = 0;
        searchAmount = 0;
        this.rand = rand;
        this.popSize = popSize;
        this.numParents = numParents;
        this.numChildren = numChildren;
        this.numberOfNodes = numberOfNodes;
        this.numberOfUniqueLinks = numberOfUniqueLinks;
        this.maxConnection = maxConnection;
        n_pr = VNDn_pr;
        this.linkLengths = linkLengths;
        this.mutatePercent = mutatePercent;
        this.preservePercent = preservePercent;
        this.sectionInheritance = sectionInheritance;
        this.plusInsteadOfComma = plusInsteadOfComma;
        this.tMax = tMax;
        this.kMax = kMax;
        this.VND = VND;
    }

    synchronized List<Genotype> combine (List<Genotype> inputList, Genotype newvalue){
        inputList.add(newvalue);
        return inputList;
    }

    public void run() {
        initialRunTime = System.nanoTime();
        // Initialise starting population
        Population pop = null;
        try {
            pop = GenerateInitialPopulation(popSize);
            do {
                // Apply recombination, mutation
                Population newpop = GenerateNewPopulation(pop, numParents, numChildren, mutatePercent, sectionInheritance);
                // Select a subset of newpop for the next pop
                if (plusInsteadOfComma) {
                    pop = UpdatePopulationPlus(pop, newpop);
                } else {
                    pop = UpdatePopulationComma(pop, newpop);
                }
                // Decide if we need to restart due to stagnation
                if (pop.hasConverged()) {
                    pop = RestartPopulation(pop, preservePercent);
                    convergeAmount++;
                }
            } while (!TerminationCriterion(tMax));
            // Find best solution in population
            List<Integer>[] popPheno = ConvertPopToPhenotype(pop);
            int index = FindBestSolution(popPheno);
            // Store best solution
            bestSolution = pop.getSolution(index);

            //finish timing program
            long finalTime = System.nanoTime();
            //Please do not remove or change the format of this output message
            System.out.println("Thread  " + this.getId() + " finished execution in " + (finalTime - initialRunTime) / 1E9 + " secs. Converged " + convergeAmount + " times.");

            addConverge(convergeAmount);
            addSearches(searchAmount);
        } catch (InterruptedException e) {
            e.printStackTrace();
        }
    }

    private synchronized void addSearches(int searchAmount) {
        Main.runSearches += searchAmount;
    }

    private synchronized void addConverge(int convergeAmount){
        Main.runConversions += convergeAmount;
    }

    /**
     * Clever method to get starting population
     * @return the starting population
     */
    private Population GenerateInitialPopulation(int popSize) throws InterruptedException {
        Population pop = new Population(popSize);

        // Generate and add the genotype encoded solutions to the initial population
        for(int i = 0; i < popSize; i++) {
            Genotype s = GenerateRandomConfiguration();

            // Preform local Search
            if(VND) {
                s = VND(s);
            }
            else {
                s = LocalSearch(s);
            }

            // Add the new solution to the population
            pop.Insert(i, s);
        }

        return pop;
    }

    /**
     * Create a random solution configuration
     * @return a new solution in genotype space
     */
    private Genotype GenerateRandomConfiguration(){
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
    private Genotype VND(Genotype currentGeno) throws InterruptedException {
        Genotype xBest = currentGeno;
        int bestScore = Evaluate(Growth(xBest));

        SynchronousSearch[] SynchronousSearch = new SynchronousSearch[n_pr];
        for (int j = 0; j < n_pr; j++) {
            SynchronousSearch[j] = new SynchronousSearch(rand, currentGeno, kMax);
        }

        int i = 0;
        while (i < linkLengths.length) {
            for (int j = 0; j < n_pr && i < linkLengths.length; j++) {
                SynchronousSearch[j].addIndex(i);
                i++;
            }
        }

        List<Callable<Genotype>> tasks = new ArrayList<Callable<Genotype>>();
        for (final SynchronousSearch search : SynchronousSearch) {
            Callable<Genotype> callable = new Callable<Genotype>() {
                @Override
                public Genotype call() throws Exception {
                    return compute(search);
                }
            };
            tasks.add(callable);
        }

        ExecutorService exec = Executors.newFixedThreadPool(n_pr);
        // some other executors we could try to see the different behaviours
        // ExecutorService exec = Executors.newCachedThreadPool();
        // ExecutorService exec = Executors.newSingleThreadExecutor();

        try {
            // Start all searches to run in parallel
            List<Future<Genotype>> results = exec.invokeAll(tasks);
            // Find best solution of all searches
            for (Future<Genotype> search : results) {

                // Wait until all threads have finished execution and get best solution
                Genotype newSolution = search.get();
                // Convert to phenotype space to then evaluate
                List<Integer> newPheno = Search.Growth(newSolution);
                int newScore = Search.Evaluate(newPheno);
                // Keep track of the solution with the best score
                if (newScore < bestScore) {
                    xBest = newSolution;
                    bestScore = newScore;
                }
            }
        } catch (ExecutionException e) {
            e.printStackTrace();
        } finally {
            exec.shutdown();
        }

        searchAmount++;
        return xBest;
    }

    private static Genotype compute(SynchronousSearch search) throws InterruptedException {
        search.run();
        return search.getResult();
    }

    private Genotype LocalSearch(Genotype currentGeno){
        List<Integer> currentPheno = Growth(currentGeno);
        do{
            Genotype newGeno = GenerateNeighbour(currentGeno);
            newGeno = Repair(newGeno);
            List<Integer> newPheno = Growth(newGeno);
            if(IsBetterThan(newPheno, currentPheno)){
                currentGeno = newGeno;
                    currentPheno = newPheno;
                }
        } while(!TerminationCriterion(tMax / 100));
        searchAmount++;
        return currentGeno;
    }

    /**
     * Generate a random solution in the neighbourhood of the given solution, does the same as applying a mutation
     * @param current solution in genotype space
     * @return a neighbour solution in genotype space
     */
    private Genotype GenerateNeighbour(Genotype current){
        int i = rand.nextInt(current.length());
        current = current.clone();
        current.flip(i);
        return current;
    }

    /**
     * Convert the genotype into an Integer list phenotype
     * @param genotype of the solution
     * @return the phenotype of that solution
     */
    static List<Integer> Growth(Genotype genotype){
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
    protected static Genotype Repair(Genotype solution) {

        Node[] nodes = new Node[numberOfNodes];
        for (int i = 0; i < numberOfNodes; i++) {
            nodes[i] = new Node(i);
        }

        int nodeIndex = 0;
        int genotypeIndex = 0;

        // Look through each link making up the nodes
        while(nodeIndex < numberOfNodes){
            // Find the range of new links connected to the node at the current nodeIndex
            int range = (numberOfNodes - 1) - nodeIndex;

            // Check each set link for the current node
            for(int i = solution.nextSetBit(genotypeIndex); i < genotypeIndex + range; i = solution.nextSetBit(i + 1)) {
                // Set bit was found
                Node thisNode = nodes[nodeIndex];
                Node otherNode = nodes[nodeIndex + 1 + i - genotypeIndex];

                // Record the link between the two nodes
                new Link(thisNode, otherNode);

                // If either nodes has too many connections then remove a random connection from that node
                if (thisNode.exceedLimit(maxConnection)) {
                    thisNode.removeRandomLink(rand);
                }
                // Check the other node in case it still has too many connections
                if (otherNode.exceedLimit(maxConnection)) {
                    otherNode.removeRandomLink(rand);
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
            for (int i = 0; i < startIndex && !infeasibleSolution; i++) {
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
    static int Evaluate(List<Integer> phenotype){
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
    private Boolean IsBetterThan(List<Integer> isPhenotype, List<Integer> thanPhenotype){
        return Evaluate(isPhenotype) < Evaluate(thanPhenotype);
    }

    /**
     * Apply Recombination that has implicit mutation to get a new population
     * @param pop the current population
     * @return the new population
     */
    private Population GenerateNewPopulation(Population pop, int numParents, int numChildren, double mutationPercent, boolean sectionInheritance) {
        Population newpop = new Population(popSize * (numChildren / numParents));

        // Integer list of indexes of solutions in the pop that have not been selected as parents yet
        List<Integer> indexes = new ArrayList<>();
        for (int i = 0; i < popSize; i++) {
            indexes.add(i);
        }

        int childIndex = 0;

        // Group all in solutions in the current population to form the groups of parents
        while (numParents <= indexes.size()) {

            // Grab a group of parent solutions
            Genotype[] parents = new Genotype[numParents];
            for (int p = 0; p < numParents; p++) {
                int parentIndex = rand.nextInt(indexes.size());
                parents[p] = pop.getSolution(indexes.get(parentIndex));
                indexes.remove(parentIndex);
            }

            // Recombine the parents to get the children solutions
            for (int c = 0; c < numChildren; c++) {

                // Initialise the child solution
                Genotype newChild = new Genotype(numberOfUniqueLinks);

                // Take a byte from each parent
                if(sectionInheritance){
                    int numberOfSections = numberOfUniqueLinks / 8;
                    for (int s = 0; s < numberOfSections; s++){
                        Byte section = parents[rand.nextInt(numParents)].getByte(s);
                        newChild.set(s, section);
                    }
                }
                // Fill the child using bits from the parents solution in the current population, pop
                else {
                    // Inherit each bit from a random parent
                    for (int b = 0; b < numberOfUniqueLinks; b++) {
                        newChild.set(b, parents[rand.nextInt(numParents)].isSet(b));
                    }
                }

                // Add mutation at random
                if (rand.nextDouble() < mutationPercent) {
                    newChild = GenerateNeighbour(newChild);
                }

                // Repair solution
                newChild = Repair(newChild);

                // Add child to the new population
                newpop.Insert(childIndex, newChild);
                childIndex++;
            }
        }
        // Add the remainder parent to the new population
        while (childIndex < newpop.Size()){
            for(int i = 0; i < indexes.size() && childIndex < newpop.Size(); i++){
                newpop.Insert(childIndex, pop.getSolution(indexes.get(i)).clone());
                childIndex++;
            }
        }

        return newpop;
    }

    /**
     * Select a subset of newpop for the next pop using the plus strategy
     * @param pop the current population
     * @param newpop the new population from recombination, mutation
     * @return a subset for the next population
     */
    private Population UpdatePopulationPlus(Population pop, Population newpop){
        Population nextpop = new Population(popSize);
        List<Integer>[] currPheno = ConvertPopToPhenotype(pop);
        List<Integer>[] newPheno = ConvertPopToPhenotype(newpop);
        // Select subset of newpop
        // In the plus strategy, we add pop and newpop together and select the best popsize solutions for the final version of nextpop
        int currBest = FindBestSolution(currPheno);
        int newBest = FindBestSolution(newPheno);
        for(int i = 0; i < popSize; i++){
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
    private Population UpdatePopulationComma(Population pop, Population newpop){
        Population nextpop = new Population(popSize);
        List<Integer>[] newPheno = ConvertPopToPhenotype(newpop);
        // Select subset of newpop
        for(int i = 0; i < popSize; i++) {
            // we select the best popsize elements from newpop
            int newBest = FindBestSolution(newPheno);
            // Add the best solutions to the next population
            nextpop.Insert(i, newpop.getSolution(newBest));
            // Get the next best solution from the new population
            newPheno[newBest] = null;
        }
        return nextpop;
    }

    /**
     * Convert Population from genotype space to phenotype space
     * @param pop
     * @return array of phenotype solutions
     */
    private List<Integer>[] ConvertPopToPhenotype(Population pop) {
        List<Integer>[] phenotypes = new List[popSize];
        for (int i = 0; i < popSize; i++){
            phenotypes[i] = Growth(pop.getSolution(i));
        }
        return phenotypes;
    }

    /**
     * Find the best solution
     * @param phenotypes
     * @return the index of solution that evaluates to the smallest score
     */
    private int FindBestSolution(List<Integer>[] phenotypes){
        int best = 0;
        for (int i = 1; i < phenotypes.length; i++){
            if(phenotypes[i] != null && (phenotypes[best] == null || IsBetterThan(phenotypes[i], phenotypes[best]))){
                best = i;
            }
        }
        return best;
    }

    /**
     * Restart the population
     * @param pop the stagnated population
     * @param preserve percent for a solutions from the current population to be put into thr nre population
     * @return a new population
     */
    private Population RestartPopulation(Population pop, double preserve) throws InterruptedException {
        Population newpop = new Population(popSize);
        List<Integer>[] currPheno = ConvertPopToPhenotype(pop);
        // preserve is the percent of the solution to keep when restarting
        int numberPreserved = (int)(popSize * preserve);
        for (int i = 0; i < numberPreserved; i++){
            int index = FindBestSolution(currPheno);
            Genotype s = pop.getSolution(index);
            newpop.Insert(i, s);
            currPheno[index] = null;
        }
        for(int i = numberPreserved; i < popSize; i++){
            Genotype s = GenerateRandomConfiguration();
            if(VND) {
                s = VND(s);
            }
            else{
                s = LocalSearch(s);
            }
            newpop.Insert(i, s);
        }
        return newpop;
    }

    /**
     * Check if the algorithm needs to stop
     * @param maxTime the max amount of time for the search to go for
     * @return true if the max time has elapsed or there have been too many convergences
     */
    private boolean TerminationCriterion(double maxTime){
        final long finalTime = System.nanoTime();
       if( (finalTime - initialRunTime) / 1E9 > maxTime) {
           return true;
       }
       else {
           return false;
       }
    }

    public Genotype getResult() {
        return bestSolution;
    }
}
