public class Main {

    /**
     * Clever method to get starting population
     * @return the starting population
     */
    private static Population GenerateInitialPopulation(){
        Population pop = new Population();
        //TODO
        return pop;
    }

    /**
     * Apply recombination and mutation to get new population
     * @param pop the current population
     * @return the new population
     */
    private static Population GenerateNewPopulation(Population pop){
        Population newpop = new Population();
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
        Population nextpop = new Population();
        //TODO
        return nextpop;
    }

    /**
     * Restart the population
     * @param pop the stagnated population
     * @return a new population
     */
    private static Population RestartPopulation(Population pop){
        Population restartedpop = new Population();
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
        // Initialise starting population
        Population pop = GenerateInitialPopulation();
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
