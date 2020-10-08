import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;
import java.util.concurrent.*;

public class Main {

    // Global variables to be read in by readCSV method
    private static int numberOfNodes;
    private static int numberOfUniqueLinks;
    private static int[] linkLengths;
    private static String[] linkIDs;

    /**
     * Read in the information from the csv file
     * @param filePath to the csv file
     */
    private static void readCSV(String filePath) {
        // Parse CSV file into BufferedReader
        try {
            BufferedReader br = new BufferedReader(new FileReader(filePath));

            String line = br.readLine();
            String[] split = line.split(",");
            numberOfNodes = Integer.parseInt(split[3]);
            numberOfUniqueLinks = (numberOfNodes * (numberOfNodes - 1)) / 2;
            linkIDs = new String[numberOfUniqueLinks];
            linkLengths = new int[numberOfUniqueLinks];
            int i = 0;

            while ("" != (line = br.readLine()) && i < numberOfUniqueLinks){
                split = line.split(",");
                linkIDs[i] = split[0];
                linkLengths[i] = Integer.parseInt(split[1]);
                i++;
            }
            //TODO
            // Check that all the link weights have been given

        } catch (FileNotFoundException e) {
            //TODO
            e.printStackTrace();
        } catch (IOException e) {
            //TODO
            e.printStackTrace();
        }
    }

    public static Genotype compute(Search search) throws InterruptedException {
        search.run();
        return search.getResult();
    }


    public static void main(String[] args) throws InterruptedException, ExecutionException {
        // Set default parameters
        int n_pr =  Runtime.getRuntime().availableProcessors();
        int popSize = 4;
        int numParents = 2;
        int numChildren = 4; // TODO Allow the number of children to be set
        int maxConnection;
        double preservePercent = 0.2;
        double mutatePercent = 0.2;
        boolean plusInsteadOfComma = true;
        boolean replicatedInsteadOfSynchronous = true;
        long tMax = 1;
        long localtMax = (long) 0.01;
        int kMax = 10;
        Random rand = new Random();

        try {
            // Read in the parameters
            if (args.length > 1) {
                // First argument should be the filepath of the lengths to connect each node in a specific csv format
                String filePath = args[0];
                readCSV(filePath);

                // Second argument should be the max number of connections for each node
                maxConnection = Integer.parseInt(args[1]);

                if (args.length > 2) {
                    //TODO
                    // Optional -Pop argument followed by an integer is the population size
                    // Optional -parents argument followed by an integer is the number of parents
                    // Optional -preserve argument followed by a double is the preserve Percent
                    // Optional -mutate argument followed by a double is the mutatePercent
                    // Optional argument of a plus or comma keyword to select the recombination method
                    // Optional argument of -replicated or -synchronous to select the parallelism method and
                    //  followed by an integer to choose the number cores to use
                    // Optional argument of -kMax followed by a long is the max time to search
                    // Optional argument of -tMax followed by a long is the max time to search
                    // Optional argument of -localtMax followed by a long is the max time to do local search
                }
                // Start the timer to record how long it takes to search for the best solution.
                final long initialTime = System.nanoTime();

                Genotype bestSolution = null;
                List<Integer> solutionPhenotype = null;
                int bestScore = 0;

                // Start Replicated Parallel Search
                if (replicatedInsteadOfSynchronous) {
                    Search[] searches = new Search[n_pr];
                    for (int i = 0; i < n_pr; i++) {
                        searches[i] = new Search(rand, popSize, numParents, numChildren, numberOfNodes, numberOfUniqueLinks, maxConnection,
                                linkLengths, mutatePercent, preservePercent, plusInsteadOfComma, tMax, localtMax, kMax);
                    }

                    List<Callable<Genotype>> tasks = new ArrayList<Callable<Genotype>>();
                    for (final Search search : searches) {
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
                            if (bestSolution == null || newScore < bestScore) {
                                bestSolution = newSolution;
                                solutionPhenotype = newPheno;
                                bestScore = newScore;
                            }
                        }
                    } finally {
                        exec.shutdown();
                    }

                } else {
                    Search search = new Search(rand, popSize, numParents, numChildren, numberOfNodes, numberOfUniqueLinks, maxConnection,
                            linkLengths, mutatePercent, preservePercent, plusInsteadOfComma, tMax, localtMax, kMax);
                    search.run();
                    search.join(0);
                    bestSolution = search.getResult();
                    solutionPhenotype = search.Growth(bestSolution);
                    bestScore = search.Evaluate(solutionPhenotype);
                }
                final double secondsToComplete = (System.nanoTime() - initialTime) / 1E9;

                System.out.println("Completed in " + secondsToComplete + " seconds");
                System.out.println("Genotype: " + bestSolution.toString());
                System.out.println("Phenotype: " + solutionPhenotype.toString());
                System.out.println(bestSolution.idLinks(linkIDs));
                System.out.println("Score: " + bestScore);
            }

        } catch (NumberFormatException | InterruptedException | ExecutionException e) {
            //TODO
            e.printStackTrace();
        }
    }
}
