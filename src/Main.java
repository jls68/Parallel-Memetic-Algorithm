import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;
import java.util.concurrent.*;

public class Main {

    /**
     * Read in the link lengths from the csv file
     * @param filePath to the csv file
     * @return an int array of each possible links' length
     */
    private static int[] readCSV(String filePath) {
        int[] links = null;

        // Parse CSV file into BufferedReader
        try {
            BufferedReader br = new BufferedReader(new FileReader(filePath));

            String line = br.readLine();
            String[] split = line.split(",");
            links = new int[split.length];
            for (int i = 0; i < split.length; i++) {
                links[i] = Integer.parseInt(split[i]);
            }
        } catch (FileNotFoundException e) {
            //TODO
            e.printStackTrace();
        } catch (IOException e) {
            //TODO
            e.printStackTrace();
        }

        return links;
    }

    public static Genotype compute(Search search) throws InterruptedException {
        search.run();
        //search.join();
        return search.getResult();
    }


    public static void main(String[] args) throws InterruptedException, ExecutionException {
        // Set default parameters
        int n_pr =  Runtime.getRuntime().availableProcessors();
        int popSize = 4;
        int numParents = 2;
        int numberOfNodes;
        int numberOfUniqueLinks;
        int maxConnection;
        int[] linkLengths;
        double preservePercent = 0.2;
        double mutatePercent = 0.2;
        boolean plusInsteadOfComma = true;
        boolean replicatedInsteadOfSynchronous = true;
        long tMax = 1;
        Random rand = new Random();

        try {
            // Read in the parameters
            if (args.length > 1) {
                // First argument should be the filepath of the lengths to connect each node in a specific csv format
                String filePath = args[0];
                linkLengths = readCSV(filePath);
                numberOfNodes = linkLengths.length;
                numberOfUniqueLinks = (numberOfNodes * (numberOfNodes - 1)) / 2;
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
                }
            }
            // If no parameter given then use hard coded test parameters
            else {
                numberOfNodes = 4;
                numberOfUniqueLinks = (numberOfNodes * (numberOfNodes - 1)) / 2;
                linkLengths = new int[numberOfUniqueLinks];
                for (int i = 0; i < numberOfUniqueLinks; i++) {
                    linkLengths[i] = i + 1;
                }
                maxConnection = 2;
            }

            Genotype bestSolution = null;
            List<Integer> solutionPhenotype = null;
            int bestScore = 0;

            // Start Replicated Parallel Search
            if(replicatedInsteadOfSynchronous) {
                Search[] searches = new Search[n_pr];
                for (int i = 0; i < n_pr; i++) {
                    searches[i] = new Search(rand, popSize, numParents, numberOfNodes, numberOfUniqueLinks, maxConnection,
                            linkLengths, mutatePercent, preservePercent, true, tMax);
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
                // some other executors you could try to see the different behaviours
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

            }
            else{
                Search search = new Search(rand, popSize, numParents, numberOfNodes, numberOfUniqueLinks, maxConnection,
                        linkLengths, mutatePercent, preservePercent, plusInsteadOfComma, tMax);
                search.run();
                search.join(0);
                bestSolution = search.getResult();
                solutionPhenotype = search.Growth(bestSolution);
                bestScore = search.Evaluate(solutionPhenotype);
            }




            //TODO
            // Output result
            System.out.println("Genotype: " + bestSolution.toString());
            System.out.println("Phenotype: " + solutionPhenotype.toString());
            System.out.println("Score: " + bestScore);

        } catch (NumberFormatException | InterruptedException | ExecutionException e) {
            //TODO
            e.printStackTrace();
        }
    }
}
