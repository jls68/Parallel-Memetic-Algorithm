import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.*;
import java.util.concurrent.*;

public class Main {

    private static int IMPOSSIBLE_CONNECTION = 1000;
    private static int RUNS = 10;

    // Global variables to be read in by readCSV method
    private static int numberOfNodes;
    private static int maxConnection;
    private static int numberOfUniqueLinks;
    private static int[] linkLengths;
    private static String[] linkIDs;

    // Global variables to record performance
    private static Genotype[] runSolutions;
    private static int[] runScores;
    private static int[] conversions;
    protected static int runConversions;

    // Global variables to be read from argument
    private static Random rand;
    private static int replicated_n_pr;
    private static int synchronous_n_pr;
    private static int popSize;
    private static int numParents;
    private static int numChildren;
    private static double preservePercent;
    private static double mutatePercent;
    private static boolean sectionInheritance;
    private static boolean plusInsteadOfComma;
    private static double tMax;
    private static int kMax;

    /**
     * Read in the information from the csv file
     * @param filePath to the csv file
     */
    private static void readCSV(String filePath) {
        // Parse CSV file into BufferedReader
        try {
            BufferedReader br = new BufferedReader(new FileReader(filePath));

            numberOfNodes = Integer.parseInt(br.readLine().split(",")[1]);
            numberOfUniqueLinks = (numberOfNodes * (numberOfNodes - 1)) / 2;
            linkIDs = new String[numberOfUniqueLinks];
            linkLengths = new int[numberOfUniqueLinks];

            int i = 0;

            // Read the information for the first possible link
            String line = br.readLine();
            String[] split = line.split(",");
            String[] nodeIDs = split[0].split("-");

            for(int n = 1; n <= numberOfNodes; n++) {
                for (int l = n + 1; l <= numberOfNodes; l++) {

                    int firstNode = Integer.parseInt(nodeIDs[0]);
                    int secondNode = Integer.parseInt(nodeIDs[1]);

                    // Check if we have the information for the next link
                    if (line != "" && firstNode == n && secondNode == l) {
                        // Add the link information
                        linkIDs[i] = split[0];
                        linkLengths[i] = Integer.parseInt(split[1]);

                        // Read in the information on the next line
                        line = br.readLine();
                        if(line != "" && line != null) {
                            split = line.split(",");
                            nodeIDs = split[0].split("-");
                        }
                    }
                    // Else this link is impossible
                    else {
                        // Add impossible link placeholder
                        linkIDs[i] = n + "x" + l;
                        linkLengths[i] = IMPOSSIBLE_CONNECTION;
                    }
                    i++;
                }
            }
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

/**Usage: [filePath] [maxDegree] -S [randomSeed] -Pop [popSize] -parents [parentSize]
 * -preserve [percentagePreserved] -mutate [mutationProbability] {-plus, -comma} {-replicated, -synchronous
 * -tMax [maxRunTime]");
**/

    public static void main(String[] args) throws InterruptedException, ExecutionException {
        // Set default parameters
        int maxProcessors = Runtime.getRuntime().availableProcessors();
        replicated_n_pr = 1;
        synchronous_n_pr = 1;
        popSize = 8;
        numParents = 2;
        numChildren = 4;
        preservePercent = 0.2;
        mutatePercent = 0.2;
        sectionInheritance = false;
        plusInsteadOfComma = true;
        tMax = 1;
        kMax = 1;
        rand = new Random();

        try {
            if (args.length == 0) {
                System.out.println("Usage: [filePath] [maxDegree] -S [randomSeed] -Pop [popSize] -parents [parentSize]" +
                        " -preserve [percentagePreserved] -mutate [mutationProbability] {-plus, -comma} {-replicated, -synchronous} -tMax [maxRunTime]");
            }
            // Read in the parameters
            if (args.length > 1) {
                // First argument should be the filepath of the lengths to connect each node in a specific csv format
                String filePath = args[0];
                readCSV(filePath);

                // Second argument should be the max number of connections for each node
                maxConnection = Integer.parseInt(args[1]);

                if (args.length > 2) {
                    for (int i = 2; i < args.length; i++) {
                        // Optional -S argument followed by an integer is the number for the random seed
                        if (args[i].equals("-S")) {
                            int seed = Integer.parseInt(args[i + 1]);
                            rand = new Random(seed);
                            i++;
                        }
                        // Optional -Pop argument followed by an integer is the population size
                        else if (args[i].equals("-Pop")) {
                            popSize = Integer.parseInt(args[i + 1]);
                            i++;
                        }
                        // Optional -parents argument followed by an integer is the number of parents
                        else if (args[i].equals("-parents")) {
                            numParents = Integer.parseInt(args[i + 1]);
                            i++;
                        }
                        // Optional -children argument followed by an integer is the number of parents
                        else if (args[i].equals("-children")) {
                            numChildren = Integer.parseInt(args[i + 1]);
                            i++;
                        }
                        // Optional -preserve argument followed by a double is the preserve Percent
                        else if (args[i].equals("-preserve")) {
                            preservePercent = Double.parseDouble(args[i + 1]);
                            i++;
                        }
                        // Optional -mutate argument followed by a double is the mutatePercent
                        else if (args[i].equals("-mutate")) {
                            mutatePercent = Double.parseDouble(args[i + 1]);
                            i++;
                        }
                        // Optional argument of a -plus or -comma keyword to select the recombination method
                        else if (args[i].equals("-plus")) {
                            plusInsteadOfComma = true;
                        } else if (args[i].equals("-comma")) {
                            plusInsteadOfComma = false;
                        }
                        // Optional argument of -replicated or -synchronous to select the parallelism method and
                        //  Also optionally followed by an integer to choose the number cores to use instead of the max
                        else if (args[i].equals("-replicated")) {
                            try {
                                replicated_n_pr = Integer.parseInt(args[i + 1]);
                                i++;
                            }
                            catch (Exception e){
                                replicated_n_pr = maxProcessors;
                            }
                        } else if (args[i].equals("-synchronous")) {
                            try {
                            synchronous_n_pr = Integer.parseInt(args[i + 1]);
                            i++;
                            }
                            catch (Exception e){
                                synchronous_n_pr = maxProcessors;
                            }
                        }
                        // Optional argument of -blockInherit to select children to inherit bytes instead of bits, a karger group from each parent
                        else if (args[i].equals("-blockInherit")) {
                            sectionInheritance = true;
                        }
                        // Optional argument of -tMax followed by a long is the max time to search
                        else if (args[i].equals("-tMax")) {
                            tMax = Double.parseDouble(args[i + 1]);
                            i++;
                        }
                        // Optional argument of -kMax followed by a long is the max time to do local search
                        else if (args[i].equals("-kMax")) {
                            kMax = Integer.parseInt(args[i + 1]);
                            i++;
                        }
                    }
                }

                conversions = new int[RUNS];
                runSolutions = new Genotype[RUNS];
                runScores = new int[RUNS];
                int bestSolutionIndex = 0;
                for (int i = 0; i < RUNS; i++) {
                    runConversions = 0;
                    System.out.println("=======================");
                    System.out.println("Run number " + (i + 1) + ":");
                    System.out.println("=======================");
                    runSolutions[i] = doRun();
                    int score = Search.Evaluate(Search.Growth(runSolutions[i]));
                    runScores[i] = score;
                    if(i > 0 && score < runScores[bestSolutionIndex]){
                        bestSolutionIndex = i;
                    }
                    conversions[i] = runConversions;
                }
                System.out.println("============== Done ==============");
                System.out.println("Dataset = " + filePath + ", number of nodes = " + numberOfNodes + ", max degree = " + maxConnection +
                        ", Replicated = " + replicated_n_pr + ", Synchronous = " + synchronous_n_pr + ", pop size = " + popSize + ", parents = " + numParents +
                        ", children = " + numChildren + ", preserve = " + preservePercent + ", mutate = " + mutatePercent + ", section = " +
                        sectionInheritance + ", plus = " + plusInsteadOfComma + ", tMax = " + tMax + ", kmax = " + kMax);
                System.out.println(RUNS + " run average scores:\t" + average(runScores) + "\t" + Arrays.toString(runScores));
                System.out.println("best solution\t" + runScores[bestSolutionIndex] + "\t" + runSolutions[bestSolutionIndex]);
                System.out.println("best solution lengths\t " + runSolutions[bestSolutionIndex].idLinks(linkIDs));
                System.out.println("max number of conversions\t" + max(conversions));
                System.out.println("avg number of conversions\t" + average(conversions) + "\t" + Arrays.toString(conversions));
            }

        } catch (NumberFormatException | InterruptedException | ExecutionException e) {
            e.printStackTrace();
        }
    }

    private static double average(int[] array){
        double total = 0;
        for(int i=0; i<array.length; i++){
            total = total + array[i];
        }
        return total / array.length;
    }

    private static int max(int[] array){
        int max = array[0];
        for(int i=1; i<array.length; i++){
            if(array[i] > max){
                max = array[i];
            }
        }
        return max;
    }

    public static Genotype doRun()
            throws InterruptedException, ExecutionException {
        // Start the timer to record how long it takes to search for the best solution.
        final long initialTime = System.nanoTime();

        Genotype bestSolution = null;
        List<Integer> solutionPhenotype = null;
        int bestScore = 0;

        Search[] searches = new Search[replicated_n_pr];
        for (int i = 0; i < replicated_n_pr; i++) {
            searches[i] = new Search(rand, popSize, numParents, numChildren, numberOfNodes, numberOfUniqueLinks, maxConnection, synchronous_n_pr,
                    linkLengths, mutatePercent, preservePercent, sectionInheritance, plusInsteadOfComma, tMax, kMax);
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

        ExecutorService exec = Executors.newFixedThreadPool(replicated_n_pr);
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
        final double secondsToComplete = (System.nanoTime() - initialTime) / 1E9;

        System.out.println("Completed in " + secondsToComplete + " seconds");
        System.out.println("Genotype: " + bestSolution.toString());
        System.out.println("Phenotype: " + solutionPhenotype.toString());
        System.out.println(bestSolution.idLinks(linkIDs));
        System.out.println("Score: " + bestScore);

        return bestSolution;
    }
}
