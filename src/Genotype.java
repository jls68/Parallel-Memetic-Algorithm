import java.util.Arrays;
import java.util.List;

public class Genotype {

    boolean[] encodedSolution;

    public Genotype(int length){
        encodedSolution = new boolean[length];
    }

    public int length(){
        return encodedSolution.length;
    }

    public void flip(int index){
        encodedSolution[index] = !encodedSolution[index];
    }

    public boolean isSet(int index){
        return encodedSolution[index] == true;
    }

    /**
     * Get a bit of the genotype
     * @param index
     * @return
     */
    public boolean getBit(int index) {
        return encodedSolution[index];
    }

    /**
     * Get a subset of the genotype
     * @param startIndex inclusive
     * @param endIndex exclusive
     * @return
     */
    public Genotype getSubset(int startIndex, int endIndex) {
        Genotype subset = new Genotype(endIndex - startIndex);
        int i = 0;
        for(int j = startIndex; j < endIndex; j++){
            subset.set(i, encodedSolution[j]);
            i++;
        }
        return subset;
    }

    public void set(int index, boolean value){
        encodedSolution[index] = value;
    }

    public Genotype clone(){
        Genotype clone = new Genotype(encodedSolution.length);
        for(int i = 0; i < encodedSolution.length; i++){
            clone.set(i, encodedSolution[i]);
        }
        return clone;
    }

    public int nextSetBit(int i) {
        for(i = i; i < encodedSolution.length; i++){
            if(isSet(i)){
                return i;
            }
        }
        // Return encodedSolution.length if no next set bit
        return i;
    }

    /**
     * Finds the ids of the links that make up this solution
     * @param linkIDs
     * @return
     */
    public String idLinks(String[] linkIDs){
        String outputMessage = "Links=";
        for (int i = 0; i < encodedSolution.length; i++){
            if(encodedSolution[i]){
                outputMessage += " " + linkIDs[i];
            }
        }
        return outputMessage;
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;
        Genotype genotype = (Genotype) o;
        return Arrays.equals(encodedSolution, genotype.encodedSolution);
    }

    @Override
    public int hashCode() {
        return Arrays.hashCode(encodedSolution);
    }

    @Override
    public String toString() {
        return "Genotype{" +
                "encodedSolution=" + Arrays.toString(encodedSolution) +
                '}';
    }
}
