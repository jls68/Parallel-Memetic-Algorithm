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
     * Get a subset of the genotype
     * @param startIndex inclusive
     * @param endIndex exclusive
     * @return
     */
    public Genotype get(int startIndex, int endIndex) {
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
}
