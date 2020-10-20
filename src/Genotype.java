import java.util.Arrays;
import java.util.List;

public class Genotype {

    byte[] encodedSolution;

    public Genotype(int length){
        // Calculate how many bytes will be needed to encode the state of all possible connections
        int byteCount = (length + 7) / 8;
        encodedSolution = new byte[byteCount];
    }

    public int length(){
        return encodedSolution.length * 8;
    }

    private int getMask(int index){
        return 1 << (index % 8);
    }

    public boolean isSet(int index){
        byte b = encodedSolution[index / 8];
        return (b & getMask(index)) > 0;
    }

    public void flip(int index){
        encodedSolution[index / 8] = (byte) (encodedSolution[index / 8] ^ getMask(index));
    }

    public byte getByte(int byteIndex){
        return encodedSolution[byteIndex];
    }

    /**
     * Get a subset of the genotype
     * @param startIndex inclusive
     * @param endIndex exclusive
     * @return
     */
    public Genotype getSubset(int startIndex, int endIndex) {
        Genotype subset = new Genotype(endIndex - startIndex);
        for(int j = this.nextSetBit(startIndex); j < endIndex; j = this.nextSetBit(j + 1)){
            subset.set(j - startIndex, true);
        }
        return subset;
    }

    /**
     * Set the a bit to a specific value
     * @param index of the bit
     * @param thereIsLink when true sets the bit to one else it is set to zero
     */
    public void set(int index, boolean thereIsLink){
        if(thereIsLink) {
            // Set bit to one
            encodedSolution[index / 8] |= getMask(index);
        }
        else {
            // Clear bit so it is zero
            encodedSolution[index / 8] &= ~getMask(index);
        }
    }

    public void set(int byteIndex, byte newByte){
        encodedSolution[byteIndex] = newByte;
    }

    public Genotype clone(){
        Genotype clone = new Genotype(this.length());
        for(int i = 0; i < encodedSolution.length; i++){
            clone.set(i, encodedSolution[i]);
        }
        return clone;
    }

    public int nextSetBit(int i) {
        for(i = i; i < this.length(); i++){
            if( i < 0){
                return this.length();
            }
            if(isSet(i)){
                return i;
            }
        }
        // Return this.length() if no next set bit
        return i;
    }

    /**
     * Finds the ids of the links that make up this solution
     * @param linkIDs
     * @return
     */
    public String idLinks(String[] linkIDs){
        String outputMessage = "Links=";
        for (int i = 0; i < this.length(); i++){
            if(this.isSet(i)){
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
        String geno = "";
        for (int b = 0; b < encodedSolution.length; b++) {
            String newByte = String.format("%8s", Integer.toBinaryString(encodedSolution[b] & 0xFF)).replace(' ', '0');
            geno += newByte + " ";
        }
        return "Genotype{" +
                "encodedSolution= " + geno +
                '}';
    }
}
