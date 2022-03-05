import java.io.*;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.*;

import static java.lang.System.gc;

class InvalidKMerException extends Exception {
    /**
     * An exception used to report that the k-mer is invalid.
     */
    public InvalidKMerException(String str) {
        super(str);
    }
}

class kmer_count {
    /**
     * The KMer class storing a k-mer in 2k bits (in reality ceiling(2k/64)*64 bits)
     */
    private static class KMer implements Comparable<KMer> {
        /**
         * Bits storing the information in the form of bits
         * A: 00, C: 01, G: 10, T: 11
         */
        public BitSet bits;

        /**
         * Constructor of the KMer object
         * @param k length of k-mer
         */
        public KMer(int k) {
            bits = new BitSet(2 * k);
        }

        /**
         * Set the content of this to be equal to another k-mer
         * @param that the other k-mer
         */
        public void set(KMer that) {
            bits = that.bits;
        }

        /**
         * Take in a DNA string and convert it into bit array. Set itself to its reverse
         * complement if the reverse complement is lexicographically smaller (canonical
         * form).
         * @param dna the DNA string
         * @throws InvalidKMerException indicating that the length of the string is incorrect
         * or the string contains illegal character (other than a/A, c/C, g/G, t/T)
         */
        public void scan(String dna) throws InvalidKMerException {
            bits = new BitSet(2 * k);
            // check string length
            if (dna.length() != k) {
                bits = null;
                throw new InvalidKMerException("Incorrect k-mer length");
            }
            // set bit1 and bit2 according to the dna
            for (int i = 0; i < k; i++) {
                switch (dna.charAt(i)) {
                    case 'A':
                    case 'a':
                        // 0, 0
                        break;
                    case 'C':
                    case 'c':
                        // 0, 1
                        bits.set(2*k-2-2*i);
                        break;
                    case 'G':
                    case 'g':
                        // 1, 0
                        bits.set(2*k-1-2*i);
                        break;
                    case 'T':
                    case 't':
                        // 1, 1
                        bits.set(2*k-2-2*i);
                        bits.set(2*k-1-2*i);
                        break;
                    default:
                        // e.g. encountering 'N', ignore the k-mer
                        bits = null;
                        throw new InvalidKMerException("KMer contains illegal character " + dna.charAt(i));
                }
            }
            // Check if reverse complement is smaller
            KMer comp = reverseComplement();
            if (!lessThan(comp)) {
                this.set(comp);
            }
        }

        /**
         * Compute the reverse complement of this.
         * @return the KMer object that shows the reverse complement.
         */
        public KMer reverseComplement() {
            KMer comp = new KMer(k);
            for (int i = 0; i < k; i++) {
                if (!bits.get(2*k-2-2*i)) {
                    comp.bits.set(2*i);
                }
                if (!bits.get(2*k-1-2*i)) {
                    comp.bits.set(2*i+1);
                }
            }
            return comp;
        }

        /**
         * Check if this is equal to another k-mer. Used in Hashtable object.
         * @param that the other k-mer
         * @return the boolean value indicating whether this == that
         */
        @Override
        public boolean equals(Object that) {
            KMer t = (KMer) that;
            return bits.equals(t.bits);
        }

        /**
         * Convert this to the DNA string. All characters are in lower case alphabet.
         * @return the original DNA string.
         */
        public String toString() {
            if (bits == null) {
                return "";
            }
            StringBuilder dna = new StringBuilder();
            for (int i = k-1; i >= 0; i--) {
                if (!bits.get(2*i+1) && !bits.get(2*i)) {
                    dna.append('a');
                } else if (!bits.get(2*i+1) && bits.get(2*i)) {
                    dna.append('c');
                } else if (bits.get(2*i+1) && !bits.get(2*i)) {
                    dna.append('g');
                } else dna.append('t');
            }
            return dna.toString();
        }

        /**
         * @return The hash code of the k-mer. Used in HashMap/Hashtable objects.
         */
        @Override
        public int hashCode() {
            return bits.hashCode();
        }

        /**
         * @param that Another k-mer
         * @return Check if this is lexicographically smaller than that.
         */
        public boolean lessThan(KMer that) {
            for (int i = 2*k-1; i >= 0; i--) {
                if (bits.get(i) ^ that.bits.get(i)) {
                    return that.bits.get(i);
                }
            }
            return false;
        }

        /**
         * The compareTo method of a comparable object. Can be used to
         * sort a bunch of KMer objects.
         * @param that the other k-mer
         * @return an integer that is positive if this > that, 0 if this == that
         * and negative otherwise.
         */
        @Override
        public int compareTo(KMer that) {
            if (lessThan(that)) return -1;
            else if (equals(that)) return 0;
            else return 1;
        }
    }


    private static class BFCounter {
        // bloom filter counter
        private static ArrayList<Short> counter;

        // parameters for hash functions for bloom filter
        private static ArrayList<Integer> hashParam1;
        private static ArrayList<Integer> hashParam2;
        private static int hashPrime;

        /**
         * Constructor of the bloom filter counter.
         * @param counterSize the size of bloom filter.
         * @param hashSize number of hash functions.
         * @param prime the prime number for hash functions. Should be a prime
         *              that is much larger than counterSize.
         */
        public BFCounter(int counterSize, int hashSize, int prime) {
            counter = new ArrayList<>();
            for (int i = 0; i < counterSize; i++) {
                counter.add((short) 0);
            }
            hashParam1 = new ArrayList<>();
            hashParam2 = new ArrayList<>();
            hashPrime = prime;

            Random rand = new Random();
            for (int i = 0; i < hashSize; i++) {
                hashParam1.add(rand.nextInt(prime));
                hashParam2.add(rand.nextInt(prime));
            }
        }

        /**
         * Insert the kmer into the bloom filter counter and return the estimated count
         * of the k-mer.
         * @param kmer the KMer object to be counted.
         * @return the estimated count of the k-mer.
         */
        public int insert(KMer kmer) {
            // insert kmer into the BFCounter and return the count of the kmer
            short min = 0x7FFF;
            int slot;
            for (int i = 0; i < hashParam1.size(); i++) {
                slot = (hashParam1.get(i) * kmer.hashCode() + hashParam2.get(i)) % hashPrime % counter.size();
                if (slot < 0) {
                    slot += counter.size();
                }
                counter.set(slot, (short) (counter.get(slot) + 1));
                if (counter.get(slot) < min) {
                    min = counter.get(slot);
                }
            }
            return min;
        }

        /**
         * clear all counts in the bloom filter.
         */
        public void reset(int counterSize) {
            for (int i = 0; i < counterSize; i++) {
                counter.set(i, (short) 0);
            }
        }
    }

    /**
     * Static parameters
     */

    // bloom filter counter for filtering out the k-mers with high counts
    private static BFCounter bfCounter;

    // a hash table storing the counts of k-mers
    private static HashMap<KMer, Short> counter;

    // the length of k-mers
    private static int k;

    // the threshold of k-mer frequency that is to be counted
    private static int q;

    // total number of k-mers to be counted
    private static int m;

    // number of splits of memory bin
    private static int numBins;
    private static int currBin;

    interface KMerOperation {
        void operation(KMer kmer);
    }

    /**
     * Count the total number of k-mer in the file
     */
    private static final KMerOperation count = (KMer kmer) -> m += 1;

    /**
     * Insert the k-mer into BFCounter. If the returned estimated count
     * is larger or equal to q, then store the k-mer in counter.
     * We only process those k-mers whose hashCode % numBins == currBin
     * @param kmer The KMer object to be counted.
     */
    private static final KMerOperation insertCounter = (KMer kmer) -> {
        int bin = kmer.hashCode() % numBins;
        if (bin < 0) bin += numBins;
        if (bin == currBin) {
            // insert kmer into the counter and return the count of the kmer
            int count = bfCounter.insert(kmer);
            if (count >= q) counter.put(kmer, (short) 0);
        }
    };

    /**
     * Count k-mers that are already present in the counter.
     */
    private static final KMerOperation countCounter = (KMer kmer) -> {
        if (counter.containsKey(kmer)) {
            counter.put(kmer, (short) (counter.get(kmer) + 1));
        }
    };

    /**
     * Read the fasta file and process all the k-mers in the file.
     * @param fileName the name of fasta file
     * @param op the operation to be done on each k-mer.
     */
    private static void readFile(String fileName, KMerOperation op) {
        StringBuilder kmerString;
        KMer kmer;
        try (BufferedReader reader = new BufferedReader(new FileReader(fileName))) {
            char next;
            while (reader.readLine() != null) {
                kmerString = new StringBuilder(k);
                System.gc();
                while (reader.ready()) {
                    next = (char) reader.read();
                    if (next == '\n') break;
                    kmerString.append(next);
                    // keep kmerString length == k
                    if (kmerString.length() < k) continue;
                    while (kmerString.length() > k) kmerString.deleteCharAt(0);

                    try {
                        // count the k-mer and put the value in counter
                        kmer = new KMer(k);
                        kmer.scan(kmerString.toString());
                        op.operation(kmer);
                    } catch (InvalidKMerException ignored) {}
                }
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    /**
     * Print out the result of counting.
     */
    private static void print() {
        for (Map.Entry<KMer, Short> entry : counter.entrySet()) {
            if (entry.getValue() >= q) {
                System.out.println(entry.getValue() + " " + entry.getKey());
            }

        }
    }

    private static void clear() {
        for (Map.Entry<KMer, Short> entry : counter.entrySet()) {
            entry.setValue((short) 0);
        }
    }

    private static void reset() {
        counter = new HashMap<>();
        gc();
    }

    public static void main(String[] args) {
        // prime number used for the hash functions in bloom filter
        // preferred to be much larger than BLOOM_FILTER_SIZE
        int BLOOM_FILTER_PRIME = 1982627;

        // Number of hash functions used in bloom filter
        int BLOOM_FILTER_HASH = 2;
        // Size of the bloom filter (length of the array)
        int BLOOM_FILTER_SIZE = 100000;


        if (args.length != 3) {
            System.out.println("Usage: java kmer_count [fasta_file_directory] [k] [q]");
        }
        String file = args[0];
        k = Integer.parseInt(args[1]);
        q = Integer.parseInt(args[2]);


        bfCounter = new BFCounter(BLOOM_FILTER_SIZE, BLOOM_FILTER_HASH, BLOOM_FILTER_PRIME);
        try {
            // use the size of the file to estimate the number of k-mers to be parsed
            int bytes = (int) Files.size(Paths.get(file));

            // set the number of bins
            numBins = bytes * BLOOM_FILTER_HASH / (q/2) / BLOOM_FILTER_SIZE + 1;
            //System.out.println(bytes + " " + numBins);
            for (currBin = 0; currBin < numBins; currBin++) {
                reset();                             // reset the counter
                readFile(file, insertCounter);       // use bloom filter to filter out frequent k-mers
                readFile(file, countCounter);        // count the exact appearance of each k-mer
                print();                             // print the result of counting
                bfCounter.reset(BLOOM_FILTER_SIZE);  // reset the bloom filter
            }
        } catch (IOException e) {
            e.printStackTrace();
        }

    }




}
