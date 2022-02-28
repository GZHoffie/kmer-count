import java.io.*;
import java.util.*;

class InvalidKMerException extends Exception {
    public InvalidKMerException(String str) {
        super(str);
    }
}

class kmer_count {

    private static class KMer implements Comparable<KMer> {

        // store the information in the form of bits
        // A: 00, C: 01, G: 10, T: 11
        public BitSet bits;

        public KMer(int k) {
            bits = new BitSet(2 * k);
        }

        public void set(KMer that) {
            bits = that.bits;
        }

        // REQUIRES: dna is a string of length k
        public void scan(String dna, int k) throws InvalidKMerException {
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

        @Override
        public boolean equals(Object that) {
            KMer t = (KMer) that;
            return bits.equals(t.bits);
        }

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

        @Override
        public int hashCode() {
            return bits.hashCode();
        }

        public boolean lessThan(KMer that) {
            for (int i = 2*k-1; i >= 0; i--) {
                if (bits.get(i) ^ that.bits.get(i)) {
                    return that.bits.get(i);
                }
            }
            return false;
        }

        @Override
        public int compareTo(KMer that) {
            if (lessThan(that)) return -1;
            else if (equals(that)) return 0;
            else return 1;
        }
    }

    private static ArrayList<Short> BFCounter;
    private static ArrayList<Integer> hashParam1;
    private static ArrayList<Integer> hashParam2;
    private static int hashPrime;
    private static Hashtable<KMer, Short> counter;
    private static int k;
    private static int q;

    public static void init(int counterSize, int hashSize, int prime) {
        BFCounter = new ArrayList<>();
        for (int i = 0; i < counterSize; i++) {
            BFCounter.add((short) 0);
        }
        counter = new Hashtable<>();
        hashParam1 = new ArrayList<>();
        hashParam2 = new ArrayList<>();
        hashPrime = prime;

        Random rand = new Random();
        for (int i = 0; i < hashSize; i++) {
            hashParam1.add(rand.nextInt(prime));
            hashParam2.add(rand.nextInt(prime));
        }
    }

    private static int insert(KMer kmer) {
        // insert kmer into the BFCounter and return the count of the kmer
        short min = 0x7FFF;
        int slot;
        for (int i = 0; i < hashParam1.size(); i++) {
            slot = (hashParam1.get(i) * kmer.hashCode() + hashParam2.get(i)) % hashPrime % BFCounter.size();
            if (slot < 0) {
                slot += BFCounter.size();
            }
            BFCounter.set(slot, (short) (BFCounter.get(slot) + 1));
            if (BFCounter.get(slot) < min) {
                min = BFCounter.get(slot);
            }
        }
        return min;
    }

    private static void readFile(String fileName) {
        StringBuilder kmerString = new StringBuilder(k);
        KMer kmer;
        int count;
        try (BufferedReader reader = new BufferedReader(new FileReader(fileName))) {
            char next;
            while (reader.readLine() != null) {
                kmerString.setLength(0);
                kmer = new KMer(k);
                while (reader.ready()) {
                    next = (char) reader.read();
                    if (next == '\n') break;
                    kmerString.append(next);
                    // keep kmerString length == k
                    if (kmerString.length() < k) continue;
                    while (kmerString.length() > k) kmerString.deleteCharAt(0);

                    try {
                        kmer.scan(kmerString.toString(), k);
                        count = insert(kmer);
                        if (count >= q) {
                            KMer temp = new KMer(k);
                            temp.set(kmer);
                            counter.put(temp, (short) 0);
                        }
                    } catch (InvalidKMerException ignored) {}
                }
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    private static void count(String fileName) {
        StringBuilder kmerString = new StringBuilder(k);
        KMer kmer;
        try (BufferedReader reader = new BufferedReader(new FileReader(fileName))) {
            char next;
            while (reader.readLine() != null) {
                kmerString.setLength(0);
                kmer = new KMer(k);
                while (reader.ready()) {
                    next = (char) reader.read();
                    if (next == '\n') break;
                    kmerString.append(next);
                    // keep kmerString length == k
                    if (kmerString.length() < k) continue;
                    while (kmerString.length() > k) kmerString.deleteCharAt(0);

                    try {
                        kmer.scan(kmerString.toString(), k);
                        if (counter.containsKey(kmer)) {
                            counter.put(kmer, (short) (counter.get(kmer) + 1));
                        }
                    } catch (InvalidKMerException ignored) {}
                }
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    private static void print() {
        for (Map.Entry<KMer, Short> entry : counter.entrySet()) {
            //if (entry.getValue() >= q) {
            System.out.println(entry.getValue() + " " + entry.getKey());
            //}
        }
    }

    public static void main(String[] args) {
        k = 10;
        q = 4000;
        init(30000, 2, 1982627);
        //KMer k = new KMer(4);
        String file = "../../../ass2_dataset/chr2L_dm6.fa";
        readFile(file);
        count(file);
        print();

        //System.out.println(k.toString());
    }




}
