# K-mer Counting Algorithms

This directory contains implementation of several k-mer counting algorithms in Java. They take `.fa` files as input
and output the count of frequently-appearing k-mers.

The implementation is designed to be using as little memory as possible. For a dataset with about 23,510,000 20-mers,
we can use as little as 200KB of memory (theoretically) to count those 20-mers that appear at least 100 times. Speed
is, however, sacrificed for low memory usage. The memory and speed analysis is provided in the results part.

## Usage

To compile the code, you may use the following

```shell
cd src/main/java
javac kmer_count.java
```

We allow user to specify the value of `k` (length of k-mer) and `q` (the threshold of minimum k-mer frequency to be 
counted). To run the code, you may use the following,

```shell
java -Xmx3m kmer_count [fasta_file_directory] [k] [q] > [output_file_directory]
```

where the `-Xmx3m` flag is used to limit heap space to 3MB (which can be deleted for faster running).

## Method

I tried two methods of k-mer counting.

1. *Bloom filter counter* (`src/main/java/kmer_count.java`): which uses bloom filter to filter out the k-mers that 
   might be appearing frequently. Those k-mers are then inserted into a hash table (`HashMap` object), whose exact frequency are later
   counted.
2. *Misra-Gries* (`src/main/java/kmer_count_misra_gries.java`): We use [Misra-Gries summary](https://en.wikipedia.org/wiki/Misra%E2%80%93Gries_summary), which only use one single
   counter in the whole process. A counter of size `m/q` (where `m` is the total number of k-mers) would be enough to count
   all k-mers whose frequency is above `q`. Has good performance when `m/q` is small.

## Implementation

To use as little memory as possible, we implement the following Java class for convenience.

### `KMer` class

We use Java `BitSet` object to store each k-mer, and try to use only `2 * k` bits to store each k-mer (in reality
it takes `ceiling(2 * k / 64) * 64` bits for each k-mer). For a 20-mer, which usually use 20 bytes of memory, now only needs
40 bits (64 bits or 8 bytes in reality) to store.

This class has a `scan()` function that takes in a DNA string as input, and map `a` to `00`, `c` to `01`, `g` to `10` and `t` to `11`.
It automatically converts a k-mer into its canonical form. That is, if the reverse complement of a k-mer is lexicographically
smaller, than it is changed to its reverse complement.

### `BFCounter` class

We use Java `ArrayList<Short>` to implement the bloom filter counter. User can specify the size `s` of the bloom filter (length
of the array), number of hash functions and a prime `p` for the hash functions. `p` should be much larger than `s`.

Using the idea of *prime field*, the set of hash functions will be automatically generated. For example, if we set the number of hash functions to 2, then
we randomly chose parameters `a1`, `a2`, `b1`, `b2` in `[1, p]`. And the hash functions will be

```
h1(x) = (a1 * x + b1) % p % s
h2(x) = (a2 * x + b2) % p % s
```

If we call `BFCounter.insert(x)`, then the corresponding slots (`h1(x)`, `h2(x)`, etc.) in the arraylist will be incremented.
The smallest value among those slots will be returned as count.


## Results

### Space and Time Complexity

1. For *bloom filter counter*, we will require the expected value in each slot of bloom filter counter to be
   less than `q / 2` so that we can filter out as few k-mers potentially with count above `q` as possible.
   For this reason, the size of bloom filter must be above `m * h / (q / 2)`, where `m` is the total number
   of k-mers and `h` is the number of hash functions. Suppose we use `short` (16 bits) for each slot in bloom filter,
   we will need `4 * m * h / q` bytes in total for the bloom filter.
   
   We can reduce the bloom filter size by dividing the `m` k-mers into `b` bins (in my implementation, I divide them by their
   hash values) so that in each bin the expected number of k-mers is just `m / b`. In this case we only need
   the bloom filter to be of size `4 * m * h / q / b` bytes.

   For the time complexity, we will need to go through the file 2 times for each bin. In total we will need
   `O(2 * b * m)` time.

2. For *Misra-Gries*, we need the counter to be at least of size `m / q` so that we don't miss any k-mer with frequency
   above `q`. Suppose we use 128 bits to store a k-mer and 16 bits to store the count, we will need in total
   `18 * m / q` bytes.

   We can also make this smaller by dividing k-mers into `b` bins and only use `18 * m / q / b` bytes.
   
   For time complexity, it should be the same as previous case, `O(2 * b * m)`.

### Theoretical Memory Usage

I find that the real memory consumption is usually a lot more than what my private fields should take, so I am doing
a theoretical memory usage analysis here.

1. For *bloom filter counter*, I use a fixed size for bloom filter 100,000, and use `short` for the counting. To reduce
   the expected value of count, we just use 2 hash functions. The bloom filter should take 200,000 bytes or
   195KB in total. For the k-mer counter (`HashMap` object), we should have less than 100 k-mers to be counted 
   for each bin, so the counter should take at most `100 * (128 + 16 + 8)` bits, which is neglectable.

   This can be run under the `-Xmx3m` flag and the memory consumption is indeed small.

2. For *Misra-Gries*, we also use a fixed size for `k-mer` counter. In my case I used 20,000, which should only
   take `20000 * (128 + 16 + 8)` bits or 371KB in total.

   This method will fail under the `-Xmx3m` flag, which I think may be due to the `KMer` object is taking up
   space more than 128 bits.

### Wall-clock running time

I'm only testing the bloom filter counter method for now because it's using less memory. For the wall-clock running time I use the command

```bash
time java [-Xmx3m] kmer_count [fasta_file_directory] [k] [q]
```

and record the time as the following.

|Test case|Time|Time with `-Xmx3m` flag|
|---------|----|-----------------------|
|15mer_ken|0m1.030s|0m1.003s|
|50mer_ken|0m0.970s|0m0.919s|
|10mer_chr2L|0m16.992s|*Fail*|
|20mer_chr2L|3m2.516s|8m45.933s|

The output of my program is stored under the `results` folder.

I'm not sure why the third case fail under `-Xmx3m` flag (it should be using even less memory than the 4th case).
This is a bug that I could not solve for now.
