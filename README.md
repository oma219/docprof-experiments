# Experiments related to Document Array Profiles


### Experiment 1: Compare listings to access correctness of data-structure

- This experiment is aimed at making sure that the document array profiles data-structure is computing the correct listings. To do this, I
compare the listings produced by the data-structure versus those using the $r$-index.

### Experiment 2: Compare the query time and index size to the the $r$-index

- This experiment simulates reads, extracts MEMs and computes the listings using the $r$-index and the document array profiles. It compares the
query time and index size for each data-structure.

### Experiment 3: Assess read classification of AMR genes using SPUMONI 2 and document array profiles

- This experiment simulates E. coli reads with AMR genes inserted into them, the goal is classify what sub-class of the gene class does the read contain. 
We compare the results from using SPUMONI 2's sampled document array to using the document array profiles.

### Experiment 4: Assess read classification of real nanopore mock community reads

- This experiment uses real nanopore mock community reads, we split the reads into seperate classes, and attempt to classify using both SPUMONI 2
and the document array profiles.
