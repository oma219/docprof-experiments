# Experiments related to Document Array Profiles

### Experiment 1: Compare listings to access correctness of data-structure

- This experiment is aimed at making sure that the document array profiles data-structure is computing the correct listings. To do this, I
compare the listings produced by the data-structure versus those using the $r$-index.

### Experiment 2: Compare the query time and index size to the the $`r`$-index :rotating_light: (Used in Paper) :rotating_light:

- This experiment simulates reads, extracts MEMs and computes the listings using the $r$-index and the document array profiles. It compares the
query time and index size for each data-structure.

### Experiment 3: Assess read classification of AMR genes using SPUMONI 2 and document array profiles

- This experiment simulates E. coli reads with AMR genes inserted into them, the goal is classify what sub-class of the gene class does the read contain. 
We compare the results from using SPUMONI 2's sampled document array to using the document array profiles.

### Experiment 4: Assess read classification of real nanopore mock community reads

- This experiment uses real nanopore mock community reads, we split the reads into seperate classes, and attempt to classify using both SPUMONI 2
and the document array profiles.

### Experiment 5: Assess read classification of on any general dataset :rotating_light: (Used in Paper) :rotating_light:

- This experiment is a general experiment that simulates reads from different datasets and performs classification using SPUMONI 2 and the document array profiles.

### Experiment 6: Assess strain-level classification based on real mock community reads :rotating_light: (Used in Paper) :rotating_light:

- This experiment takes the real nanopore reads from a dataset, and performs strain-level classificaton for each bacterial species.

### Experiment 7: Computes sequence homology across different datasets using ANI :rotating_light: (Used in Paper) :rotating_light:

- Takes in a dataset, computes all-pairs ANI values across datasets and summarize in csv file

### Experiment 8: Analyzes the monotonic increases of the DAP in both directions

### Experiment 9: Analyzes the size of the taxonomic DAP with different size columns

### Experiment 10: Verifies the output results from the taxonomic DAP matches the r-index

### Experiment 11: Tests out the taxonomic DAP using a small taxonomy and plots the results

### Experiment 12: Compares the construction of two-pass, heursitic, and no-heuristic of the DAP

### Experiment 13: Builds progressively larger DAPs over subsets of genomes using two-pass construction and profiles it

### Experiment 14: Extract genera from SILVA and build taxonomic compressed DAP (uses rank file from SILVA)
