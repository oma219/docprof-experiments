########################################################
# Randomly generates reads from input FASTA
# file with a certain read length and certain error rate.
#########################################################

import sys
import random
import argparse

def main(input_file, total_reads_needed, percent_err, sim_rna, read_length):
    # determine the mutation string options
    if sim_rna:
        mut_str = "ACUG"
    else:
        mut_str = "ACTG"

    # concatenate all the sequence together
    all_seqs = ""
    with open(input_file, "r") as input_fd:
        for line in input_fd:
            if '>' not in line:
                all_seqs += line.strip()

    if len(all_seqs) < 2*read_length:
        print("Error: reference sequence is too small to simulate the requested reads.")
        exit(1)

    # simulate all the reads and print them to stdout
    for read_num in range(total_reads_needed):
        index = random.randrange(0, len(all_seqs) - read_length)
        curr_read = all_seqs[index:index+read_length]

        # mutate each character with probability p
        mut_curr_read = ""
        for i, base in enumerate(curr_read):
            if random.random() < percent_err:
                pot_mutations = [ch for ch in mut_str if base != ch]
                new_base = random.choice(pot_mutations)
            else:
                new_base = base
            mut_curr_read += new_base
        
        print(f">read_{read_num}\n{mut_curr_read}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="This script is meant to generate reads.")
    parser.add_argument("-i", "--input", dest="input_file", required=True, type=str)
    parser.add_argument("-e", "--error", dest="percent_err", required=True, type=float)
    parser.add_argument("-n", "--num", dest="num_reads", required=True, type=int)
    parser.add_argument("-l", "--length", dest="read_length", required=True, type=int)
    parser.add_argument("--RNA", dest="sim_rna", action="store_true", default=False)
    args = parser.parse_args()

    input_file = sys.argv[1]
    main(args.input_file, args.num_reads, args.percent_err, args.sim_rna, args.read_length)