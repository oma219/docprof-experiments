########################################################
# Randomly generates Illumina reads from input FASTA
# file with no errors.
#########################################################

import sys
import random
import argparse

def output_reads(seq, num_reads_so_far, total_reads_needed, percent_err, mut_str, num_reads_per_seq=10000):
    read_length = 100
    if len(seq) > 200:
        for i in range(num_reads_per_seq):
            index = random.randrange(0, len(seq) - read_length)
            curr_read = seq[index:index+read_length]

            # Mutate each character with probability p
            mut_curr_read = ""
            for i, base in enumerate(curr_read):
                if random.random() < percent_err:
                    pot_mutations = [ch for ch in mut_str if base != ch]
                    new_base = random.choice(pot_mutations)
                else:
                    new_base = base
                mut_curr_read += new_base
            
            print(f">read_{num_reads_so_far}\n{mut_curr_read}")
            num_reads_so_far += 1

            if num_reads_so_far == total_reads_needed:
                return num_reads_so_far
    return num_reads_so_far

def gen_reads(input_file, total_reads_needed, percent_err, sim_rna):
    # Initialize variables
    num_seqs = 0
    num_reads_so_far = 0

    curr_seq_header = ""
    curr_seq = ""

    # Determine the mutation string options
    if sim_rna:
        mut_str = "ACUG"
    else:
        mut_str = "ACTG"

    # Read through the input file
    with open(input_file, "r") as input_fd:
        for line in input_fd:
            if '>' in line and num_seqs == 0:
                num_seqs += 1
                curr_seq_header = line.strip()
            elif '>' in line and num_seqs > 0:
                num_reads_so_far = output_reads(curr_seq, num_reads_so_far, total_reads_needed, percent_err, mut_str)
                num_seqs += 1
                curr_seq_header = line.strip()
                curr_seq = ""
            else:
                curr_seq += line.strip()
            
            # Check if we have found all the reads
            if num_reads_so_far == total_reads_needed:
                break
        
        if num_reads_so_far == total_reads_needed:
            return
        num_reads_so_far = output_reads(curr_seq, num_reads_so_far, total_reads_needed, percent_err, mut_str)
            
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="This script is meant to generate reads.")
    parser.add_argument("-i", "--input", dest="input_file", required=True, type=str)
    parser.add_argument("-e", "--error", dest="percent_err", required=True, type=float)
    parser.add_argument("-n", "--num", dest="num_reads", required=True, type=int)
    parser.add_argument("--RNA", dest="sim_rna", action="store_true", default=False)
    args = parser.parse_args()

    input_file = sys.argv[1]
    gen_reads(args.input_file, args.num_reads, args.percent_err, args.sim_rna)