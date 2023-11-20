"""
Purpose: Script to remove PCR duplicates from a 
         FASTQ file. 

         Motivation was that I noticed
         in the 16S rRNA dataset there were many
         duplicate reads that weren't really 
         helping me test the classification power
         so I decided to remove them.
"""

import os
import sys
import argparse


def go_through_reads_and_remove_duplicate_reads(mate1_file, mate2_file, output_dir, sample_name):
    # Generate output file paths
    mate1_file_out = args.output_dir + os.path.basename(mate1_file)
    mate2_file_out = args.output_dir + os.path.basename(mate2_file)

    # Load in all the lines for both files
    with open(mate1_file, "r") as mate1_fd, open(mate2_file, "r") as mate2_fd:
        mate1_lines = [x.strip() for x in  mate1_fd.readlines()]
        mate2_lines = [x.strip() for x in  mate2_fd.readlines()]
    
    # Open the output files, and print out only read sequences that are
    # unique in the mate1 file
    with open(mate1_file_out, "w") as mate1_out, open(mate2_file_out, "w") as mate2_out:
        def print_read_in_fastq(header, seq, qual, out_fd):
            out_fd.write(f"{header}\n{seq}\n+\n{qual}\n")

        header_m1 = ""; seq_m1 = ""; qual_m1 = ""; pos = 0;
        header_m2 = ""; seq_m2 = ""; qual_m2 = ""; pos = 0;
        seq_list = set()

        for mate1_line, mate2_line in zip(mate1_lines, mate2_lines):
            if mate1_line.startswith(f"@{sample_name}"):
                assert mate2_line.startswith(f"@{sample_name}")
                header_m1 = mate1_line; header_m2 = mate2_line 
                pos += 1
            elif pos == 1:
                seq_m1 = mate1_line; seq_m2 = mate2_line
                pos += 1
            elif pos == 2:
                pos += 1
            elif pos == 3:
                qual_m1 = mate1_line; qual_m2 = mate2_line
                pos = 0
                
                # key statement: make sure sequence has not
                # occurred in the read set
                if seq_m1 not in seq_list:
                    seq_list.add(seq_m1)
                    print_read_in_fastq(header_m1, seq_m1, qual_m1, mate1_out)
                    print_read_in_fastq(header_m2, seq_m2, qual_m2, mate2_out)

def check_args(args):
    if not os.path.isfile(args.mate1):
        print("Error: input file 1 provided is not valid.")
        exit(1)
    if not os.path.isfile(args.mate2):
        print("Error: input file 2 provided is not valid.")
        exit(1)
    if not os.path.isdir(args.output_dir):
        print("Error: output directory is not valid.")
        exit(1)
    elif args.output_dir[-1] != "/":
        args.output_dir += "/"

def parse_arguments():
    parser = argparse.ArgumentParser(description="parse read files and output directory...")
    parser.add_argument("--mate1", dest="mate1", help="path to mate1 file", required=True)
    parser.add_argument("--mate2", dest="mate2", help="path to mate2 file", required=True)
    parser.add_argument("--output-dir", dest="output_dir", help="path to output_dir", required=True)
    parser.add_argument("--sample", dest="sample_name", help="name of sample (e.g. A500_V12)", required=True)
    args = parser.parse_args()
    return args

if __name__ == "__main__":
    args = parse_arguments()
    check_args(args)
    go_through_reads_and_remove_duplicate_reads(args.mate1, 
                                                args.mate2, 
                                                args.output_dir, 
                                                args.sample_name)
