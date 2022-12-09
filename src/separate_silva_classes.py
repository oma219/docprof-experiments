#!/usr/bin/env python3

# Name: separate_silva_classes.py
# Description: This is a python script that will separate the 
#              SILVA database based on a given level of the
#              taxonomy.
#
# Date: December 4th, 2022

import os
import sys
import argparse

def main(args):
    """ main method of the script """
    
    # Step 1: Extract the lines from the SILVA database
    all_lines = []
    with open(args.input_file, "r") as in_fd:
        all_lines = [x.strip() for x in in_fd.readlines()]

    # Step 2: Focus only on sequences that have an "understandable" taxonomy
    filtered_lines = []
    tax_id_set = {}
    counts = [0 for i in range(20)]
    total_num_seqs = 0

    for i, line in enumerate(all_lines):
        if '>' in line:
            taxonomy = line.split()[1].split(";")
            total_num_seqs += 1
            if taxonomy[-1] == taxonomy[-2] and "uncultured" != taxonomy[-1]:
                counts[len(taxonomy)] += 1
                filtered_lines.append(line); filtered_lines.append(all_lines[i+1])
                #tax_id_set.add(taxonomy[-1])
                if taxonomy[-1] not in tax_id_set:
                    tax_id_set[taxonomy[-1]] = 1
                else:
                    tax_id_set[taxonomy[-1]] += 1

    print(f"[log] input file had {total_num_seqs} seqeunces in it.")
    print(f"[log] we filtered it down to {sum(counts)} sequences\n")
    print(f"[log] number of taxonomy levels: {counts}")
    print(f"[log] number of genera found: {len(tax_id_set)}\n")

    # Sort based on number of sequences in each genus
    sorted_tax_id_set = dict(sorted(tax_id_set.items(), key=lambda item: item[1], reverse=True))
    
    # Step 3: Write out the sequences in each group separately
    print(f"[log] writing out the separate FASTA files for top 1000 genera")
    for i, key in enumerate(sorted_tax_id_set):
        # Only focus on top 1000 classes
        if i >= 100:
            break
        # Write out all the sequences for this genus
        with open(args.output_dir + f"tax_group_{i+1}.fna", "w") as out_fd:
            for j, line in enumerate(filtered_lines):
                if '>' in line:
                    taxonomy = line.split()[1].split(";")
                    if taxonomy[-2] == key:
                        out_fd.write(f"{line}\n{filtered_lines[j+1]}\n")
    
    # Step 4: Write out the filelist
    with open(args.output_dir + "filelist.txt", "w") as out_fd:
        for i in range(100):
            out_fd.write(f"{args.output_dir}tax_group_{i+1}.fna {i+1}\n")



def parse_arguments():
    """ Defines the command-line argument parser, and return arguments """
    parser = argparse.ArgumentParser(description="This script allows users to separate a SILVA database based on a "
                                                  "specific level of the taxonomy. ")
    parser.add_argument("-i", "--input", dest="input_file", required=True, help="input SILVA database")
    parser.add_argument("-o", "--output-dir", dest="output_dir", required=True, help="output directory for files")
    args = parser.parse_args()
    return args

def check_arguments(args):
    """ Validate the command-line arguments """
    if not os.path.isfile(args.input_file):
        print("Error: the path for the provided input file is not valid.")
        exit(1)
    if not os.path.isdir(args.output_dir):
        print("error: the output directory path is not valid.")
        exit(1)
    elif args.output_dir[-1] != '/':
        args.output_dir += "/"

if __name__ == "__main__":
    args = parse_arguments()
    check_arguments(args)
    main(args)