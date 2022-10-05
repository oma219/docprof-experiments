#!/usr/bin/env python3

# Name: separate_gene_classes.py
# Description: This is a python script is used by experiment 3 in order
#              to take in a list of genes from different sub-classes of
#              AMR genes and then separate them into separate FASTA files.
#
# Date: October 5, 2022

import os
import argparse
import sys

def main(args):
    """ main method: generates those separate FASTA files for each class """

    # Declare some variables for parsing input FASTA
    header_lines = [] # Each is formatted as a list [">seq_0_class_A", "A"] with class name in second position
    seq_lines = []
    
    class_set = set()
    num_items_per_class = {}
    length_per_class = {}

    # IMPORTANT: this loop below is assuming each sequence is located into 1 line, it 
    #            is not spread over multiple lines, that is how I found the file to be formatted.

    # Open the full gene file and verify the number of classes
    with open(args.input_file, "r") as input_fd:
        go_for_extract = False # Some genes are not labeled
        class_name = ""

        for line in input_fd:
            if ">" in line:
                go_for_extract = True
                header_info = ["", ""]
                header_info[0] = line.strip()

                class_name = line.split("|")[3].split("_")[1]
                header_info[1] = line.split("|")[3].split("_")[1]

                # Makes sure this is a labeled gene
                if len(class_name) == 1:
                    class_set.add(class_name)
                    header_lines.append(header_info)

                    # Keep track of the number of items in each class
                    if class_name not in num_items_per_class:
                        num_items_per_class[class_name] = 1
                    else:
                        num_items_per_class[class_name] += 1
                else:
                    go_for_extract = False
            else:
                if go_for_extract:
                    seq_lines.append(line.strip())
                    if class_name not in length_per_class:
                        length_per_class[class_name] = len(line.strip())
                    else:
                        length_per_class[class_name] += len(line.strip())

    assert len(header_lines) == len(seq_lines), "issue in parsing the complete gene set"
    assert len(class_set) == args.num_datasets, "number of datasets in file does not match parameters"

    print(f"[log] num_items_per_class = {num_items_per_class}")
    print(f"[log] length_per_class = {length_per_class}")

    # Determine the smallest class by number of bases (we want to balance the class sizes)
    min_class_bp = min([length_per_class[x] for x in length_per_class])
    output_length_per_class = {}

    print(f"[log] smallest class size in bp = {min_class_bp}")

    class_num = 1
    for class_name in sorted(class_set):
        total_length = 0
        pos = 0
        
        with open(args.output_dir + f"class_{class_num}.fna", "w") as out_fd:
            while pos < len(header_lines) and total_length < min_class_bp:
                if header_lines[pos][1] == class_name:
                    total_length += len(seq_lines[pos])
                    out_fd.write(f"{header_lines[pos][0]}\n{seq_lines[pos]}\n")
                pos += 1
            
            output_length_per_class[class_name] = total_length
            class_num += 1
    
    print(f"[log] output_length_per_class = {output_length_per_class}")


def parse_arguments():
    """ Defines the command-line argument parser, and return arguments """
    parser = argparse.ArgumentParser(description="This script helps to separate out the genes from each class into a separate file.")
    parser.add_argument("-n", "--num", dest="num_datasets", required=True, help="number of different sub-classes of genes", type=int)
    parser.add_argument("-i", dest="input_file", required=True, help="input file with all the genes")
    parser.add_argument("-o", dest="output_dir", required=True, help="output directory for all the FASTA file")
    args = parser.parse_args()
    return args

def check_arguments(args):
    """ Checks for invalid arguments """
    if not os.path.isdir(args.output_dir):
        print("Error: output directory does not exist.")
        exit(1) 
    if args.num_datasets < 1:
        print("Error: number of datasets must be greater than 0")
        exit(1)
    if not os.path.isfile(args.input_file):
        print("Error: path to input file is not valid")
        exit(1)

if __name__ == "__main__":
    args = parse_arguments()
    check_arguments(args)
    main(args)