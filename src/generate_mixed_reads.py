#!/usr/bin/env python3

# Name: separate_gene_classes.py
# Description: This is a python script is used by experiment 3 in order
#              to take generate nanopore reads that are mix of bacterial
#              sequence and AMR genes. The goal is for each read to be 
#              a 50:50 mixture of the two.
#
# Date: October 5, 2022

import os
import argparse
import sys
import random

def main(args):
    """ main method: take in bacterial reads, and randomly insert AMR genes into them """
    
    for class_num in range(1, args.num_datasets+1):
        # Load all the gene reads into a list
        gene_reads = []
        with open(args.input_gene_dir + f"class_{class_num}_reads.fna", "r") as input_fd:
            read = ["", ""]
            for line in input_fd:
                if ">" in line:
                    read[0] = line.strip()
                else:
                    read[1] = line.strip()
                    gene_reads.append(read)
                    read = ["", ""]

        print(f"[log] For class {class_num}, there were {len(gene_reads)} gene reads found.")

        # Iterate through the bacteria reads and output chimeras of the two
        start_pos = (class_num - 1) * args.num_reads_per_class
        end_pos = class_num * args.num_reads_per_class
        
        out_fd = open(args.output_dir + f"class_{class_num}_reads.fna", "w")
        output_ratio = [0, 0]

        with open(args.input_reads, "r") as input_fd:
            header = ""
            pos = 0
            num_iterations = 0
            for line in input_fd:
                num_iterations += 1
                if ">" in line:
                    header = line.strip()
                else:
                    seq = line.strip()
                    if pos >= start_pos and pos < end_pos:
                        random_gene = random.choice(gene_reads)[1]
                        out_fd.write(f"{header}\n{seq}{random_gene}\n")
                        output_ratio[0] += len(seq)
                        output_ratio[1] += len(random_gene)
                    elif pos >= end_pos:
                        break
                    pos += 1
        out_fd.close()
        print(f"[log] For class {class_num}, bacterial reads from pos={start_pos} to pos={end_pos} were used.")
        print(f"[log] For class {class_num}, {output_ratio[0]/sum(output_ratio): .3f}% is bacteria, {output_ratio[1]/sum(output_ratio): .3f}% is gene\n")

def parse_arguments():
    """ Defines the command-line argument parser, and return arguments """
    parser = argparse.ArgumentParser(description="This script helps to generate mixed reads between bacteria and AMR genes.")
    parser.add_argument("-n", "--num", dest="num_datasets", required=True, help="number of different sub-classes of genes", type=int)
    parser.add_argument("--num-reads", dest="num_reads_per_class", type=int, help="number of reads for each class")
    parser.add_argument("-r", dest="input_reads", required=True, help="path to input bacteria reads")
    parser.add_argument("-i", dest="input_gene_dir", required=True, help="path to directory with the different classes of genes")
    parser.add_argument("-o", dest="output_dir", required=True, help="output directory for the different output reads")
    args = parser.parse_args()
    return args

def check_arguments(args):
    """ Checks for invalid arguments """
    if not os.path.isdir(args.output_dir):
        print("Error: output directory does not exist.")
        exit(1) 
    if not os.path.isdir(args.input_gene_dir):
        print("Error: input gene directory does not exist.")
        exit(1) 
    if args.num_datasets < 1:
        print("Error: number of datasets must be greater than 0")
        exit(1)
    if not os.path.isfile(args.input_reads):
        print("Error: path to input file is not valid")
        exit(1)

if __name__ == "__main__":
    args = parse_arguments()
    check_arguments(args)
    main(args)