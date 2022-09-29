#!/usr/bin/env python3

# Name: analyze_sam.py
# Description: This is a python script is used by experiment 1 to 
#              analyze the SAM files from r-index and the output
#              when using the document array profiles.
#
# Date: September 12th, 2022

import os
import argparse
import math
import csv
import pysam

def calculate_accuracy_values(confusion_matrix, num_datasets):
    """ Calculates accuracy score given confusion matrix """
    accuracies = []
    for pivot in range(num_datasets):
        tp = confusion_matrix[pivot][pivot]
        fp = fn = tn = 0
        for row in range(num_datasets):
            for column in range(num_datasets):
                curr = confusion_matrix[row][column]
                if(column == pivot and row != pivot):
                    fp += curr
                elif(row == pivot and column != pivot):
                    fn += curr
                elif(row != pivot):
                    tn += curr
        accuracies.append([pivot,tp,tn,fp,fn])
    return accuracies

def main(args):
   
    # Go through alignments for the pivot with respect to each
    # database, and create a mappings dictionary
    read_mappings = {}
    input_entries = 0
    print("\n[log] building a dictionary of the read alignments for pivot {pivot}".format(pivot=args.curr_pivot_num))

    # Iterates through each SAM file (one for each database)
    for j in range(args.num_datasets):
        curr_sam = pysam.AlignmentFile(args.sam_dir+"pivot_{pivot}_align_dataset_{dataset}.sam".format(pivot=args.curr_pivot_num, dataset=j+1), "r")

        for read in curr_sam.fetch():
            query_length = int(read.query_name.split("_")[3])
            if read.query_name not in read_mappings:
                read_mappings[read.query_name] = [query_length]
            read_mappings[read.query_name].append(j)
        curr_sam.close()
        print("[log] For pivot {pivot_num}, finished processing alignments with respect to database {num}".format(pivot_num=args.curr_pivot_num, num=j+1))
    
    # Process the document listing file into a dictionary
    doclist_dict = {}
    with open(args.doclist_file, "r") as input_fd:
        curr_name = ""
        for line in input_fd:
            if ">" in line:
                curr_name = line.strip()[1:]
            else:
                line_split = line.split()
                assert len(line_split) == 2, "Issue with doclisting result"

                curr_doclist = [int(x) for x in line_split[1][1:-2].split(",")]
                doclist_dict[curr_name] = curr_doclist
    
    # Verify for each read the document listings are equalivant
    matching_lists = 0
    unmatched_lists = 0
    for key in read_mappings:
        if read_mappings[key][1:] != doclist_dict[key]:
            unmatched_lists += 1
        else:
            matching_lists += 1
    print("[log] finished comparing to document array profiles, writing to output file")

    # Write results to output file
    out_fd = open(args.output_dir + "output.txt", "w")
    out_fd.write(f"[log] Out of {len(read_mappings)} reads, there were {matching_lists} reads with same document\n"
                 f"      listing from the r-index and document array profiles. While, there were {unmatched_lists}\n"
                 f"      reads with different listings.\n")
    out_fd.close()

def parse_arguments():
    """ Defines the command-line argument parser, and return arguments """
    parser = argparse.ArgumentParser(description="This script helps to analyze the SAM file from experiment 5"
                                                 "in order to form a confusion matrix.")
    parser.add_argument("-n", "--num", dest="num_datasets", required=True, help="number of datasets in this experiment", type=int)
    parser.add_argument("-c", "--current-pivot", dest="curr_pivot_num", help="current pivot that we are aligning", type=int, required=True)
    parser.add_argument("-s", "--sam_file", dest="sam_dir", required=True, help="path to directory with SAM files to be analyzed")
    parser.add_argument("-o", "--output_path", dest = "output_dir", required=True, help="path to output directory")
    parser.add_argument("--doclist", dest="doclist_file", required=True, help="path to document listing output file")
    #parser.add_argument("--half_mems", action="store_true", default=False, dest="half_mems", help="sam corresponds to half-mems")
    #parser.add_argument("--mems", action="store_true", default=False, dest="mems", help="sam corresponds to mems (it can either be mems or half-mems, not both)")
    #parser.add_argument("-t", "--threshold", dest = "t", required=False, default=0, help="optional threshold value for experiment 8", type = int)
    #parser.add_argument("-l", "--length-text", dest = "fasta_index_file", required=False, default="", help="path to .fai file with total reference length")
    args = parser.parse_args()
    return args

def check_arguments(args):
    """ Checks for invalid arguments """
    if not os.path.isdir(args.sam_dir):
        print("Error: sam directory does not exist.")
        exit(1) 
    if args.num_datasets < 1:
        print("Error: number of datasets must be greater than 0")
        exit(1)
    if not os.path.isfile(args.doclist_file):
        print("Error: path to document listing file is not valid")
        exit(1)


    #if (args.half_mems and args.mems) or (not args.half_mems and not args.mems):
    #    print("Error: exactly one type needs to be chosen (--half-mems or --mems)")
    #    exit(1)

if __name__ == "__main__":
    args = parse_arguments()
    check_arguments(args)
    main(args)

