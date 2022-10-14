#!/usr/bin/env python3

########################################################################
# Name: classify_reads_sam.py
#
# Description: Take in a SAM file to determine which reads are in each
#              class (truly from bacteria genomes vs. those from the
#              yeast genome).
# 
# Date: June 19th, 2022
#########################################################################

import os
import argparse
import pysam

def build_reference_name_sets(all_refs):
    """ Builds a list of sets to be used to check for header names """
    ref_list = []
    
    # add the sets for each class
    for pos_ref_file in all_refs:
        curr_set = set()
        with open(pos_ref_file, "r") as input_fd:
            all_headers = [curr_set.add(x.strip()) for x in input_fd.readlines()]
        ref_list.append(curr_set)

    return ref_list

def find_class_of_reference_name(ref_list, ref_name):
    """ Finds the class that a particular reference name occurs in """
    class_ids = []
    for i, name_set in enumerate(ref_list):
        if ref_name in name_set:
            class_ids.append(i+1)
    if len(class_ids) != 1:
        print(f"Error: this read hits {len(class_ids)} which is not expected.")
        exit(1)
    return class_ids[0]

def main(args):
    """ main method for script """
    # Open SAM file and create reference name sets
    sam_file = pysam.AlignmentFile(args.sam_file, "r")
    ref_list = build_reference_name_sets(args.class_refs)

    # add all the alignments to dictionary
    read_count = 0
    read_mappings = {}
    print("\n[log] building a dictionary of the read alignments")

    for read in sam_file.fetch():
        class_num = find_class_of_reference_name(ref_list, read.reference_name)
        if read.query_name in read_mappings:
            read_mappings[read.query_name].append(class_num)
        else:
            read_mappings[read.query_name] = [class_num]
        read_count += 1

    sam_file.close()

    # remove alignments that are ambiguous (multiple hits to different species)
    ambiguous_reads = 0
    reads_to_remove = []
    print(f"[log] found {len(read_mappings)} reads and a total of {read_count} alignments")
    print("\n[log] now will start to remove reads that have ambiguous placement")
    
    for key in read_mappings:
        assert len(read_mappings[key]) > 0
        if len(read_mappings[key]) > 1 and len(set(read_mappings[key])) > 1:
            reads_to_remove.append(key)
            ambiguous_reads += 1
    
    for key in reads_to_remove: # remove those reads
        read_mappings.pop(key)
    
    print(f"[log] removed {ambiguous_reads} reads, now we have {len(read_mappings)} reads left.")

    # quantify the number of alignments in each class
    num_classes = len(args.class_refs)

    read_distribution = [0 for i in range(num_classes)]
    for key in read_mappings:
        curr_set = set(read_mappings[key])
        assert len(curr_set) == 1

        read_distribution[read_mappings[key][0]-1] += 1
    normalized_read_dist = [round(x/sum(read_distribution), 3) for x in read_distribution]

    print(f"\n[log] reads found for each class: {read_distribution}")
    print(f"[log] percentages found for each class: {normalized_read_dist}")

    # print out reads for each class separately
    # IMPORTANT: I grab these reads from the original fastq from SRA instead of the 
    # SAM file since that was using the reads with the 720 bp prefix removed.

    raw_read_file = open(args.raw_reads, "r")
    raw_reads = [x.strip() for x in raw_read_file.readlines()]
    assert len(raw_reads) % 4 == 0, "Assertion Error: the number of lines in fastq is not multiple of 4"

    reads_written = set()
    for curr_class_num in range(1, num_classes+1):

        # Open file to write the reads from this current class
        curr_file = open(args.output_dir + f"class_{curr_class_num}_reads.fa", "w")
    
        for i in range(0, len(raw_reads), 4):
            read_name = raw_reads[i].split()[0][1:]
            read_seq = raw_reads[i+1]

            if read_name in read_mappings and read_name not in reads_written:
                # IMPORTANT: this list could have mulitple items, but all they are the same
                class_num = read_mappings[read_name][0]

                if class_num == curr_class_num:
                    curr_file.write(f">{read_name}\n{read_seq[:2000]}\n")
        
        print(f"[log] finished writing the reads from class {curr_class_num} to its own file.")
        curr_file.close()

    raw_read_file.close()
    print("\n[log] finished writing the reads to separate files for each class")

def parse_arguments():
    """ Parse the command-line arguments """
    parser = argparse.ArgumentParser(description="Take in a SAM file and separate out the reads from each class.")
    parser.add_argument("-i", dest="sam_file", help="path to SAM file we will be parsing.", required=True)
    parser.add_argument("-c", dest="class_refs", help="path to txt file with each class reference sequences", required=True, nargs='*')
    parser.add_argument("-o", dest="output_dir", help="path to output directory", required=True)
    parser.add_argument("-r", dest="raw_reads", help="path to original Fastq file with raw reads (not reads with prefix excluded)")
    args = parser.parse_args()
    return args

def check_arguments(args):
    """ Verifies certain aspects about the input parameters """
    if not os.path.isfile(args.sam_file):
        print("Error: the provided SAM file does not exist.")
        exit(-1)
    if not os.path.isdir(args.output_dir):
        print("Error: the provided output directory is not valid.")
        exit(1)
    elif args.output_dir[-1] != "/":
        args.output_dir += "/"
    if not os.path.isfile(args.raw_reads):
        print("Error: the raw read file is not valid.")
        exit(1)
    
if __name__ == "__main__":
    args = parse_arguments()
    check_arguments(args)
    main(args)
