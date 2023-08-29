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

class TreeNode(object):
    def __init__(self, data):
        self.data = data
        self.children = []

    def add_child(self, obj):
        self.children.append(obj)

def validate_subtree(tree_str):
  assert tree_str.count("(") == tree_str.count(")")

def strip_parentheses(tree_str):
  return tree_str[1:-1]

def node_split(input_str):
  """ 
  Split str based on comma but go around interior nodes so 
  ((One,Two),(Three,Four)),Five would return a list of 
  [((One,Two),(Three,Four)), Five]
  """
  out_list = []; pos = 0
  while pos < len(input_str):
    curr_str = ""
    # Case 1: interior node 
    if input_str[pos] == "(":
      curr_str += input_str[pos]
      status = 1
      pos += 1

      while status != 0:
        if input_str[pos] == "(":
          status += 1
        elif input_str[pos] == ")":
          status -= 1
        curr_str += input_str[pos]
        pos += 1
      out_list.append(curr_str)
    # Case 2: comma
    elif input_str[pos] == ",":
      pos += 1
    # Case 3: tip
    else:
      while pos < len(input_str) and input_str[pos] != ",":
        curr_str += input_str[pos]
        pos += 1
      out_list.append(curr_str)
  return out_list

def find_max_depth_of_tree(curr_node):
  """ return depth of a tree"""
  if len(curr_node.children) == 0:
    return 0
  else:
    max_depth = 0
    for child in curr_node.children:
      max_depth = max(max_depth, find_max_depth_of_tree(child))
    return (max_depth+1)

def traversal_newick_tree(curr_node, traversal):
  """ performs a pre-order traversal of tree """
  traversal.append(curr_node.data)
  if len(curr_node.children) == 0:
    return traversal
  else:
    for child in curr_node.children:
      traversal = traversal_newick_tree(child, traversal)
  return traversal

def find_nodes_at_depth_N(curr_node, N, node_list, curr_depth):
  if curr_depth == N:
    node_list.append(curr_node.data)
    return node_list
  else:
      for child in curr_node.children:
        node_list = find_nodes_at_depth_N(child, N, node_list, curr_depth+1)
      return node_list

def get_leaf_nodes_of_tree(curr_node, node_list):
  if len(curr_node.children) == 0:
    node_list.append(curr_node.data)
    return node_list
  else:
    for child in curr_node.children:
      node_list = get_leaf_nodes_of_tree(child, node_list)
    return node_list

def build_newick_tree(curr_node, tree_str):
  """ builds the tree based on newick format """
  validate_subtree(tree_str)
  if "(" not in tree_str and ")" not in tree_str:
    children = tree_str.split(",")
    for child in children:
      curr_node.add_child(TreeNode(child))
  else:
    interior_children = node_split(tree_str)
    for child in interior_children:
      if "(" not in child:
        curr_node.add_child(TreeNode(child))
      else:
        curr_node.add_child(build_newick_tree(TreeNode("Interior"), strip_parentheses(child)))
  return curr_node

def build_ordered_genera_list(traversal_with_names, genera_list):
  """ Based on the traversal, reorganize the genera list """
  final_list = []
  for taxa in traversal_with_names:
    if taxa != "Interior" and taxa != "Root":
      if taxa in genera_list and taxa not in final_list:
        final_list.append(taxa)
  assert len(final_list) == len(genera_list)
  return final_list

def get_silva_tree(input_filepath):
    """ Loads the input in Newick Tree format """
    in_fd = open(input_filepath, "r")
    lines = [x.strip() for x in in_fd.readlines()]
    assert len(lines) == 1

    silva_tree_str = lines[0][:-2]
    return silva_tree_str

def load_id_to_name_map(input_filepath):
  """ Takes in map file to connect ID numbers to taxa names """
  in_fd = open(input_filepath, "r")
  lines = [x.strip() for x in in_fd.readlines()]

  id_map = {}
  for line in lines:
    line_split = line.split()
    # build the name, since some are multiple words
    name = line_split[1]; pos = 2
    while line_split[pos] != "-1":
      name += " " + line_split[pos]
      pos += 1
    id_map[int(line_split[0])] = name
  return id_map

def load_silva_headers(silva_filepath):
  """ Loads the SILVA database into memory """
  in_fd = open(silva_filepath, "r")
  all_lines = [x.strip() for x in in_fd.readlines() if ">" in x]
  return all_lines

def get_genera_list(silva_seqs):
  """ Take in a list of headers, and return all the genera """
  genera_list = set()
  eukar_list = set()
  for header in silva_seqs:
    # 1. For Bacteria/Archaea, we focus on sequences with 6 taxa levels
    if "Archaea;" in header or "Bacteria;" in header:
      header_split = header.split(";")
      if len(header_split) == 7 and header_split[5] != "uncultured":
        genera_list.add(header_split[5])
    # 2. For Eukaryota, we focused on sequences where genus is in organism name
    elif "Eukaryota" in header:
        header_split = header.split(";")
        if header_split[-2] in header_split[-1] and header_split[-2] != "uncultured":
          genera_list.add(header_split[-2])
  return list(genera_list)

def extract_genus_from_header(header):
  """ Takes in SILVA header line and returns genus, if it exists """
  # 1. For Bacteria/Archaea, we focus on sequences with 6 taxa levels
  if "Archaea;" in header or "Bacteria;" in header:
      header_split = header.split(";")
      if len(header_split) == 7 and header_split[5] != "uncultured":
          return header_split[5]
      else:
        return ""
  # 2. For Eukaryota, we focused on sequences where genus is in organism name
  elif "Eukaryota" in header:
      header_split = header.split(";")
      if header_split[-2] in header_split[-1] and header_split[-2] != "uncultured":
          return header_split[-2]
      else:
        return ""
  else:
    return ""
      
def main_use_order(args):
    """ main method, version that uses tree order to organize documents """

    # Step 1: Load in SILVA tree structure, determine a list of the genera
    #         ordered based on the tree structure.
    full_input_tree = get_silva_tree(args.tree_path)
    id_map = load_id_to_name_map(args.tree_map_path)
    silva_seqs = load_silva_headers(args.input_file)
    genera_list = get_genera_list(silva_seqs)

    # Build the SILVA tree, perform traversal, and generate ordered list
    root = build_newick_tree(TreeNode("Root"), strip_parentheses(full_input_tree))

    traversal = traversal_newick_tree(root, [])
    traversal_with_names = [id_map[int(x)] if x != "Interior" and x != "Root" else x for x in traversal]
    ordered_genera_list = build_ordered_genera_list(traversal_with_names, genera_list)

    # Step 2: Process the headers to identify how many sequences are in each genera
    count_dict = dict.fromkeys(ordered_genera_list, 0)
    zero_count = 0
    for header in silva_seqs:
        # 1. For Bacteria/Archaea, we focus on sequences with 6 taxa levels
        if "Archaea;" in header or "Bacteria;" in header:
            header_split = header.split(";")
            if len(header_split) == 7 and header_split[5] != "uncultured":
                count_dict[header_split[5]] += 1
        # 2. For Eukaryota, we focused on sequences where genus is in organism name
        elif "Eukaryota" in header:
            header_split = header.split(";")
            if header_split[-2] in header_split[-1] and header_split[-2] != "uncultured":
                count_dict[header_split[-2]] += 1

    # with open(args.output_dir + "output.csv", "w") as out_fd:
    #     for i, genus in enumerate(ordered_genera_list):
    #         out_fd.write(f"{i+1},{count_dict[genus]}\n")

    print(f"\n[log] input file had {len(silva_seqs)} seqeunces in it.")
    print(f"[log] we filtered it down to {sum([count_dict[x] for x in count_dict.keys()])} sequences")
    print(f"[log] number of genera found: {len(ordered_genera_list)}\n")

    # Step 3: Load all the sequences from the SILVA database, not just headers this time
    all_lines = []
    with open(args.input_file, "r") as in_fd:
        all_lines = [x.strip() for x in in_fd.readlines()]

    # Step 4: Write out sequences for each genera separately
    num_genera_to_write = min(len(ordered_genera_list), max(args.num_genera, 2))
    print(f"[log] writing out the separate FASTA files for top {num_genera_to_write} genera")

    for i, key in enumerate(ordered_genera_list):
        # Only focus on top n classes
        if i >= num_genera_to_write:
            break
        # Write out all the sequences for this genus
        with open(args.output_dir + f"tax_group_{i+1}.fna", "w") as out_fd:
            for j, line in enumerate(all_lines):
                if '>' in line:
                    curr_seq_genus = extract_genus_from_header(line)
                    if curr_seq_genus == key:
                        out_fd.write(f"{line}\n{all_lines[j+1]}\n")

    # Step 5: Create the filelist
    with open(args.output_dir + "filelist.txt", "w") as out_fd:
        for i in range(num_genera_to_write):
            out_fd.write(f"{args.output_dir}tax_group_{i+1}.fna {i+1}\n")
    print(f"[log] finished writing list with the paths\n")

def main(args):
    """ main method of the script (original version, does not use tree order)"""
    
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
        if i >= 1000:
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
        for i in range(1000):
            out_fd.write(f"{args.output_dir}tax_group_{i+1}.fna {i+1}\n")

def parse_arguments():
    """ Defines the command-line argument parser, and return arguments """
    parser = argparse.ArgumentParser(description="This script allows users to separate a SILVA database based on a "
                                                  "specific level of the taxonomy. ")
    parser.add_argument("-i", "--input", dest="input_file", required=True, help="input SILVA database")
    parser.add_argument("-o", "--output-dir", dest="output_dir", required=True, help="output directory for files")
    parser.add_argument("--use-order", dest="use_genera_order", action="store_true", default=False, help="use SILVA tree order to order the documents")
    parser.add_argument("--tree", dest="tree_path", default="", required=False, help="path to *.tre file for the SILVA database")
    parser.add_argument("--tree-map", dest="tree_map_path", default="", required=False, help="path to *.map file for the SILVA database")
    parser.add_argument("-n", dest="num_genera", default=10, type=int, help="number of genera to write out to separate files")
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
    
    if args.use_genera_order:
        if not os.path.isfile(args.tree_path):
            print("Error: path to *.tre file is not valid")
            exit(1)
        if not os.path.isfile(args.tree_map_path):
            print("Error: path to *.map file is not valid")
            exit(1)
        
        if ".tre" not in args.tree_path:
            print("Error: the tree file does not have valid extension (*.tre)")
            exit(1)
        if ".map" not in args.tree_map_path:
            print("Error: the tree map file does not have valid extension (*.map)")
            exit(1)

if __name__ == "__main__":
    args = parse_arguments()
    check_arguments(args)
    if not args.use_genera_order:
        main(args)
    else:
        main_use_order(args)