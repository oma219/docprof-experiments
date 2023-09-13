#!/usr/bin/env python3

# Name: plot_lca_on_tree_exp11.py
# Description: This is a python script is used by experiment 11 
#              to take in a taxonomy, and plot the LCAs of each MEM
#              on the tree.
# Date: September 13th, 2022

import os
import argparse
import math
import csv
import pysam
import matplotlib as mpl
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np

class TreeNode:
    def __init__(self, id, name, parent=None):
        self.id = id
        self.name = name
        self.children = []
        self.parent = parent
        self.weight = 0
    def add_child(self, child):
        self.children.append(child)
    def add_weight(self, weight):
        self.weight += weight

def get_tree_node_from_list(tree_nodes, name):
    """
    Input: tree_nodes - list of tree nodes in tree
           name - name of node
    Output: TreeNode or None - depending on if that name exists in tree
    """
    for node in tree_nodes:
        if name == node.name:
            return node
    return None

def get_tree_node_from_list_with_num(tree_nodes, num):
    """
    Input: tree_nodes - list of tree nodes in tree
           num - id number of node
    Output: TreeNode or None - depending on if that name exists in tree
    """
    for node in tree_nodes:
        if num == node.id:
            return node
    return None

def build_taxonomy_from_file(path):
    """
    Input: file - path to taxonomy files
    Output: root node of tree, list of nodes
    """
    tree_nodes = []
    edge_list = []
    node_labels = {}
    id_num = 0
    with open(path, "r") as in_fd:
        for line in in_fd:
            line_split = line.strip().split("->")
            parent_num = -1
            
            # Get parent node name
            parent_name = line_split[0].strip()

            # Check if parent already exists, to avoid recreating it
            if get_tree_node_from_list(tree_nodes, parent_name) != None:
                parent_node = get_tree_node_from_list(tree_nodes, parent_name)
            else:
                parent_node = TreeNode(id_num, parent_name)
                node_labels[id_num] = parent_name
                tree_nodes.append(parent_node)
                id_num += 1
            parent_num = parent_node.id

            # Build the children nodes
            child_names = [x.strip() for x in line_split[1].split()]
            for name in child_names:
                # Check if child already exists
                if get_tree_node_from_list(tree_nodes, name) != None:
                    curr_node = get_tree_node_from_list(tree_nodes, name)
                else:
                    curr_node = TreeNode(id_num, name, parent=parent_node)
                    node_labels[id_num] = name
                    id_num += 1
                    tree_nodes.append(curr_node)
                edge_list.append((parent_num, curr_node.id))

                # Update parent field for child
                parent_node.add_child(curr_node)
    return tree_nodes, edge_list, node_labels

def build_doc_to_name_dict(doc_to_name_file):
    """
    Input: doc_to_name_file - path to file
    Output: dictionary - maps node names to doc numbers
    """
    name_dict = {}
    with open(doc_to_name_file, "r") as in_fd:
        for line in in_fd:
            line_split = [x.strip() for x in line.split()]
            assert len(line_split) == 2

            name_dict[int(line_split[1])] = line_split[0]
    return name_dict

def extact_doc_listings(lengths_file):
    """
    Input: lengths_file - path to output from from docprofiles
    Output: list of tuples - each tuple contains read length, leftmost, 
            rightmost occurrence
    """
    read_data = []
    with open(lengths_file, "r") as in_fd:
        read_name = ""
        read_length = 0
        for line in in_fd:
            if '>' in line:
                read_name = line.strip()
                read_length = int(line.split("_")[-1])
            else:
                listings_results = line.split()[1].strip().split(",")
                assert len(listings_results) == 2

                leftmost_occ = int(listings_results[0].replace("{", ""))
                rightmost_occ = int(listings_results[1].replace("}", ""))
                read_data.append([read_length, leftmost_occ, rightmost_occ])
    return read_data

def travel_to_root(curr_node):
    """
    Input: curr_node - leaf node that string was found in
    Output: list of node names up to root
    """
    node_list = []
    while curr_node != None:
        node_list.append(curr_node.name)
        curr_node = curr_node.parent
    return node_list

def compute_lca_from_traversals(left_to_root, right_to_root):
    """
    Input: two lists - traversals from leaf node up to root
    Output: name of LCA node
    """
    for name in left_to_root:
        if name in right_to_root:
            return name

    print("Error: Issue occurred during LCA computation.")
    exit(1)

def compute_lca_for_each_read(tree_nodes, doc_to_name, read_data):
    """
    Input: read_data - list of tuples with read length, and leftmost, right most occurrence
    Output: list of nodes with updated weights 
    """
    for length, left, right in read_data:
        # Get names from doc numbers
        left_name = doc_to_name[left]
        right_name = doc_to_name[right]

        # Get the nodes for each name, and traverse to root
        left_node = get_tree_node_from_list(tree_nodes, left_name)
        right_node = get_tree_node_from_list(tree_nodes, right_name)
        left_to_root = travel_to_root(left_node)
        right_to_root = travel_to_root(right_node)

        # Get the LCA node, and increment the weight
        lca_node_name = compute_lca_from_traversals(left_to_root, right_to_root)
        curr_node = get_tree_node_from_list(tree_nodes, lca_node_name)
        curr_node.add_weight(length)
    return tree_nodes

def make_graph(edge_list, node_sizes, plot_title, output_file_path):
    # Create graph, create edges, and spread out nodes
    seed = 13648
    G = nx.Graph() 
    for left, right in edge_list:
        G.add_edge(left, right)
    pos = nx.spring_layout(G, seed=seed)

    # Initialize size of nodes
    M = G.number_of_nodes()
    node_colors = node_sizes

    # Build nodes, edges
    cmap = plt.cm.plasma
    nodes = nx.draw_networkx_nodes(
        G, 
        pos, 
        node_size=300, 
        node_color=node_colors,
        cmap=cmap
        )
    edges = nx.draw_networkx_edges(
        G,
        pos,
        node_size=node_sizes,
        arrowstyle="->",
        arrowsize=10,
        edge_color="indigo",
        # edge_cmap=cmap,
        width=2,
    )

    # Save plot
    pc = mpl.collections.PatchCollection([], cmap=cmap)
    pc.set_array(node_colors)

    ax = plt.gca()
    ax.set_axis_off()
    plt.colorbar(pc, ax=ax)
    plt.show()

    plt.title(plot_title)
    plt.savefig(output_file_path, dpi=800)
    plt.clf()

def main(args):
    # Build tree, and get mapping from names to doc numbers
    tree_nodes, edge_list, node_labels = build_taxonomy_from_file(args.tax_file)
    doc_to_name = build_doc_to_name_dict(args.doc_to_name_file)

    for doc_num in range(args.num_docs):
        # Extract read lengths, and their leftmost, rightmost occurrence
        read_data = extact_doc_listings(args.input_files[doc_num])
        tree_nodes = compute_lca_for_each_read(tree_nodes, doc_to_name, read_data)

        # Generate node sizes based on weight
        node_sizes = [0 for i in range(len(tree_nodes))]
        node_size_and_name = []
        for node in tree_nodes:
            node_sizes[node.id] = node.weight * (1.0/1000000.0)
            node_size_and_name.append((node.name, node.weight))
        make_graph(edge_list, node_sizes, doc_to_name[doc_num], f"{args.output_dir}doc_{doc_num+1}_plots.png")
        
        # Write weights to file
        with open(f"{args.output_dir}doc_{doc_num+1}_results.txt", "w") as out_fd: 
            node_size_and_name.sort(key=lambda x: x[1])
            for tup in node_size_and_name:
                out_fd.write(f"{tup}\n")
        
        # Clear the weight of all the nodes
        for curr_node in tree_nodes:
            curr_node.weight = 0
        
def parse_arguments():
    """ Defines the command-line argument parser, and return arguments """
    parser = argparse.ArgumentParser(description="This script is for experiment 11, plotting LCA on taxonomy.")
    parser.add_argument("-t", "--taxonomy", dest="tax_file", required=True, type=str)
    parser.add_argument("-i", "--input", dest="input_files", required=True, type=str, nargs='*')
    parser.add_argument("-n", "--num-of-docs", dest="num_docs", required=True, type=int)
    parser.add_argument("-d", "--doc-to-name", dest="doc_to_name_file", required=True, type=str)
    parser.add_argument("-o", "--output-dir", dest="output_dir", required=True, type=str)
    args = parser.parse_args()
    return args

def check_arguments(args):
    """ Checks for invalid arguments """
    if not os.path.isfile(args.tax_file):
        print("Error: path to taxonomy is not valid")
        exit(1)
    for path in args.input_files:
        if not os.path.isfile(path):
            print("Error: at least one of the input file path is not valid.")
            exit(1)
    if args.num_docs != len(args.input_files):
        print("Error: number of documents does not match input file list")
        exit(1)
    if not os.path.isfile(args.doc_to_name_file):
        print("Error: path to doc_to_name file is not valid")
        exit(1)
    if not os.path.isdir(args.output_dir):
        print("Error: path to output directory is not valid.")
        exit(1)
    if args.output_dir[-1] != "/":
        args.output_dir += "/"

if __name__ == "__main__":
    args = parse_arguments()
    check_arguments(args)
    main(args)