"""
Purpose: Script to classify the 16S rRNA reads
         that were processed by pfp_doc64
"""

import os
import sys
import argparse
import math
from pprint import pprint

class SILVA_taxonomy:
    def __init__(self,
                 silva_taxonomy_file):
        """
        Input Fields:
            silva_taxonomy_file: SILVA file that maps every taxonomic
                                 traversal to node id
        """
        self.raw_data = []
        self.initialize_raw_data(silva_taxonomy_file)

        self.trav_to_rank = {}
        for _, trav, rank in self.raw_data:
            self.trav_to_rank[trav] = rank
        
        self.id_to_trav = {}
        for node_id, trav, _ in self.raw_data:
            self.id_to_trav[int(node_id)] = trav

    def initialize_raw_data(self, 
                            input_path):
        """
        Initializes the raw_data attribute with
        tuples containing the node_id, traversal, and
        taxonomic rank
        """
        with open(input_path, "r") as in_fd:
            for line in in_fd:
                last_pos_of_traversal, traversal = SILVA_taxonomy.get_traversal_from_full_line(line)
                
                line_split = line.split()
                node_id = int(line_split[last_pos_of_traversal+1])
                rank = line_split[last_pos_of_traversal+2]

                self.raw_data.append([node_id, traversal, rank])

    def get_traversal_from_full_line(line):
        """
        Get the traversal part of line, the presence of spaces in
        in traversal string makes it slightly less than trivial

        line: string with entire lines contents from silva taxonomy file
        """
        line_split = line.split()
        last_pos_of_traversal = max([i for i in range(len(line_split)) if ";" in line_split[i]])
        traversal = ' '.join(line_split[0:last_pos_of_traversal+1])
        
        assert (len(line_split) - last_pos_of_traversal) <= 4   
        return last_pos_of_traversal, traversal

    def get_rank_from_traversal(self,
                                traversal):
        """ 
        Returns taxonomic rank of a traversal path if
        it exists, otherwise None
        """
        if traversal not in self.trav_to_rank:
            return None 
        else:
            return self.trav_to_rank[traversal]

    def check_traversal_presence(self, 
                                 traversal):
        """
        Returns True if traversal is present in SILVA
        file, and False otherwise
        """
        if traversal in self.trav_to_rank:
            return True 
        else:
            return False

    def get_silva_dict_for_level_to_list(self):
        """
        Creates a dictionary thats takes a taxonomic rank like
        'domain', 'phylum' and returns a list of the names in the
        SILVA taxonomy. 
        """
        tree_levels = ["domain", "phylum", "class", "order", "family", "genus"]
        silva_ranks_dict = {}

        for level in tree_levels:
            silva_level_names = []
            for tup in self.raw_data:
                if tup[2] == level:
                    name = tup[1].split(";")[-2]
                    silva_level_names.append(TaxonomyTraversal(tup[1], self))
            silva_ranks_dict[level] = silva_level_names
        return silva_ranks_dict

    def get_traversal_from_id(self,
                              node_id):
        """ 
        Returns traversal path based on node id if
        it exists, otherwise None
        """
        if node_id not in self.id_to_trav:
            return None 
        else:
            return self.id_to_trav[int(node_id)]

class TaxonomyTraversal:
    def __init__(self,
                traversal,
                silva_taxonomy):
        """
        Input Fields:
            traversal: a string from the doc_to_traversal file
            silva_taxonomy_file: SILVA file that maps every traversal to rank
        """
        self.trav_str = traversal
        self._domain = None
        self._phylum = None
        self._class = None
        self._order = None
        self._family = None
        self._genus = None

        self.initialize_ranks(traversal, silva_taxonomy)
    
    def initialize_ranks(self,
                         traversal,
                         silva_taxonomy):
        """ parse the traversal and fill in taxonomic ranks """
        assert silva_taxonomy.check_traversal_presence(traversal)
        curr_str = ""
        for ch in traversal:
            curr_str += ch
            if ch == ";":
                rank = silva_taxonomy.get_rank_from_traversal(curr_str)
                if rank != None:
                    if rank == "domain":
                        self._domain = curr_str
                    elif rank == "phylum":
                        self._phylum = curr_str
                    elif rank == "class":
                        self._class = curr_str
                    elif rank == "order":
                        self._order = curr_str
                    elif rank == "family":
                        self._family = curr_str
                    elif rank == "genus":
                        self._genus = curr_str

    def get_certain_level(self,
                         level_name):
        """
        Return the specific level requested if it 
        is present, if not present return None
        """
        assert level_name in ["domain", "phylum", "class", "order", "family", "genus"]
        curr_trav = ""

        if level_name == "domain":
            curr_trav = self._domain
        elif level_name == "phylum":
            curr_trav = self._phylum
        elif level_name == "class":
            curr_trav = self._class
        elif level_name == "order":
            curr_trav = self._order
        elif level_name == "family":
            curr_trav = self._family
        elif level_name == "genus":
            curr_trav = self._genus
        
        if curr_trav is None:
            return None
        else:
            return curr_trav.split(";")[-2]

class ReadsetTraversal:
    def __init__(self,
                 traversal):
        """
        Input Fields:
            traversal: a string from the *.tab file
        """
        self.trav_str = traversal 
        self._domain = None
        self._phylum = None 
        self._class = None 
        self._order = None 
        self._family = None 
        self._genus = None

        self.initialize_ranks()
            
    def initialize_ranks(self):
        """ initialize fields for different levels of taxonomy """
        self._domain = self.extract_specific_level_from_read_set("domain")
        self._phylum = self.extract_specific_level_from_read_set("phylum")
        self._class = self.extract_specific_level_from_read_set("class")
        self._order = self.extract_specific_level_from_read_set("order")
        self._family = self.extract_specific_level_from_read_set("family")
        self._genus = self.extract_specific_level_from_read_set("genus")

    def extract_specific_level_from_read_set(self, level):
        """ parse out the specific level from traversal """
        level_to_info = {"domain": ("sk__", 0), 
                        "phylum": ("p__", 2),
                        "class": ("c__", 3),
                        "order": ("o__", 4),
                        "family": ("f__", 5),
                        "genus": ("g__", 6)}
        assert level in level_to_info
        prefix, relevant_pos = level_to_info[level]

        trav_split = self.trav_str.split(";")
        assert prefix in trav_split[relevant_pos]

        return trav_split[relevant_pos][len(prefix):]

    def get_certain_level(self,
                         level_name):
        """
        Return the specific level requested if it 
        is present, if not present return None
        """
        assert level_name in ["domain", "phylum", "class", "order", "family", "genus"]
        curr_trav = ""

        if level_name == "domain":
            return self._domain
        elif level_name == "phylum":
            return self._phylum
        elif level_name == "class":
            return self._class
        elif level_name == "order":
            return self._order
        elif level_name == "family":
            return self._family
        elif level_name == "genus":
            return self._genus

class Read:
    def __init__(self,
                 name,
                 mate1_listings,
                 mate2_listings,
                 kraken2_result,
                 correct_taxa):
        """
        Input Fields:
            name: name of the read
            mate{1,2}_listings: output from PFP_DOC64
            correct_taxa: ReadsetTraversal Object that represents the correct classification for this read
        """
        self.name = name 
        self.mate1_listings = mate1_listings 
        self.mate2_listings = mate2_listings
        self.kraken2_result = kraken2_result
        self.correct_taxa = correct_taxa
    
    def classify_with_pfpdoc(self,
                             level,
                             list_of_classes,
                             doc_id_to_traversal):
        """ 
        Classifies the read using pfp_doc results at
        a certain level of tree. 

        level: string, that tells you what level of tree we are classifying at
        list_of_classes: list of TaxonomyTraversal objects that are all possible classes
        doc_id_to_traversal: dictionary from node ids to TaxonomyTraversal objects

        Returns True if Correct, and False Otherwise
        """
        # make sure that we only have one correct class
        correct_name = self.correct_taxa.get_certain_level(level)
        classes_with_correct_name = [(x.get_certain_level(level) == correct_name) for x in list_of_classes]
        assert sum(classes_with_correct_name) == 1

        # create a dictionary of names to weights
        name_to_weight = {}
        for x in list_of_classes:
            name_to_weight[x.get_certain_level(level)] = [0, 0]
        
        # helper method to generate weight vector
        def get_weight_vector(num_values):
            if num_values == 1:
                return [1.0]
            elif num_values == 2:
                return [0.5, 0.5]
            else:
                final_array = [0.0 for i in range(num_values)]
                remaining_weight = 1.0

                # handle odd length scenarios
                if num_values % 2 == 1:
                    final_array[int(num_values/2)] = round(1.0/(num_values-1), 4)
                    remaining_weight -= 1.0/(num_values-1)

                for i in range(int(num_values/2)-1, -1, -1):
                    weight = remaining_weight * 1.0/(num_values-1)
                    final_array[i] = round(weight,4)
                    final_array[num_values-i-1] = round(weight,4)
                    remaining_weight -= 2 * weight 

                add_on_top = remaining_weight/num_values
                return [round(x + add_on_top, 4) for x in final_array]

        # helper method to add to dictionary
        def add_mate_listings_to_weight_vector(mate_listings, weight_dict):
            mate_listings_list = mate_listings.split()
            assert len(mate_listings_list) % 2 == 0

            for i in range(0, len(mate_listings_list), 2):
                # get the length of exact match
                start_and_stop = mate_listings_list[i][1:-1]
                stop = int(start_and_stop.split(",")[1])
                start = int(start_and_stop.split(",")[0])
                length = stop - start + 1
                assert length > 0

                # get the traversal of documents that are hit
                doc_listing = [int(x) for x in mate_listings_list[i+1][1:-1].split(",")]
                weight_vec = get_weight_vector(len(doc_listing))
                for doc_id in doc_listing:
                    if doc_id in doc_id_to_traversal:
                        taxa_hit = doc_id_to_traversal[doc_id].get_certain_level(level)
                        if taxa_hit in weight_dict:
                            weight_dict[taxa_hit][0] += length
                            weight_dict[taxa_hit][1] += 1
            return weight_dict
        
        # add each mate1 results to list
        name_to_weight = add_mate_listings_to_weight_vector(self.mate1_listings, name_to_weight)
        name_to_weight = add_mate_listings_to_weight_vector(self.mate2_listings, name_to_weight)

        # print(self.name)
        # print(self.mate1_listings)
        # print(self.mate2_listings)

        # count = 0
        # for k, v in sorted(name_to_weight.items(), key=lambda item: item[1], reverse=True):
        #     print(f"{k}: {v}")
        #     count += 1
        #     if count == 10:
        #         break

        # print(f"correct taxa: {correct_name}")

        # get the predicted taxa, and return results
        max_count = max(name_to_weight.values(), key=lambda x: x[0])[0]
        max_name_taxa = [k for (k, v) in name_to_weight.items() if v[0] == max_count]
        
        # handle the case of uncultured genera
        if max_name_taxa == ["uncultured"]:
            name_to_weight.pop("uncultured")
            max_count = max(name_to_weight.values(), key=lambda x: x[0])[0]
            max_name_taxa = [k for (k, v) in name_to_weight.items() if v[0] == max_count] 

        # print(max_name_taxa)
        # print(max_count)
        #predicted_taxa = max(name_to_weight, key=name_to_weight.get)
        if any([correct_name in x for x in max_name_taxa]):
            return [True, max_name_taxa[0]]
        else:
            return [False, max_name_taxa[0]]

    def classify_with_kraken(self,
                             level,
                             silva_tax_obj):
        """
        Classifies using Kraken2 result, and then
        returns whether read is TP, VP, FN or FP

        level: the level of taxonomy we want to classify at
        silva_tax_obj: SILVA_taxonomy object
        """
        kraken2_output = self.kraken2_result.split()
        if kraken2_output[0] == "U":
            return ["FN", None]

        node_id = int(kraken2_output[2])
        traversal = silva_tax_obj.get_traversal_from_id(node_id)
        assert traversal is not None

        trav_obj = TaxonomyTraversal(traversal, silva_tax_obj)
        kraken2_name = trav_obj.get_certain_level(level)
        correct_name = self.correct_taxa.get_certain_level(level)

        if kraken2_name == correct_name:
            return ["TP", correct_name]
        elif kraken2_name is not None:
            return ["FP", kraken2_name]

        # check if classification is vaguely correct (only if
        # there is a match at a level below domain)
        tree_levels = ["genus", "family", "order", "class", "phylum", "domain"]
        assert level in tree_levels
        start_index = tree_levels.index(level)

        for pos in range(start_index+1, len(tree_levels)-1):
            kraken2_name = trav_obj.get_certain_level(tree_levels[pos])
            correct_name = self.correct_taxa.get_certain_level(tree_levels[pos])
            if kraken2_name == correct_name:
                return ["VP", None]
        return ["FP", None]

class ReadSet:
    def __init__(self,
                mate1_listings,
                mate2_listings,
                doc_to_traveral_file,
                silva_traversal_file,
                readset_traversal_file,
                kraken2_read_file,
                bracken_output):
        """
        Input Fields:
            mate{1,2}_listings: output files from PFP_DOC64
            doc_to_traversal_file: index file from PFP_DOC64
            silva_traversal_file: taxonomy file from SILVA database
            readset_traversal_file: metadata file from Alexandre et al, 2018
            kraken2_read_file: classification file from Kraken2
        """
        self.paths = {}
        self.initialize_paths_variable(mate1_listings, 
                                       mate2_listings, 
                                       doc_to_traveral_file, 
                                       silva_traversal_file,
                                       readset_traversal_file,
                                       kraken2_read_file,
                                       bracken_output)
        ReadSet.log_msg("finished paths initialization.")

        self.silva_taxonomy = SILVA_taxonomy(self.paths["silva_traversal"])
        ReadSet.log_msg("finished SILVA initialization.")

        self.doc_id_to_full_traversal = {}
        self.initialize_pfpdoc_doc_to_traversal_dict(self.paths["doc_id_to_traversal"])
        ReadSet.log_msg("finished doc_id initialization.")

        self.read_name_to_traversal = {}
        self.initialize_readset_truthset_dict(self.paths["readset_truthset"])
        ReadSet.log_msg("finished truthset initialization.")

        self.silva_nodes_at_each_level = self.silva_taxonomy.get_silva_dict_for_level_to_list()
        ReadSet.log_msg("finished loading node lists.")

        self.read_list = []
        self.load_reads(self.paths["mate1_listings"],
                        self.paths["mate2_listings"],
                        self.paths["kraken2_read_class"])
        ReadSet.log_msg(f"finished loading {len(self.read_list)} read pairs")

        # # interesting reads: A500_V12_443-200K28 (3982), A500_V12_366-200K592(3325)
        # print(self.read_list[4150].classify_with_pfpdoc("order", 
        #                                           self.silva_nodes_at_each_level["order"],
        #                                           self.doc_id_to_full_traversal))
        # exit(1)
        
    def initialize_paths_variable(self,
                                  mate1_listings,
                                  mate2_listings,
                                  doc_to_traveral_file,
                                  silva_traversal_file,
                                  readset_traversal_file,
                                  kraken2_read_class,
                                  bracken_output):
        """ initialize a dictionary with all the file paths """
        self.paths["mate1_listings"] = mate1_listings
        self.paths["mate2_listings"] = mate2_listings
        self.paths["doc_id_to_traversal"] = doc_to_traveral_file
        self.paths["silva_traversal"] = silva_traversal_file
        self.paths["readset_truthset"] = readset_traversal_file
        self.paths["kraken2_read_class"] = kraken2_read_class
        self.paths["bracken_output"] = bracken_output

    def initialize_pfpdoc_doc_to_traversal_dict(self,
                                                file_path):
        """ parse the two-column file that maps doc # to genus """
        with open(file_path, "r") as in_fd:
            for line in in_fd:
                line_split = [x.strip() for x in line.split()]
                curr_traversal = " ".join(line_split[1:])
                self.doc_id_to_full_traversal[int(line_split[0])-1] = TaxonomyTraversal(curr_traversal,
                                                                                        self.silva_taxonomy)

    def initialize_readset_truthset_dict(self,
                                         input_path):
        """ parse the two-column file into dictionary mapping read name to traversal """
        with open(input_path, "r") as in_fd:
            for line in in_fd:
                line_split = [x.strip() for x in line.split()]
                assert len(line_split) == 2
                self.read_name_to_traversal[line_split[0]] = ReadsetTraversal(line_split[1])

    def load_reads(self,
                  mate1_path,
                  mate2_path,
                  kraken2_path):
        """ load the results into Read objects """
        with open(mate1_path, "r") as mate1_fd, open(mate2_path, "r") as mate2_fd, open(kraken2_path, "r") as kraken2_fd:
            mate1_lines = [x.strip() for x in  mate1_fd.readlines()]
            mate2_lines = [x.strip() for x in  mate2_fd.readlines()]
            kraken2_lines = [x.strip() for x in kraken2_fd.readlines()]
        
        header_m1 = ""; listing_m1 = ""; 
        header_m2 = ""; listing_m2 = ""; pos = 0
        curr_kraken_line = ""; count = 0

        for mate1_line, mate2_line in zip(mate1_lines, mate2_lines):
            if mate1_line.startswith(">") and mate2_line.startswith(">"):
                header_m1 = mate1_line[1:]; header_m2 = mate2_line[1:]
                pos += 1
            elif pos == 1:
                listing_m1 = mate1_line; listing_m2 = mate2_line 
                pos = 0

                # Grab the corresponding kraken2 line ...
                curr_kraken_line = kraken2_lines[count]
                count += 1

                # Save read object ...
                assert header_m1.split("/")[0] == header_m2.split("/")[0]
                assert header_m1.split("/")[0] == curr_kraken_line.split()[1]

                read_group_name = header_m1.split("-")[0]
                assert read_group_name in self.read_name_to_traversal

                self.read_list.append(Read(header_m1, 
                                           listing_m1, 
                                           listing_m2,
                                           curr_kraken_line,
                                           self.read_name_to_traversal[read_group_name]))

    def generate_sensitivity_plot(self):
        """ generate an output file with sensitivity vs major clade """
        major_clades = ["genus", "family", "order", "class", "phylum", "domain"]

        # Step 1: analyze docprofiles
        for clade in major_clades:
            results = []
            for i, read_obj in enumerate(self.read_list):
                results.append(read_obj.classify_with_pfpdoc(clade, 
                                            self.silva_nodes_at_each_level[clade],
                                            self.doc_id_to_full_traversal)[0])
            num_tp = sum(results)
            sensitivity = round(num_tp/len(results), 4)
            print(f"docprofiles,{clade},{sensitivity}")
        
        # Step 2: analyze the kraken results
        for clade in major_clades:
            results = []
            for i, read_obj in enumerate(self.read_list):
                results.append(read_obj.classify_with_kraken(clade,
                                                             self.silva_taxonomy)[0]) 
            assert len(results) == (results.count("TP") + results.count("FP") + results.count("VP") + results.count("FN"))

            num_tp = results.count("TP")
            sensitivity = round(num_tp/len(results), 4)
            print(f"kraken2,{clade},{sensitivity}")
    
    def generate_abundance_plot(self):
        """ generate an abundance plot at genus level """
        true_genera = {}
        for i, read_obj in enumerate(self.read_list):
            curr_genus = read_obj.correct_taxa.get_certain_level("genus")
            if curr_genus not in true_genera:
                true_genera[curr_genus] = 1
            else:
                true_genera[curr_genus] += 1

        total_sum = sum(true_genera.values())
        for key in true_genera:
            print(f"Truth,{key},{round(true_genera[key]/total_sum, 4)}")

        # get the genera abundances for docprofiles method
        docprofiles_generas = {"Others": 0}
        for key in true_genera:
            docprofiles_generas[key] = 0

        for i, read_obj in enumerate(self.read_list):
            curr_genus = read_obj.classify_with_pfpdoc("genus", 
                                            self.silva_nodes_at_each_level["genus"],
                                            self.doc_id_to_full_traversal)[1]
            if curr_genus in docprofiles_generas:
                docprofiles_generas[curr_genus] += 1
            else:
                docprofiles_generas["Others"] += 1

        total_sum = sum(docprofiles_generas.values())
        for key in docprofiles_generas:
            print(f"docprofiles,{key},{round(docprofiles_generas[key]/total_sum, 4)}")

        # get the genera abundances for kraken2 method
        kraken2_generas = {"Others": 0}
        for key in true_genera:
            kraken2_generas[key] = 0

        with open(self.paths["bracken_output"], "r") as in_fd:
            bracken_lines = [x.strip() for x in in_fd.readlines()][1:]
            for line in bracken_lines:
                line_split = line.split("\t")
                assert len(line_split) == 7

                curr_genus = line_split[0]
                num_reads = int(line_split[5])

                if curr_genus in kraken2_generas:
                    kraken2_generas[curr_genus] += num_reads
                else:
                    kraken2_generas["Others"] += num_reads

        total_sum = sum(kraken2_generas.values())
        for key in kraken2_generas:
            print(f"kraken2,{key},{round(kraken2_generas[key]/total_sum, 4)}")
        
        # compute the bray-curtis dissimilarity
        def compute_bray_curtis_dissimilarity(true_comp, exp_comp):
            Cij = 0
            for key in true_comp:
                if true_comp[key] > 0 and exp_comp[key] > 0:
                    Cij += min(true_comp[key], exp_comp[key])
            Si = sum(true_comp.values())
            Sj = sum(exp_comp.values())

            bc = 1.0 - (2 * Cij)/(Si + Sj)
            return round(bc, 4)
            
        bc_dist1 = compute_bray_curtis_dissimilarity(true_genera, docprofiles_generas)
        bc_dist2 = compute_bray_curtis_dissimilarity(true_genera, kraken2_generas)

        print(f"Bray-Curtis dissimilarity (truth-to-docprofiles) = {bc_dist1}")
        print(f"Bray-Curtis dissimilarity (truth-to-kraken2) = {bc_dist2}")

    def log_msg(msg):
        """ prints out log message """
        print("\033[0;32m[log]\033[00m " + msg)

################################################
# Main method
################################################

def main(args):
    # Create the readset object ...
    rs = ReadSet(args.mate1_listings,
                 args.mate2_listings,
                 args.doc_to_traversal_file,
                 args.silva_tax_file,
                 args.readset_truthset,
                 args.kraken2_read_class,
                 args.bracken_output)
    
    # rs.generate_sensitivity_plot()
    rs.generate_abundance_plot()

################################################
# Parse command-line options and check them ...
################################################

def parse_arguments():
    parser = argparse.ArgumentParser(description="classify the 16S rRNA classifications ...")
    parser.add_argument("--mate1-listings", dest="mate1_listings", help="path to mate1 listings file", required=True)
    parser.add_argument("--mate2-listings", dest="mate2_listings", help="path to mate2 listings file", required=True)
    parser.add_argument("--doc-id-to-traversal", dest="doc_to_traversal_file", help="path to doc_to_traversal.txt file", required=True)
    parser.add_argument("--silva-tax-ranks", dest="silva_tax_file", help="path to silva taxonomy file", required=True)
    parser.add_argument("--readset-truthset", dest="readset_truthset", help="path to truth ids for readset", required=True)
    parser.add_argument("--output-dir", dest="output_dir", help="path to output_dir", required=True)
    parser.add_argument("--kraken2-read-class", dest="kraken2_read_class", help="path to kraken2 read classifications", required=True)
    parser.add_argument("--bracken-output", dest="bracken_output", help="path to bracken output", required=True)
    args = parser.parse_args()
    return args

def check_args(args):
    if not os.path.isfile(args.mate1_listings):
        print("Error: input file 1 listings provided is not valid.")
        exit(1)
    if not os.path.isfile(args.mate2_listings):
        print("Error: input file 2 listings provided is not valid.")
        exit(1)
    if not os.path.isdir(args.output_dir):
        print("Error: output directory is not valid.")
        exit(1)
    elif args.output_dir[-1] != "/":
        args.output_dir += "/"
    if not os.path.isfile(args.doc_to_traversal_file):
        print("Error: doc_to_traversal provided is not valid.")
        exit(1)
    if not os.path.isfile(args.silva_tax_file):
        print("Error: silva taxonomy file provided is not valid.")
        exit(1)
    if not os.path.isfile(args.readset_truthset):
        print("Error: read set truthset provided is not valid.")
        exit(1)
    if not os.path.isfile(args.kraken2_read_class):
        print("Error: kraken2 read classification file is not valid.")
        exit(1)
    if not os.path.isfile(args.bracken_output):
        print("Error: bracken read classification file is not valid.")
        exit(1)

if __name__ == "__main__":
    args = parse_arguments()
    check_args(args)
    main(args)