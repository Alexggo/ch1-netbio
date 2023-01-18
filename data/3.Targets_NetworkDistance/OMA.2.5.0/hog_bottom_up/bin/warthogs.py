#!/usr/bin/env python
import sys
import argparse
import time as time
from Bio import Phylo
from gethogs import file_manager, lib, version
from gethogs.settings import Settings
from gethogs.statistic_tracker import StatisticTracker


def init_configuration():

    parser = argparse.ArgumentParser(description="warthogs - running GETHOGs",
                                     epilog="""
IMPORTANT:
  - Pairwise orthologous/paralogous relations files (standalone): No *-* inside you species name.
  - Pairwise orthologous/paralogous relations files (oma): No *.* inside you species name.
  - Parameter_1 (pair + update): 0 < value <= 100. 
"""
                                     )
    parser.add_argument('-i', '--input', required=True,
                        help='Path to input folder containing the pairwise orthologs')
    parser.add_argument('-k', '--paralogs', required=False,
                        help='Path to input folder containing the pairwise paralogs')
    parser.add_argument('-t', '--type', choices=("standalone", "oma"), default="standalone",
                        help="Type of structure of the orthologous pairwise relations folder, "
                             "*standalone* if all files are contained in a *PairwiseOrthologs* "
                             "folder (files: species_name1-species_name2.extension) or *oma* if "
                             "all files are contained in a nested oma-type folder structures "
                             "(dir:species_name1, files inside:species_name2.extension)")
    parser.add_argument('-m', '--merge_method', default="pair", choices=("pair", "update"),
                        help="The method used to clean the orthology graph during the "
                             "reconstruction of the ancestral HOGs. You can select *pair* "
                             "for cleaning the edges pairs of HOGs by pairs (-p required to "
                             "specify the minimal percent of relations mandatory of be kept) or "
                             "*update* if you want to update the orthology graph each time "
                             "you modify a pair of hogs (-p required to specify the minimal "
                             "percent of relations mandatory of be kept)")
    parser.add_argument('-p', '--parameter', type=int, required=True,
                        help="Parameter used by the merge_method")
    parser.add_argument('-o', '--output', required=True, type=str,
                        help="Output file name of orthoxml output file")
    parser.add_argument('-s', '--species_tree', required=True, type=str,
                        help="Path to a newick file with the species tree use as skeleton "
                             "for traversal. The leaves must match the genome names of the "
                             "--input data")
    parser.add_argument('-g', '--genome_info', type=str,
                        help="Genome Info file (optional) with order and number of genes")
    parser.add_argument('-u', '--unmerged',
                        help="maximum number of unmerged before freezing an HOGs")
    parser.add_argument('-d', '--dynamic',
                        help="merge threshold of the root used by the dynamic propagation")
    parser.add_argument('--version', action="version", version=version())
    conf = parser.parse_args()
    for opt, value in vars(conf).items():
        if value is None:
            continue
        if opt == 'input':
            Settings.set_pairwise_folder(str(value))
        elif opt == "paralogs":
            Settings.set_paralogs_folder(str(value))
        elif opt == "type":
            Settings.set_input_type(str(value))
        elif opt == "merge_method":
            Settings.set_method_merge(str(value))
        elif opt == "parameter":
            Settings.set_parameter_1(str(value))
        elif opt == "output":
            Settings.set_output_file(str(value))
        elif opt == "species_tree":
            Settings.set_input_tree(str(value))
        elif opt == "dynamic":
            Settings.set_dynamic_target_threshold(str(value))
        elif opt == "unmerged":
            Settings.set_unmerged_threshold(str(value))
        elif opt == "genome_info":
            Settings.set_genome_info(str(value))
        else:
            raise KeyError("Unknown argument '{}'. Please check manual or report".format(opt))


def run_gethogs():
    ############################
    start_time = time.time()
    ############################

    print("\n --- WARTHOGs:")
    print("\t- Start at " + str(time.strftime("%H:%M on %Y-%m-%d ")))
    print("\t- Orthology relations folder:" + str(Settings.pairwise_folder) +" ("+ str(Settings.input_type)+" format)")
    print("\t- Method use to merge HOGS: " + str(Settings.method_merge))
    print("\t- Output file: " + str(Settings.output_file))
    print("\n")

    Settings.check_required_argument()
    Settings.check_consistency_argument()
    Settings.set_xml_manager(file_manager.XML_manager())

    backbone_tree = Phylo.read(Settings.input_tree, "newick")

    if Settings.dynamic_treshold:
        lib.propagate_dynamic_threshold(backbone_tree.root, 0)
    lib.draw_tree(backbone_tree)
    lib.recursive_traversal(backbone_tree.root)

    Settings.xml_manager.finish_xml()

    # print StatisticTracker.frozen_per_level
    # print "cc_per_level", StatisticTracker.cc_per_level
    # print "hog_per_level",StatisticTracker.hog_per_level
    # print "levels",StatisticTracker.levels
    # print "notmerged_per_level",StatisticTracker.notmerged_per_level
    # print "nr_genes",StatisticTracker.nr_genes
    # print "nr_genes_per_genome",StatisticTracker.nr_genes_per_genome
    # print "merged_per_level",StatisticTracker.merged_per_level
    # print "time_per_level",StatisticTracker.time_per_level

    ############################
    end_time = time.time()
    print("--- %s seconds (total Time)" % (end_time - start_time))
    print("\n Done \\0/ \n")
    ############################


if __name__ == "__main__":
    Settings = Settings()
    StatisticTracker = StatisticTracker()
    init_configuration()
    run_gethogs()
