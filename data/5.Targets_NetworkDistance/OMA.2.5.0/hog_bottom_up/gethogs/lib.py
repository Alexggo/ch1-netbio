import time
from . import entity, statistic_tracker
from Bio import Phylo
from .settings import Settings


def propagate_dynamic_threshold(node, depth):
    children = [child for child in node]

    if not children:
        return [depth]
    else:
        list_descendant = []
        for child in node:
            child_depth = depth + 1
            list_child = propagate_dynamic_threshold(child, child_depth)
            list_descendant = list_descendant + list_child
        node.threshold = calculate_internal_score(list_descendant, depth)
        return list_descendant


def calculate_step(depth, leave_score, root_score):
    x = (int(leave_score) - int(root_score)) / (int(depth) - 1)
    return x


def calculate_internal_score(array, int_depth):
    score = 0
    for depth in array:
        step = calculate_step(depth, Settings.parameter_1, Settings.dynamic_treshold)
        score += int(Settings.parameter_1) - ((depth - 1 - int_depth) * step)
    return score / len(array)


def recursive_traversal(node):
    if node.is_terminal():
        start_time = time.time()
        node.genome = entity.Genome()
        node.genome.init_extent_genomes(node)
        end_time = time.time()
        print("<- %s seconds." % (end_time - start_time))
        print("\n")
    else:
        for child in node:
            recursive_traversal(child)
        start_time = time.time()
        node.genome = entity.Genome()
        node.genome.init_ancestral_genomes(node)
        end_time = time.time()
        statistic_tracker.StatisticTracker.set_time_per_level(node.genome.taxon, end_time - start_time)
        print("<- %s seconds." % (end_time - start_time))
        print("\n")


def draw_tree(tree):
    tree.ladderize()
    Phylo.draw_ascii(tree)


def get_list_species_name_genomes(list_genomes_object):
    """
    return the list of all species contained in the given genomes (all leaves in the subtree rooted by the genomes given)
    :param list_genomes_object:
    :return:
    """
    species_name = []
    for genome in list_genomes_object:
        for species in genome.species:
            species_name.append(str(species))
    return species_name


def get_percentage_orthologous_relations(hog_1, hog_2, extent_relations):
    """
    compute the % of extant orthologous relations / total possible relations between two hogs
    :param hog_1: 
    :param hog_2:
    :param extent_relations:
    :return:
    """
    n1 = sum(map(len, hog_1.genes.values()))
    n2 = sum(map(len, hog_2.genes.values()))
    #maximum_relations = len(hog_1.genes) * len(hog_2.genes)
    maximum_relations = n1 * n2
    score = float(extent_relations) / float(maximum_relations)
    return score * 100


def get_percentage_relations_between_two_set_of_HOGs(hogs_1, hogs_2, graph):
    """
    compute the % of extant orthologous relations / total possible relations between two set of hogs
    :param node1:
    :param node2:
    :param graph:
    :return:
    """

    number_pairwise_relations = get_all_pr_between_two_set_of_HOGs(graph, hogs_1, hogs_2)

    nbr_genes_hog1 = 0
    nbr_genes_hog2 = 0
    for hog1 in hogs_1:
        for g in hog1.genes.values():
            for e in g:
                nbr_genes_hog1 += 1
    for hog2 in hogs_2:
        for g in hog2.genes.values():
            for e in g:
                nbr_genes_hog2 += 1
    maximum_relations = nbr_genes_hog1 * nbr_genes_hog2
    score = float(number_pairwise_relations) / float(maximum_relations)
    score = score * 100
    return score


def get_all_pr_between_two_set_of_HOGs(graph, list_hogs1, list_hogs2):
    """
    return the numbers of extant orthologous relations between two set of hogs
    :param graph:
    :param list_hogs1:
    :param list_hogs2:
    :return:
    """
    pairwise_relations = 0
    for hog1 in list_hogs1:
        for hog2 in list_hogs2:
            try:
                pairwise_relations += graph[(hog1, hog2)]

            except KeyError:
                try:
                    pairwise_relations += graph[(hog2, hog1)]
                except KeyError:
                    pass

    return pairwise_relations
