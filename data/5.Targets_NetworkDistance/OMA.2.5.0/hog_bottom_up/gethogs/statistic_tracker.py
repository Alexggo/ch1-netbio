__author__ = 'admin'


class StatisticTracker(object):

    '''
    StatisticTracker obj used in the warthogs to save various information
    '''

    # genome
    nr_genes = 0
    nr_genome = 0
    nr_genes_per_genome = {}
    zoo = []

    # level
    levels = []
    time_per_level = {}
    cc_per_level = {}
    merged_per_level = {}
    notmerged_per_level = {}
    hog_per_level = {}
    frozen_per_level = {}

    @classmethod
    def add_extent_genome_stat(cls, genome_name, nr_genes):
         cls.nr_genes += nr_genes
         cls.nr_genome += 1
         cls.zoo.append(genome_name)
         cls.nr_genes_per_genome[genome_name] = nr_genes

    @classmethod
    def add_taxanomic_range(cls, tr_name):
        cls.levels.append(tr_name)

    @classmethod
    def set_time_per_level(cls, level_name, time):
        cls.time_per_level[level_name] = time

    @classmethod
    def add_level_stat(cls, level_name, cc_per_level, merged_per_level, notmerged_per_level, hog_per_level):
        cls.levels.append(level_name)
        cls.cc_per_level[level_name] = cc_per_level
        cls.merged_per_level[level_name] = merged_per_level
        cls.notmerged_per_level[level_name] = notmerged_per_level
        cls.hog_per_level[level_name] = hog_per_level




