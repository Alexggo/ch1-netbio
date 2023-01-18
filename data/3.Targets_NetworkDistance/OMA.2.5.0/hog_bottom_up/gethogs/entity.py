from .settings import Settings
from .statistic_tracker import StatisticTracker
from . import genome_merger, lib


class Genome(object):
    IdCount = 0
    zoo = {}  # key: species_name & value: extent genomes
    retirement_house = {}  # key: taxon & value: ancestral genomes

    def __init__(self):
        self.HOGS = []  # list of HOGs related to this genome
        self.species = []  # list of all species contained within taxon
        self.children = []  # list of chidren obj:Genome
        self.genes = {}  # key: external id from pairwise files, value: Gene obj
        self.taxon = None  # taxon name
        self.type = None  # Ancestral or extent
        self.UniqueId = Genome.IdCount
        Genome.IdCount += 1

    @classmethod
    def get_extent_genomes(cls):
        if Settings.genome_info is not None:
            it = Settings.genome_info.keys()
        else:
            it = cls.zoo.keys()
        return [cls.zoo[species_name] for species_name in it
                if species_name in cls.zoo and cls.zoo[species_name].type == "extent"]

    # Actual genomes function
    def init_extent_genomes(self, node):
        """
        init the Genome obj with extent genomes specification
        :param node:
        :return:
        """
        self.species = [node.name]  # att:species only composed of the str:species_name
        print('-- Creation of ' +  str(self.species[0]) + " genome:")
        self.children = [self] # att:children only composed of the obj:Genome (trick to handle ancestral and extent genomes together)
        self.type = "extent"
        nr_genes = Settings.inputfile_handler.get_number_of_proteins(node.name)
        gene_range = range(1, nr_genes+1)
        self.create_genes_hogs_extent_genomes(gene_range)
        Genome.zoo[node.name] = self
        print('-> Genome of {} (composed of {:d} genes) created'.format(self.species[0], len(self.genes.keys())))
        StatisticTracker.add_extent_genome_stat(self.species[0], len(self.genes.keys()))

    def get_gene_by_ext_id(self, ext_id):
        try:
            return self.genes[ext_id]
        except KeyError as e:
            raise

    def create_genes_hogs_extent_genomes(self, gene_range):
        """
        for each genes of this species in pairwise folder, create the related obj:gene and obj:hog
        :param gene_range:
        :return:
        """

        for gene_ext_id in gene_range:
            gene = Gene(gene_ext_id, self)
            hog = HOG()
            hog.init_solohog(gene, self)
            gene.hog[self] = hog
            self.HOGS.append(hog)
            self.genes[gene.ext_id]=gene

    # Ancestral genomes function
    def init_ancestral_genomes(self, node):
        if Settings.dynamic_treshold:
            self.taxon_treshold = node.threshold
        self.children = [child.genome for child in node]
        for genome in self.children:
            for species in genome.species:
                self.species.append(species)
        self.set_taxon_name(node)
        StatisticTracker.add_taxanomic_range(self.taxon)
        StatisticTracker.frozen_per_level[self.taxon] = []
        print('-- Creation of the ancestral genome at ' +  str(self.taxon) + ":")
        merge = genome_merger.Merge_ancestral(self)
        self.HOGS = merge.newHOGs
        Genome.retirement_house[self.taxon] = self
        print('-> Genome at '+ str(self.taxon) + " (composed of " + str(len(self.HOGS)) +  " HOGs) created.") # to fill

    def set_taxon_name(self,node):
        if node.name:
            self.taxon = node.name.replace('__dc__',':').replace('__sc__',',')\
                         .replace('__po__', '(').replace('__pc__', ')')\
                         .replace('_', ' ')
        else:
            species_list = lib.get_list_species_name_genomes(self.children)
            taxon_name = ''
            for species in species_list:
                taxon_name = taxon_name + "/" + str(species)
            self.taxon = taxon_name


class HOG(object):
    IdCount = 0

    def __init__(self):
        self.genes = {} # key: obj:Genome and value: list obj:Gene
        self.xml = None # related orthologous groups in the xml
        self.topspecie = None # oldest obj:Genome where this hog is found (for taxon updating matters)
        self.id = HOG.IdCount
        self.unmerged_count = 0
        HOG.IdCount += 1

    def init_solohog(self, gene, species):
        self.topspecie = species
        self.genes[species] = [gene]
        Settings.xml_manager.create_xml_solohog(self)

    def merge_with(self, hog_to_merge_with):
        for key, value in hog_to_merge_with.genes.items():
            if key in self.genes:
                for g in value:
                    if g not in self.genes[key]:
                        self.genes[key].append(g)
            else:
                self.genes[key] = hog_to_merge_with.genes[key]

    def update_top_species_and_genes(self, new_top_genome):
        self.topspecie = new_top_genome
        for genome, genes in self.genes.items():
            for gene in genes:
                gene.hog[new_top_genome] = self


class Gene(object):
    IdCount = 0

    def __init__(self, ext_id, species):
        self.hog = {}  # key: obj:Genome  and value: obj:HOG
        self.ext_id = ext_id  # Id the paiwise files
        self.species = species  # obj:Genome that contains this gene
        if Settings.genome_info is not None:
            self.int_id = Settings.genome_info[species.species[0]].offset + ext_id
        else:
            self.int_id = Gene.IdCount  # unique id use by warthogs
            Gene.IdCount += 1

    def get_hog(self, query_genome):
        return self.hog[query_genome]





