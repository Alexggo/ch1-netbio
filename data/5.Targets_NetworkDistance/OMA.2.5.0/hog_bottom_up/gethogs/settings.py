from Bio import Phylo
import sys
import collections
from gethogs import file_manager

GenomeInfo = collections.namedtuple('GenomeInfo', ['name', 'offset', 'nr_genes', 'taxid', 'sciname', 'dbname', 'release', 'date'])
GenomeInfo.__new__.__defaults__ = (0,) + ('', ) * (len(GenomeInfo._fields)-4)


class Settings(object):
    '''
    Settings obj used in the warthogs different obj and file
    '''

    pairwise_folder = None
    paralogs_folder = None
    input_type = None
    inputfile_handler = None
    method_merge = None
    input_tree = None
    parameter_1 = None
    output_file = None
    xml_manager = None
    list_species = None
    dynamic_treshold = False
    unmerged_treshold = False
    genome_info = None
    oma_id_format = False

    @classmethod
    def check_required_argument(cls):
        if cls.pairwise_folder == None or cls.input_type == None or cls.method_merge == None or cls.input_tree == None or cls.parameter_1 == None or cls.output_file == None:
            print(
            'There is a missing parameter, please check that you use all of the following arguments: -i -t -m -p -o -s')
            print('found: {0.pairwise_folder}, {0.input_type}, {0.method_merge}, {0.input_tree}, {0.parameter_1}, {0.output_file}'.format(cls))
            sys.exit(1)

    @classmethod
    def check_consistency_argument(cls):

        if cls.input_type not in ["oma", "standalone"]:
            print('The type of the pairwise folder you specified is not valid \n')
            print('Type provided:', cls.input_type)
            print('Types available: oma, standalone')
            sys.exit(1)

        if cls.method_merge not in ["update", "pair"]:
            print('The type of method to merge HOGs you specified is not valid \n')
            print('Type provided:', cls.method_merge)
            print('Types available: pair, update')
            sys.exit(1)

        cls.inputfile_handler = file_manager.inputfile_handler_factory()

        # Check consistency between species tree and pairwise files
        species_tree = Phylo.read(cls.input_tree, "newick")
        species_input_tree = [sp.name for sp in species_tree.get_terminals()]
        species_orthologs_file = cls.inputfile_handler.get_list_of_species()

        if not set(species_input_tree).issubset(set(species_orthologs_file)):
            print('Inconsistency between species in the species tree and the pairwise folder \n')
            print('Species in the tree:', species_input_tree, ' \n')
            print('Species in the pairwise folder:', species_orthologs_file, ' \n')
            intersection = set(species_input_tree).intersection(set(species_orthologs_file))
            print('Intersection ', intersection, ' \n')
            print('Unique species tree ', set(species_input_tree) - intersection, ' \n')
            print('Unique file ', set(species_orthologs_file) - intersection, ' \n')
            sys.exit(1)

        cls.list_species = species_input_tree

    @classmethod
    def set_pairwise_folder(cls, parameter):
        cls.pairwise_folder = parameter

    @classmethod
    def set_paralogs_folder(cls, parameter):
        cls.paralogs_folder = parameter

    @classmethod
    def set_input_type(cls, parameter):
        cls.input_type = parameter

    @classmethod
    def set_dynamic_target_threshold(cls, parameter):
        cls.dynamic_treshold = parameter

    @classmethod
    def set_unmerged_threshold(cls, parameter):
        cls.unmerged_treshold = parameter

    @classmethod
    def set_method_merge(cls, parameter):
        cls.method_merge = parameter

    @classmethod
    def set_parameter_1(cls, parameter):
        cls.parameter_1 = parameter

    @classmethod
    def set_output_file(cls, parameter):
        cls.output_file = parameter

    @classmethod
    def set_input_tree(cls, parameter):
        cls.input_tree = parameter

    @classmethod
    def set_xml_manager(cls, parameter):
        cls.xml_manager = parameter

    @classmethod
    def set_genome_info(cls, fname):
        """specify a order of the genomes in the final orthoxml and a predefined number of genes.

        The file needs to have the following format:
        Class <tab> nr_genomes <tab> total_nr_genes\n
        GenomeID <tab> offset\n
        GenomeID <tab> offset\n (i.e. sum of genes in all previous genomes)
        ...

        :param fname: path to the file containing this information."""
        info = collections.OrderedDict()
        with open(fname) as fh:
            class_, nr_genomes, tot_nr_genes = fh.readline().split('\t')
            nr_genomes, tot_nr_genes = int(nr_genomes), int(tot_nr_genes)
            g1, off1, cls_nr, taxid1, sciname1, dbname1, release1, date1 = fh.readline().strip('\n').split('\t')
            off1 = int(off1)
            for line in fh:
                try:
                    g, off, cls_nr, taxid, sciname, dbname, release, date = line.strip('\n').split('\t')
                except ValueError:
                    raise ValueError('Cannot split "{}" into 7 parts'.format(line))
                off = int(off)
                info[g1] = GenomeInfo(g1, off1, off - off1, int(taxid1), sciname1, dbname1, release1, date)
                g1, off1, taxid1, sciname1, dbname1, release1, date1 = g, off, taxid, sciname, dbname, release, date
            info[g1] = GenomeInfo(g1, off1, tot_nr_genes - off1, int(taxid1), sciname1, dbname1, release1, date1)
        cls.genome_info = info
        cls.oma_id_format = True
