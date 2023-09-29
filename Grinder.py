import arguments
from Bio import SeqIO
from misc import *
import random
import re

class Grinder:
    def __init__(self, *args):
        # Initialize the Grinder object with any necessary parameters
        self.args = args
        self.parsed_args = self.argparse()
        self.initialized_args = self.initialize()
        print(self.parsed_args)  

    def argparse(self):
        """
        Process arguments
        """
        parser = arguments.create_parser()
        args = parser.parse_args(self.args)    
        if args.profile_file:
            self.process_profile_file(args)
        return args

    def process_profile_file(self, args):
        """
        Find profile file in arguments and read the profiles. The profile file
        only contains Grinder arguments, and lines starting with a '#' are
        comments.
        """
        try:
            with open(args.profile_file, 'r') as file:
                for line in file.readlines():
                    line = line.strip()
                    if not line or line.startswith('#'):
                        continue # Skip empty lines and comments
                    arg_name, arg_value_str = line.split(' ', 1)
                    arg_values = arg_value_str.split()
                    if hasattr(args, arg_name):
                        setattr(args, arg_name, arg_values[0] if
                                len(arg_values) == 1 else arg_values)
                    else:
                        print(f"Warning: {arg_name} is not a recognized argument.")
        except FileNotFoundError:
            print(f"Error: Could not read file '{args.profile_file}'")   

    def initialize(self):
        args = self.parsed_args
        # Parameter processing - read_dist  
        self.split_list('read_dist', 'read_length', 0, 100, 'int', None)
        self.split_list('read_dist', 'read_model', 1, 'uniform', 'str', ['uniform', 'normal'])
        self.split_list('read_dist', 'read_delta', 2, 0, 'int', None)
        delattr(args, 'read_dist')

        # Parameter processing - insert_dist
        self.split_list('insert_dist', 'mate_length', 0, 0, 'int', None)
        self.split_list('insert_dist', 'mate_model', 1, 'uniform', 'str', ['uniform', 'normal'])
        self.split_list('insert_dist', 'mate_delta', 2, 0, 'int', None)
        delattr(args, 'insert_dist')

        # Parameter processing - abundance_model
        self.split_list('abundance_model', 'distrib', 0, 'uniform', 'str',
                        ['uniform','linear','powerlaw','logarithmic','exponential'])
        self.split_list('abundance_model', 'param', 1,  1, 'float', None)
        delattr(args, 'abundance_model')

        # Parameter processing - mutation_dist
        self.split_list('mutation_dist', 'mutation_model', 0, 0, 'str', ['uniform','linear','poly4'])
        self.split_list('mutation_dist', 'mutation_para1', 1, 0, 'int', None)
        self.split_list('mutation_dist', 'mutation_para2', 2, 0, 'int', None)
        delattr(args, 'mutation_dist')

        # Parameter processing - mutation_ratio
        mutation_ratio = getattr(args, 'mutation_ratio', [])
        mutation_ratio.append(0) if len(mutation_ratio) == 1 else mutation_ratio[1]
        mutation_ratio_sum = mutation_ratio[0] + mutation_ratio[1]
        if mutation_ratio_sum == 0:
            mutation_ratio[0] = mutation_ratio[1] = 50
        else:
            mutation_ratio[0] = round(mutation_ratio[0] * 100 / mutation_ratio_sum, 1)
            mutation_ratio[1] = round(mutation_ratio[1] * 100 / mutation_ratio_sum, 1)
        
        # Parameter processing - chimera_dist
        chimera_dist = getattr(args, 'chimera_dist', [])
        chimera_perc = getattr(args, 'chimera_perc', [])
        chimera_dist_total = sum(chimera_dist)
        if chimera_dist_total == 0:
            chimera_dist = None
        else:
            chimera_dist = normalize(chimera_dist, chimera_dist_total)
            # Calculate cdf
            if chimera_perc:
                chimera_dist_cdf = proba_cumul(chimera_dist)
                setattr(args, 'chimera_dist_cdf', chimera_dist_cdf) 
        setattr(args, 'chimera_dist', chimera_dist) 

        # Parameter processing - fastq_output required qual_levels
        fastq_output = getattr(args, 'fastq_output', [])
        qual_levels = getattr(args, 'qual_levels', [])
        if fastq_output and (not qual_levels or len(qual_levels) == 0):
            raise ValueError("Error: <qual_levels> needs to be specified to output FASTQ reads")
        
        # Random number generator: seed or be auto-seeded
        random_seed = getattr(args, 'random_seed', [])
        if random_seed is not None:
            random.seed(random_seed)
        else:
            random_seed = random.seed()
        
        # Sequence length check
        read_length = getattr(args, 'read_length', [])
        read_delta = getattr(args, 'read_delta', [])
        mate_length = getattr(args, 'mate_length', [])
        mate_delta = getattr(args, 'mate_delta', [])
        max_read_length = read_length + read_delta  # approximation
        if mate_length:  # Check if mate_length is not zero
            min_mate_length = mate_length - mate_delta
            if max_read_length > min_mate_length:
                raise ValueError("Error: The mate insert length cannot be "
                                 "smaller than read length. Try increasing the "
                                 "mate insert length or decreasing the read "
                                 "length")
        
        # Pre-compile regular expression to check if reads are valid
        exclude_chars = getattr(args, 'exclude_chars', [])
        if exclude_chars is not None:
            exclude_re = re.compile(f"[{exclude_chars}]", re.IGNORECASE)  # Match any of the chars
            setattr(args, 'exclude_re', exclude_re) 
        
        # Read MIDs
        multiplex_ids = getattr(args, 'multiplex_ids', [])
        num_libraries = getattr(args, 'num_libraries', [])
        if multiplex_ids is not None:
            multiplex_ids = self.read_multiplex_id_file(multiplex_ids, num_libraries)






        








    def split_list(self, arg, new_arg, index, def_value, type, options):
        args = self.parsed_args
        arg_value = getattr(args, arg, [])
        value = arg_value[index] if len(arg_value) > index else def_value
        if type == 'int':
            if is_int(value):
                setattr(args, new_arg, int(value)) 
        elif type == 'float':
            if is_float(value):
                setattr(args, new_arg, int(value)) 
        elif type == 'str':
            if is_option(value, options):
                setattr(args, new_arg, value)

    def read_multiplex_id_file(self, file, nof_indep):
        mids = []
        # Read FASTA file containing the MIDs
        with open(file, 'r') as in_file:
            for record in SeqIO.parse(in_file, 'fasta'):
                mids.append(str(record.seq))
        # Sanity check
        nof_mids = len(mids)
        if nof_mids < nof_indep:
            raise ValueError(f"Error: {nof_indep} communities were requested but the MID file had only {nof_mids} sequences.")
        elif nof_mids > nof_indep:
            print(f"Warning: {nof_indep} communities were requested but the MID file contained {nof_mids} sequences. Ignoring extraneous MIDs...")      
        return mids

        
    