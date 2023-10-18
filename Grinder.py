import arguments
from AmpliconSearch import AmpliconSearch
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature
from Bio.Seq import Seq
from misc import *
import random
import re
import sys

class Grinder:
    def __init__(self, *args):
        # Initialize the Grinder object with any necessary parameters
        self.args = args
        self.parsed_args = self.argparse()
        self.initialized_args = self.initialize()
        #print(self.parsed_args)  

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
        
        # Import reference sequences
        if hasattr(args, 'chimera_dist_cdf'):
            # Each chimera needs >= 1 bp. Use # sequences required by largest chimera.
            min_seq_len = len(chimera_dist) + 1
        else:
            min_seq_len = 1
        reference_file = getattr(args, 'reference_file', [])
        unidirectional = getattr(args, 'unidirectional', [])
        forward_reverse = getattr(args, 'forward_reverse', [])
        abundance_file = getattr(args, 'abundance_file', [])
        delete_chars = getattr(args, 'delete_chars', [])
        maximum_length = getattr(args, 'maximum_length', [])
        database = self.database_create(reference_file,
                                        unidirectional, 
                                        forward_reverse,
                                        abundance_file,
                                        delete_chars,
                                        min_seq_len,
                                        maximum_length)
            


        





        








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
    
    def database_create(self, fasta_file, unidirectional,
                        forward_reverse_primers=None, abundance_file=None,
                        delete_chars=None, min_len=1, maximum_length=10000):
        """
        Read and import sequences
        Parameters:
        * FASTA file containing the sequences or '-' for stdin. REQUIRED
        * Sequencing unidirectionally? 0: no, 1: yes forward, -1: yes reverse
        * Amplicon PCR primers (optional): Should be provided in a FASTA file
        and use the IUPAC convention. If a primer sequence is given, any
        sequence that does not contain the primer (or its reverse complement
        for the reverse primer) is skipped, while any sequence that matches is
        trimmed so that it is flush with the primer sequence
        * Abundance file (optional): To avoid registering sequences in the
        database unless they are needed
        * Delete chars (optional): Characters to delete form the sequences.
        * Minimum sequence size: Skip sequences smaller than that
        """
        # Input filehandle
        if not fasta_file:
            raise ValueError('Error: No reference sequences provided')
        if fasta_file == '-':
            in_handle = sys.stdin
        else:
            in_handle = open(fasta_file, 'r')
        
        # Get list of all IDs with a manually-specified abundance
        ids_to_keep = {}
        if abundance_file:
            ids, abs = self.community_read_abundances(abundance_file)
            for comm in ids:    
                for id in comm:
                    ids_to_keep[id] = None
        
        # Read FASTA file containing the primers   
        if forward_reverse_primers:
            amplicon_search = AmpliconSearch()
            primer_dict = amplicon_search.create_deg_primer_dict(primer_file=forward_reverse_primers)
        else:
            amplicon_search = 0

        # Process database sequences
        seq_db = {}     # sequence objects (all amplicons)
        seq_ids = {}    # reference sequence IDs and IDs of their amplicons
        mol_types = {}  # count of molecule types (dna, rna, protein)

        for ref_seq in SeqIO.parse(in_handle, "fasta"):
            # Skip empty sequences
            if not ref_seq.seq:
                continue
            # Record molecule type
            seq_type = identify_sequence_type(ref_seq.seq)
            mol_types[seq_type] = mol_types.get(seq_type, 0) + 1
            # Skip unwanted sequences
            if ids_to_keep and ref_seq.id not in ids_to_keep:
                continue
            # If we are sequencing from the reverse strand, reverse complement now
            if unidirectional == -1:
                ref_seq = SeqRecord(ref_seq.seq.reverse_complement(),
                                    id=ref_seq.id,
                                    name=ref_seq.name,
                                    description=ref_seq.description)
            # Extract amplicons if needed  
            if amplicon_search:
                amplicon_result = amplicon_search.find_amplicons(ref_seq, primer_dict)
                if len(amplicon_result) > 0:
                    for primer_id, amp_seq in amplicon_result.items():
                        # Remove forbidden chars
                        if delete_chars:
                            clean_seq = amp_seq
                            for char in delete_chars:
                                clean_seq = clean_seq.replace(char, '')
                            # Update sequence with cleaned sequence string
                            amp_seq = clean_seq
                        # Skip the sequence if it is too small
                        if len(amp_seq) < min_len:
                            continue
                        # Skip the sequence if it is too long
                        if len(amp_seq) > maximum_length:
                            continue
                        barcode = ref_seq.id + "_" + primer_id
                        seq_db[barcode] = amp_seq
                        if ref_seq.id not in seq_ids:
                            seq_ids[ref_seq.id] = {}
                        seq_ids[ref_seq.id][barcode] = None
            else:
                seq_db[ref_seq.name] = str(ref_seq.seq)
        # YOU ARE HERE AND MATCHING LINE IS 3063

    def community_read_abundances(self, file):
        # Read abundances of genomes from a file
        ids = []  # genome IDs
        abs = []  # genome relative abundance 
        with open(file, 'r') as io:
            for line_num, line in enumerate(io, 1):
                # Ignore comment or empty lines
                if line.strip() == '' or line.startswith('#'):
                    continue
                # Read abundance info from line
                parts = line.split()
                if parts:
                    id_ = parts[0]
                    rel_abs = parts[1:]
                    for comm_num, rel_ab in enumerate(rel_abs):
                        while len(ids) <= comm_num:
                            ids.append([])
                        while len(abs) <= comm_num:
                            abs.append([])
                        ids[comm_num].append(id_)
                        abs[comm_num].append(rel_ab)
                else:
                    print(f"Warning: Line {line_num} of file '{file}' has an \
                          unknown format. Skipping it...")

        return ids, abs