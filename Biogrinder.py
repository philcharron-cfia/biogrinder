import arguments
from AmpliconSearch import AmpliconSearch
from KmerCollection import KmerCollection
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature
from Bio.Seq import Seq
from misc import *
import math
import random
import re
import sys

class Biogrinder:
    def __init__(self, *args):
        # Initialize the Biogrinder object with any necessary parameters
        self.args = args
        self.args = self.argparse()
        print(self.args)
        self.initialize() 

    def argparse(self):
        """
        Process arguments
        """
        parser = arguments.create_parser()
        args = parser.parse_args(self.args)    
        if args.profile_file:
            self.process_profile_file(args)
        else:
            # Set instance attributes from parsed arguments
            for arg, value in vars(args).items():
                setattr(self, arg, value)
        return args 
    
    def community_calculate_amplicon_abundance(self, r_spp_abs, r_spp_ids, seq_ids):
        '''
        Convert abundance of species into abundance of their amplicons because
        there can be multiple amplicon per species and the amplicons have a
        different ID from the species. The r_spp_ids and r_spp_abs lists are
        the ID and abundance of the species, sorted by decreasing abundance.
        '''
        i = 0
        while i < len(r_spp_ids):
            species_ab = r_spp_abs[i]
            species_id = r_spp_ids[i]
            amplicon_ids = list(seq_ids[species_id].keys())
            nof_amplicons = len(amplicon_ids)
            amplicon_abs = [species_ab / nof_amplicons] * nof_amplicons
            
            r_spp_abs = r_spp_abs[:i] + amplicon_abs + r_spp_abs[i+1:]
            r_spp_ids = r_spp_ids[:i] + amplicon_ids + r_spp_ids[i+1:]
            
            i += nof_amplicons

        return r_spp_abs, r_spp_ids
    
    def community_calculate_diversities(self, c_structs):
        overall_diversity = 0
        perc_shared = 0
        perc_permuted = 0

        # Calculate diversity (richness) based on given community abundances
        nof_libs = len(c_structs)
        all_ids = {}
        richnesses = []

        for c_struct in c_structs:
            richness = 0
            for i in range(len(c_struct['ids'])):
                id = c_struct['ids'][i]
                ab = c_struct['abs'][i]

                if not ab:
                    continue

                richness += 1

                if id in all_ids:
                    all_ids[id] += 1
                else:
                    all_ids[id] = 1
            
            richnesses.append(richness)

        overall_diversity = len(all_ids)

        # Calculate percent shared
        nof_non_shared = sum(1 for id, count in all_ids.items() if count < nof_libs)
        perc_shared = (overall_diversity - nof_non_shared) * 100 / overall_diversity

        return richnesses, overall_diversity, perc_shared, perc_permuted
    
    def community_calculate_species_abundance(self, distrib, param, diversity):
        '''
        Calculate relative abundance based on a distribution and its parameters.
        Input is a model, its 2 parameters, and the number of values to generate
        Output is a reference to a list of relative abundance. The abundance
        adds up to 1
        '''
        rel_ab = []
        total = 0

        if distrib == 'uniform':
            val = 1 / diversity
            for index in range(diversity):
                rel_ab.append(val)
            total = 1

        elif distrib == 'linear':
            slope = 1 / diversity
            for index in range(diversity):
                value = 1 - slope * index
                rel_ab.append(value)
                total += value

        elif distrib == 'powerlaw':
            if param is None:
                raise ValueError("Error: The powerlaw model requires an input parameter.")
            for index in range(diversity):
                value = (index + 1) ** -param
                rel_ab.append(value)
                total += value

        elif distrib == 'logarithmic':
            if param is None:
                raise ValueError("Error: The logarithmic model requires an input parameter.")
            for index in range(diversity):
                value = math.log(index + 2) ** -param
                rel_ab.append(value)
                total += value

        elif distrib == 'exponential':
            if param is None:
                raise ValueError("Error: The exponential model requires an input parameter.")
            for index in range(diversity):
                value = math.exp(- (index + 1) * param)
                rel_ab.append(value)
                total += value

        else:
            raise ValueError(f"Error: {distrib} is not a valid rank-abundance distribution.")

        # Normalize to 1 if needed
        if total != 1:
            rel_ab = [x/total for x in rel_ab]

        return rel_ab
    
    def community_given_abundances(self, file, seq_ids):
        """
        Read a file of genome abundances. The file should be space or
        tab-delimited. The first column should be the IDs of genomes, and the
        subsequent columns are for their relative abundance in different
        communities. An optional list of valid IDs can be provided. Then the
        abundances are normalized so that their sum is 1.
        """
        # Read abundances
        ids, abs = self.community_read_abundances(file)

        # Remove genomes with unknown IDs and calculate cumulative abundance
        totals = [0] * len(ids)
        
        for comm_num, id_list in enumerate(ids):
            i = 0
            while i < len(id_list):
                id = id_list[i]
                ab = int(abs[comm_num][i])
                if not seq_ids or id in seq_ids:
                    totals[comm_num] += ab
                    i += 1
                    
                else:
                    del id_list[i]
                    del abs[comm_num][i]
                    print(f"Requested reference sequence '{id}' in file "
                          f"'{file}' does not exist in the matching database.")
        
        # Process the communities
        c_structs = []
        for comm_num, (comm_ids, comm_abs) in enumerate(zip(ids, abs)):
            comm_total = totals[comm_num]
            if comm_total == 0:
                print("Warning: The abundance of all the genomes for community "
                      f"{comm_num + 1} was zero. Skipping this community.")
                continue

            # Normalize the abundances
            comm_abs = [int(ab) / comm_total for ab in comm_abs]

            # Sort relative abundances by decreasing 
            sorted_pairs = sorted(zip(comm_abs, comm_ids), reverse=True)
            comm_abs, comm_ids = zip(*sorted_pairs)
            
            # Save community structure
            c_structs.append({'ids': list(comm_ids), 'abs': list(comm_abs)})
        return c_structs
    
    def community_permuted(self, c_ids, perc_permuted):
        '''
        Change the abundance rank of species in all but the first community.
        The number of species changed in abundance is determined by the percent
        permuted, i.e. a given percentage of the most abundant species in this
        community.
        '''
        nof_indep = len(c_ids)
    
        # Leave the first community alone, but permute the ones after
        for c in range(2, nof_indep + 1):
            ids = c_ids[c-1]
            diversity = len(ids)
            
            # Number of top genomes to permute
            # Percent permuted is relative to diversity in this community
            nof_permuted = int((perc_permuted / 100) * diversity + 0.5) # round number

            idxs = []
            if nof_permuted > 0:
                # Add shuffled top genomes
                permuted_idxs = random.sample(range(nof_permuted), nof_permuted)
                idxs.extend(permuted_idxs)
            
            if diversity - nof_permuted > 0:
                # Add other genomes in the same order
                non_permuted_idxs = list(range(nof_permuted, diversity))
                idxs.extend(non_permuted_idxs)
            
            ids[:] = [ids[i] for i in idxs]

        return c_ids, perc_permuted

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

    def community_shared(self, seq_ids, nof_indep, perc_shared, diversities):
        '''
        Randomly split a library of sequences into a given number of groups
        that share a specified percent of their genomes.
        The % shared is the number of species shared / the total diversity in
        all communities
        Input:  arrayref of sequence ids
                number of communities to produce
                percentage of genomes shared between the communities
                diversity (optional, will use all genomes if not specified) 
        Return: arrayref of IDs that are shared
                arrayref of arrayref with the unique IDs for each community
        '''
    
        # If diversity is not specified (is '0'), use the maximum value possible
        nof_refs = len(seq_ids)
        min_diversity = float('inf')

        for i in range(len(diversities)):
            if diversities[i] == 0:
                diversities[i] = nof_refs / (perc_shared/100 + nof_indep*(1-perc_shared/100))
                diversities[i] = int(diversities[i])
                if i > 0 and diversities[i-1] != diversities[i]:
                    raise ValueError("Error: Define either all the diversities "
                                     "or none.")
            
            if diversities[i] < min_diversity:
                min_diversity = diversities[i]

        if min_diversity == 0:
            raise ValueError(f"Error: Cannot make {nof_indep} libraries sharing "
                             f"{perc_shared} % species from {nof_refs} references")

        nof_shared = int(min_diversity * perc_shared / 100)
        perc_shared = nof_shared * 100 / min_diversity

        nof_uniques = []
        sum_not_uniques = 0
        for diversity in diversities:
            nof_unique = diversity - nof_shared
            sum_not_uniques += nof_unique
            nof_uniques.append(nof_unique)

        overall_diversity = nof_shared + sum_not_uniques
        if nof_refs < overall_diversity:
            raise ValueError("Error: The number of reference sequences "
                             f"available ({nof_refs}) is not large enough to "
                             f"support the requested diversity ({overall_diversity} "
                             f"genomes overall with {perc_shared} % genomes "
                             f"shared between {nof_indep} libraries)")

        ids = list(seq_ids.keys())
        shared_ids = []
        for _ in range(nof_shared):
            rand_offset = random.randint(0, nof_refs-1)
            rand_id = ids.pop(rand_offset)
            nof_refs = len(ids)
            shared_ids.append(rand_id)

        unique_ids = [[] for _ in range(nof_indep)]
        for lib_num in range(nof_indep):
            nof_unique = nof_uniques[lib_num]
            for _ in range(nof_unique):
                rand_offset = random.randint(0, nof_refs-1)
                rand_id = ids.pop(rand_offset)
                nof_refs = len(ids)
                unique_ids[lib_num].append(rand_id)

        shared_ranks = random.sample(range(1, min_diversity+1), nof_shared)

        c_ranks = []
        for lib_num in range(nof_indep):
            diversity = diversities[lib_num]
            ranks = [None] * diversity
            
            for i in range(nof_shared):
                id_ = shared_ids[i]
                rank = shared_ranks[i]
                ranks[rank-1] = id_

            ids = unique_ids[lib_num]
            for rank in range(1, diversity+1):
                if ranks[rank-1] is None:
                    ranks[rank-1] = ids.pop()
            
            c_ranks.append(ranks)

        return c_ranks, overall_diversity, diversities, perc_shared
    
    def community_structures(self, seq_ids, abundance_file, distrib, param,
                             nof_indep, perc_shared, perc_permuted, diversities,
                             forward_reverse):
        """
        Create communities with a specified structure, alpha and beta-diversity.
        """

        # Calculate community structures
        c_structs = []
        if abundance_file:
            # Sanity check
            if len(diversities) > 1 or diversities[0]:
                print("Warning: Diversity cannot be specified when an "
                      "abundance file is specified. Ignoring it.")
            if perc_shared > 0 or perc_permuted < 100:
                print("Warning: Percent shared and percent permuted cannot be "
                      "specified when an abundance file is specified. Ignoring "
                      "them.")
            # One or several communities with specified rank-abundances
            c_structs = self.community_given_abundances(abundance_file, seq_ids)

            # Calculate number of libraries
            got_indep = len(c_structs)
            
            if nof_indep > got_indep:
                raise ValueError(f"Error: {nof_indep} communities were "
                                    "requested but the abundance file specified "
                                    f"the abundances for only {got_indep}.")
            elif nof_indep < got_indep:
                print(f"Warning: {nof_indep} communities were requested by "
                        "the abundance file but it specified the abundances "
                        f"for {got_indep}. Ignoring extraneous communities "
                        "specified in the file.")
            nof_indep = got_indep
            self.num_libraries = nof_indep

            # Calculate diversities based on given community abundances
            (self.diversity, self.overall_diversity, self.shared_perc,
             self.permuted_perc) = self.community_calculate_diversities(c_structs)
 
        else:
            # One or several communities with rank-abundance to be calculated
            # Sanity check
            if nof_indep == 1:  # 1 is the default value
                nof_indep = len(diversities)  
            
            if nof_indep != len(diversities):
                if len(diversities) == 1:
                    # Use same diversity for all libraries
                    diversity = diversities[0]
                    for _ in range(1, nof_indep):
                        diversities.append(diversity)
                else:
                    raise ValueError("Error: The number of richness values "
                                     f"provided ({len(diversities)}) did not "
                                     "match the requested number of libraries "
                                     f"({nof_indep}).")
            
            self.num_libraries = nof_indep
            # Select shared species
            c_ids = None
            overall_diversity = 0
            (c_ids, overall_diversity, diversities,
             perc_shared) = self.community_shared(seq_ids,
                                             nof_indep,
                                             perc_shared,
                                             diversities)
            
            # Shuffle the abundance-ranks of the most abundant genomes
            c_ids, perc_permuted = self.community_permuted(c_ids, perc_permuted)

            # Update values in self object 
            self.overall_diversity = overall_diversity
            self.diversity = diversities
            self.shared_perc = perc_shared
            self.permuted_perc = perc_permuted

            # Create a community structure list
            c_structs = []
            for c in range(1, nof_indep + 1):
                # Assign a random parameter if needed
                comm_param = param if param else randig(1, 0.05)
                # Calculate relative abundance of the community members
                diversity = self.diversity[c-1]
                c_abs = self.community_calculate_species_abundance(distrib, comm_param, diversity)
                c_id = c_ids[c-1]
                c_struct = {
                    'ids': c_id,
                    'abs': c_abs,
                    'param': comm_param,
                    'model': distrib
                }
                c_structs.append(c_struct)

        for c_struct in c_structs:
            # Convert sequence IDs to object IDs
            c_struct['abs'], c_struct['ids'] = self.community_calculate_amplicon_abundance(c_struct['abs'], c_struct['ids'], seq_ids)
        return c_structs

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
        seq_db = []     # sequence objects (all amplicons)
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
                    #print(amplicon_result)
                    for result in amplicon_result:
                        amp_seq = result.seq
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
                        #barcode = ref_seq.id + "_" + primer_id
                        result.seq = amp_seq
                        seq_db.append(result)
                        if ref_seq.id not in seq_ids:
                            seq_ids[ref_seq.id] = {}
                        seq_ids[ref_seq.id][result.id] = None
            else:
                seq_db.append(ref_seq)
                if ref_seq.id not in seq_ids:
                    seq_ids[ref_seq.id] = {}
                seq_ids[ref_seq.id][ref_seq.id] = None
                print(seq_ids)
            
        # Error if no usable sequences in the database
        if len(seq_ids) == 0:
            raise Exception("Error: No genome sequences could be used. If you " 
                            "specified a file of abundances for the genome "
                            "sequences, make sure that their ID match the ID "
                            "in the FASTA file. If you specified amplicon "
                            "primers, verify that they match some genome "
                            "sequences.")

        # Determine database type: dna, rna, protein
        db_alphabet = self.database_get_mol_type(mol_types)
        
        self.alphabet = db_alphabet
        # Error if using amplicon on protein database
        if db_alphabet == 'protein' and forward_reverse_primers is not None:
            raise ValueError("Error: Cannot use amplicon primers with proteic "
                             "reference sequences")

        # Error if using wrong direction on protein database
        if db_alphabet == 'protein' and unidirectional != 1:
            raise ValueError(f"Error: Got <unidirectional> = {unidirectional} "
                             "but can only use <unidirectional> = 1 with "
                             "proteic reference sequences")
        
        database = {'db': seq_db, 'ids': seq_ids}
        return database

    def database_get_children_seq(self, refseqid):
        """
        Retrieve all the sequences object made from a reference sequence based
        on the ID of the reference sequence.
        """
        children = []
        for child_oid in self.database['ids'][refseqid]:
            seq_obj = self.database_get_seq(child_oid)
            if seq_obj:
                children.append(seq_obj)
        return children

    def database_get_mol_type(self, mol_types):
        """
        Given a count of the different molecule types in the database, determine
        what molecule type it is.
        """
        # Determine the molecule type with the highest count
        max_type = max(mol_types, key=mol_types.get)
        max_count = mol_types[max_type]

        # Calculate the count of the other molecule types
        other_count = sum(count for type_, count in mol_types.items() if type_ != max_type)

        # Check if the max count is less than the count of the other molecule types
        if max_count < other_count:
            raise Exception("Error: Cannot determine what type of molecules "
                            f"the reference sequences are. Got {max_count} "
                            f"sequences of type '{max_type}' and {other_count} "
                            "others.")

        # Check if the molecule type is recognized
        if max_type not in ['dna', 'rna', 'protein']:
            raise Exception("Error: Reference sequences are in an unknown "
                            f"alphabet '{max_type}'")

        return max_type
    
    def database_get_parent_id(self, oid):
        """
        Based on a sequence object ID, retrieve the ID of the reference 
        sequence it came from
        """
        seq_id = self.database_get_seq(oid) 
        return seq_id.name

    def database_get_seq(self, oid):
        """
        Retrieve a sequence object from the database based on its object ID.
        """
        db = self.database["db"]
        seq_obj = next((record for record in db if record.id == oid), None)
        if seq_obj is None:
            print(f"Warning: Could not find sequence with object ID '{oid}' in the database")
        return seq_obj
    
    def initialize(self):
        # Parameter processing - read_dist  
        self.read_length = is_int(self.read_dist[0] if len(self.read_dist) > 0 else 100)
        self.read_model = is_option(self.read_dist[1] if len(self.read_dist) > 1 else 'uniform',
                                    ['uniform', 'normal'])
        self.read_delta = is_int(self.read_dist[2] if len(self.read_dist) > 2 else 0)

        # Parameter processing - insert_dist
        self.mate_length = is_int(self.insert_dist[0] if len(self.insert_dist) > 0 else 0)
        self.mate_model = is_option(self.insert_dist[1] if len(self.insert_dist) > 1 else 'uniform',
                                    ['uniform', 'normal'])
        self.mate_length = is_int(self.insert_dist[2] if len(self.insert_dist) > 2 else 0)

        # Parameter processing - abundance_model
        self.distrib = is_option(self.abundance_model[0] if len(self.abundance_model) > 0 else 'uniform',
                                    ['uniform','linear','powerlaw','logarithmic','exponential'])
        self.param = is_float(self.abundance_model[1] if len(self.abundance_model) > 1 else 1)

        # Parameter processing - mutation_dist
        self.mutation_model = is_option(self.mutation_dist[0] if len(self.mutation_dist) > 0 else 'uniform',
                                        ['uniform','linear','poly4'])
        self.mutation_para1 = is_int(self.mutation_dist[1] if len(self.mutation_dist) > 1 else 0)
        self.mutation_para2 = is_int(self.mutation_dist[2] if len(self.mutation_dist) > 2 else 0)
        

        # Parameter processing - mutation_ratio
        self.mutation_ratio.append(0) if len(self.mutation_ratio) == 1 else self.mutation_ratio[1]
        self.mutation_ratio_sum = self.mutation_ratio[0] + self.mutation_ratio[1]
        if self.mutation_ratio_sum == 0:
            self.mutation_ratio[0] = self.mutation_ratio[1] = 50
        else:
            self.mutation_ratio[0] = round(self.mutation_ratio[0] * 100 / self.mutation_ratio_sum, 1)
            self.mutation_ratio[1] = round(self.mutation_ratio[1] * 100 / self.mutation_ratio_sum, 1)
        
        # Parameter processing - chimera_dist
        self.chimera_dist_total = sum(self.chimera_dist)
        if self.chimera_dist_total == 0:
            self.chimera_dist = None
        else:
            self.chimera_dist = normalize(self.chimera_dist, self.chimera_dist_total)
            # Calculate cdf
            if self.chimera_perc:
                self.chimera_dist_cdf = self.proba_cumul(self.chimera_dist)

        # Parameter processing - fastq_output required qual_levels
        if self.fastq_output and (not self.qual_levels or len(self.qual_levels) == 0):
            raise ValueError("Error: <qual_levels> needs to be specified to output FASTQ reads")
        
        # Random number generator: seed or be auto-seeded
        if self.random_seed is not None:
            random.seed(self.random_seed)
        else:
            self.random_seed = random.seed()
        
        # Sequence length check
        self.max_read_length = self.read_length + self.read_delta  # approximation
        if self.mate_length:  # Check if mate_length is not zero
            self.min_mate_length = self.mate_length - self.mate_delta
            if self.max_read_length > self.min_mate_length:
                raise ValueError("Error: The mate insert length cannot be "
                                 "smaller than read length. Try increasing the "
                                 "mate insert length or decreasing the read "
                                 "length")
        
        # Pre-compile regular expression to check if reads are valid
        if self.exclude_chars is not None:
            self.exclude_re = re.compile(f"[{self.exclude_chars}]", re.IGNORECASE)  # Match any of the chars
        
        # Read MIDs
        if self.multiplex_ids is not None:
            self.multiplex_ids = self.read_multiplex_id_file(self.multiplex_ids, self.num_libraries)
        
        # Import reference sequences
        if self.chimera_dist_cdf:
            # Each chimera needs >= 1 bp. Use # sequences required by largest chimera.
            self.min_seq_len = len(self.chimera_dist) + 1
        else:
            self.min_seq_len = 1

        self.database = self.database_create(self.reference_file,
                                             self.unidirectional, 
                                             self.forward_reverse,
                                             self.abundance_file,
                                             self.delete_chars,
                                             self.min_seq_len,
                                             self.maximum_length)
        self.initialize_alphabet(self.alphabet)

        if (self.alphabet == 'protein' and 
            self.mate_length != 0 and 
            self.mate_orientation != 'FF'):
            raise Exception("Error: Can only use <mate_orientation> FF with "
                            "proteic reference sequences")
        
        # Genome relative abundance in the different independent libraries to create
        self.c_structs = self.community_structures(self.database['ids'],
                                                   self.abundance_file,
                                                   self.distrib,
                                                   self.param,
                                                   self.num_libraries,
                                                   self.shared_perc,
                                                   self.permuted_perc,
                                                   self.diversity,
                                                   self.forward_reverse)
        
        # Count kmers in the database if we need to form kmer-based chimeras
        if self.chimera_perc and self.chimera_kmer:
            # Get all wanted sequences (not all the sequences in the database)
            ids_dict = {}
            ids = []
            seqs = []

            for c_struct in self.c_structs:
                for id in c_struct["ids"]:
                    if id not in ids_dict:
                        ids_dict[id] = None
                        ids.append(id)
                        seqs.append(self.database_get_seq(id))                     
            ids_dict.clear()
            # Now create a collection of kmers
            self.chimera_kmer_col = KmerCollection(k=self.chimera_kmer,
                                                   seqs=seqs,
                                                   ids=ids).filter_shared(2)
        # Markers to keep track of computation progress
        self.cur_lib = 0
        self.cur_read = 0

    def initialize_alphabet(self, alphabet):
        """
        Store the characters of the alphabet to use and calculate their cdf so that
        we can easily pick them at random later.
        """
        self.alphabet_dict = {}
        
        if alphabet == 'dna':
            self.alphabet_dict = {
                'A': None,
                'C': None,
                'G': None,
                'T': None,
                'N': None,
                "-": None
            }
        elif alphabet == 'rna':
            self.alphabet_dict = {
                'A': None,
                'C': None,
                'G': None,
                'U': None,
                'N': None,
                "-": None
            }
        elif alphabet == 'protein':
            self.alphabet_dict = {
                'A': None, 'R': None, 'N': None, 'D': None, 'C': None,
                'Q': None, 'E': None, 'G': None, 'H': None, 'I': None,
                'L': None, 'K': None, 'M': None, 'F': None, 'P': None,
                'S': None, 'T': None, 'W': None, 'Y': None, 'V': None
                # 'B': None, # D or N
                # 'Z': None, # Q or E
                # 'X': None  # any amino-acid
                # J, O and U are the only unused letters
            }
        else:
            raise Exception(f"Error: unknown alphabet '{alphabet}'")
        
        num_chars = len(self.alphabet_dict)
        # CDF for this alphabet
        self.alphabet_complete_cdf = self.proba_cumul([1/num_chars] * num_chars)
        self.alphabet_truncated_cdf = self.proba_cumul([1/(num_chars-1)] * (num_chars-1))         

    def lib_coverage(self, c_struct):
        """
        Calculate number of sequences needed to reach a given coverage.
        If the number of sequences is provided, calculate the coverage.
        """
        coverage = self.coverage_fold
        nof_seqs = self.total_reads
        read_length = self.read_length

        # Calculate library length and size
        ref_ids = c_struct['ids']
        diversity = len(ref_ids)
        lib_length = 0

        for ref_id in ref_ids:
            seqobj = self.database_get_seq(ref_id)
            seqlen = len(seqobj) 
            lib_length += seqlen

        # Calculate number of sequences to generate based on desired coverage.
        # If both number of reads and coverage fold were given, coverage has
        # precedence.
        if coverage:
            nof_seqs = (coverage * lib_length) / read_length
            nof_seqs = int(nof_seqs) + (nof_seqs % 1 > 0)  # ceiling
        coverage = (nof_seqs * read_length) / lib_length

        # Sanity check
        if nof_seqs < diversity:
            print("Warning: The number of reads to produce is lower than the "
                  "required diversity. Increase the coverage or number of reads "
                  " to achieve this diversity.")
            self.diversity[self.cur_lib - 1] = nof_seqs
        return nof_seqs, coverage

    def next_lib(self):
        self.cur_lib += 1
        self.cur_read = 0
        self.cur_total_reads = 0
        self.cur_coverage_fold = 0
        self.next_mate = None
        self.positions = None
        
        if 0 <= self.cur_lib - 1 < len(self.c_structs):
            c_struct = self.c_structs[self.cur_lib - 1]
            # Create probabilities of picking genomes from community structure
            self.positions = self.proba_create(c_struct, self.length_bias, self.copy_bias)

            # Calculate needed number of sequences based on desired coverage
            self.cur_total_reads, self.cur_coverage_fold = self.lib_coverage(c_struct)
            
            # If chimeras are needed, update the kmer collection with sequence abundance
            kmer_col = self.chimera_kmer_col
            #print(weights)
            if kmer_col:
                weights = {}
                for i in range(len(c_struct['ids'])):
                    id = c_struct['ids'][i]
                    weight = c_struct['abs'][i]
                    weights[id] = weight
                kmer_col.weights = weights
                kmers, freqs = kmer_col.counts(None, 1, 1)
                self.chimera_kmer_arr = kmers
                self.chimera_kmer_cdf = self.proba_cumul(freqs)
        else:
            c_struct = None
        return c_struct            

    def proba_bias_dependency(self, c_struct, size_dep, copy_bias):
        '''
        Affect probability of picking a species by considering genome length
        or gene copy number bias
        '''
        # Calculate probability
        probas = []
        totproba = 0
        diversity = len(c_struct['ids'])
        for i in range(diversity):
            proba = c_struct['abs'][i]

            if self.forward_reverse:
                # Gene copy number bias
                if copy_bias:
                    refseq_id = self.database_get_parent_id(c_struct['ids'][i])
                    nof_amplicons = len(self.database_get_children_seq(refseq_id))
                    proba *= nof_amplicons
            else:
                # Genome length bias
                if size_dep:
                    id = c_struct['ids'][i]
                    seq = self.database_get_seq(id)
                    _len = len(seq)  # Assuming the seq is a string in Python
                    proba *= _len
                    

            probas.append(proba)
            totproba += proba

        # Normalize if necessary
        if totproba != 1:
            probas = normalize(probas, totproba)
        return probas

    def proba_create(self, c_struct, size_dep, copy_bias):
        # Calculate size-dependent, copy number-dependent probabilities
        probas = self.proba_bias_dependency(c_struct, size_dep, copy_bias)
        # Generate proba starting position
        positions = self.proba_cumul(probas)
        return positions
    
    def proba_cumul(self, probas):
        sum_val = 0
        cumul_probas = [0]
        for prob in probas:
            sum_val += prob
            cumul_probas.append(sum_val)
        return cumul_probas
    
    def process_profile_file(self, args):
        """
        Find profile file in arguments and read the profiles. The profile file
        only contains Biogrinder arguments, and lines starting with a '#' are
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

    def read_multiplex_id_file(self, file, nof_indep):
        self.mids = []
        # Read FASTA file containing the MIDs
        with open(file, 'r') as in_file:
            for record in SeqIO.parse(in_file, 'fasta'):
                self.mids.append(str(record.seq))
        # Sanity check
        self.nof_mids = len(self.mids)
        if self.nof_mids < nof_indep:
            raise ValueError(f"Error: {nof_indep} communities were requested "
                             f"but the MID file had only {self.nof_mids} sequences.")
        elif self.nof_mids > nof_indep:
            print(f"Warning: {nof_indep} communities were requested but the MID "
                  f"file contained {self.nof_mids} sequences. Ignoring extraneous MIDs.")      
        return self.mids