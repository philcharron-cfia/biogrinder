from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Data import IUPACData
from pydna.common_sub_strings import common_sub_strings


class AmpliconSearch:
    def __init__(self):
        # Initialize the AmpliconSearch object with any necessary parameters
        self.primer_file = None        
        #primer_dict = self.create_primer_dict(primer_file)
        #primer_deg_dict = self.create_deg_primer_dict(primer_dict)
        #return primer_deg_dict

    
    def create_deg_primer_dict(self, primer_file):
        self.primer_file = primer_file
        self.primer_deg_dict = {}
        
        self.primer_dict = self.create_primer_dict(primer_file)
        self.primer_items = iter(self.primer_dict.items())
        for primer_name, primer_seq in self.primer_items:
            self.primer_deg_dict[primer_name] = self.expand_degenerate_primers(primer_seq)
        return self.primer_deg_dict

    def create_primer_dict(self, primer_file):
        primer_dict = {}
        with open(primer_file, 'r') as file:
            for record in SeqIO.parse(file, 'fasta'):
                primer_name = record.id
                primer_sequence = str(record.seq)
                primer_dict[primer_name] = primer_sequence
        file.close()
        if len(primer_dict) % 2 != 0:
            raise ValueError("A primer in the provided file does not have a "
                             "matching pair. Ensure each primer in the file has "
                             "a corresponding forward and reverse primer.")
        return primer_dict
    
    def degenerate_bases(self):
        """Return a dictionary of degenerate bases and their expansions."""
        # Filter out non-degenerate bases
        return {k: v for k, v in IUPACData.ambiguous_dna_values.items() if len(v) > 1}

    def expand_degenerate_primers(self, primer):
        """Expand a degenerate primer into all possible sequences."""
        bases = self.degenerate_bases()
        sequences = ['']

        for nucleotide in primer:
            if nucleotide in bases:
                sequences = [prefix + base for prefix in sequences for base in bases[nucleotide]]
            else:
                sequences = [prefix + nucleotide for prefix in sequences]
        return sequences

    def create_mismatch_primer_dict(self, primer_dict, mismatch):
        self.primer_dict = primer_dict
        self.primer_mismatch_dict = {}

        #self.primer_dict = self.create_primer_dict(primer_file)
        self.primer_items = iter(self.primer_dict.items())
        for primer_name, primer_seq_list in self.primer_items:
            self.primer_mismatch_dict[primer_name] = []
            for primer_seq in primer_seq_list:
                if mismatch == '1':
                    self.primer_mismatch_dict[primer_name].extend(self.expand_mismatch_primers_1(primer_seq))
                elif mismatch == '2':
                    self.primer_mismatch_dict[primer_name].extend(self.expand_mismatch_primers_2(primer_seq))
                elif mismatch == '3':
                    self.primer_mismatch_dict[primer_name].extend(self.expand_mismatch_primers_3(primer_seq))
        return self.primer_mismatch_dict

    def expand_mismatch_primers_1(self, primer):
        """Expand a primer into all possible sequences with up to 1 mismatch"""
        nucleotides = {'A', 'T', 'C', 'G'}
        sequences = set()
        sequences.add(primer)
        for i, nucleotide in enumerate(primer):
            for replacement in nucleotides - {nucleotide}:
                mutated_sequence = primer[:i] + replacement + primer[i+1:]
                sequences.add(mutated_sequence)
        
        return list(sequences)

    def expand_mismatch_primers_2(self, primer):
        """Expand a primer into all possible sequences with up to 2 mismatches"""
        nucleotides = {'A', 'T', 'C', 'G'}
        sequences = set()
        sequences.add(primer)

        # First mismatch
        for i, nucleotide1 in enumerate(primer):
            for replacement1 in nucleotides - {nucleotide1}:
                mutated_sequence1 = primer[:i] + replacement1 + primer[i+1:]

                # Add the sequence with one mismatch
                sequences.add(mutated_sequence1)

                # Second mismatch
                for j, nucleotide2 in enumerate(mutated_sequence1):
                    if j != i:  # Ensure the second mismatch is at a different position
                        for replacement2 in nucleotides - {nucleotide2}:
                            mutated_sequence2 = mutated_sequence1[:j] + replacement2 + mutated_sequence1[j+1:]
                            sequences.add(mutated_sequence2)

        return list(sequences)

    def expand_mismatch_primers_3(self, primer):
        """Expand a primer into all possible sequences with up to 3 mismatches"""
        nucleotides = {'A', 'T', 'C', 'G'}
        sequences = set()
        sequences.add(primer)

        # First mismatch
        for i, nucleotide1 in enumerate(primer):
            for replacement1 in nucleotides - {nucleotide1}:
                mutated_sequence1 = primer[:i] + replacement1 + primer[i+1:]

                # Add the sequence with one mismatch
                sequences.add(mutated_sequence1)

                # Second mismatch
                for j, nucleotide2 in enumerate(mutated_sequence1):
                    if j != i:  # Ensure the second mismatch is at a different position
                        for replacement2 in nucleotides - {nucleotide2}:
                            mutated_sequence2 = mutated_sequence1[:j] + replacement2 + mutated_sequence1[j+1:]

                            # Add the sequence with two mismatches
                            sequences.add(mutated_sequence2)

                            # Third mismatch
                            for k, nucleotide3 in enumerate(mutated_sequence2):
                                if k != i and k != j:  # Ensure the third mismatch is at a different position
                                    for replacement3 in nucleotides - {nucleotide3}:
                                        mutated_sequence3 = mutated_sequence2[:k] + replacement3 + mutated_sequence2[k+1:]
                                        sequences.add(mutated_sequence3)
        
        return list(sequences)

    def find_amplicons(self, sequence, primer_dict, maximum_length):
        """Find amplicons given potentially degenerate primers."""
        primer_items = iter(primer_dict.items())
        #sequence = sequence.upper() 
        amplicons = []
        for primer_name1, primer_seq1 in primer_items:
            primer_name2, primer_seq2 = next(primer_items)
            primer_combo = primer_name1 + "+" + primer_name2
            primer_seq1_rc = [str(Seq(seq).reverse_complement()) for seq in primer_seq1]
            primer_seq2_rc = [str(Seq(seq).reverse_complement()) for seq in primer_seq2]
            sequence_str = str(sequence.seq)
            match_F_array= [common_sub_strings(F, sequence_str, limit=len(F)) for F in primer_seq1]
            match_R_array= [common_sub_strings(R, sequence_str, limit=len(R)) for R in primer_seq2]
            match_F_rc_array= [common_sub_strings(F, sequence_str, limit=len(F)) for F in primer_seq1_rc]
            match_R_rc_array= [common_sub_strings(R, sequence_str, limit=len(R)) for R in primer_seq2_rc]          
            i=0
            for match_F in match_F_array:
                for f in match_F:
                    for match_R_rc in match_R_rc_array:
                        for r in match_R_rc:             
                            length = r[1] + r[2] - f[1]
                            if length > 0 and length < maximum_length:
                                barcode = sequence.id + "_" + primer_combo + "_LEN"  + str(length) + "_" + str(i)
                                amplicon = Seq(sequence_str[f[1]:r[1]+r[2]])
                                amplicon_r = SeqRecord(amplicon,
                                                       id = barcode,
                                                       name = sequence.id)
                                amplicons.append(amplicon_r)
                                #amplicons[combo_name] = sequence_str[f[1]:r[1]+r[2]]
                                i += 1
            i=0            
            for match_F_rc in match_F_rc_array:
                for f in match_F_rc:
                    for match_R in match_R_array:
                        for r in match_R:
                            length = f[1] + f[2] - r[1]
                            if length > 0 and length < maximum_length:
                                barcode = sequence.id + "_" + primer_combo + "_rev_LEN" + str(length) + "_" + str(i)
                                #amplicons[combo_name] = sequence_str[r[1]:f[1]+f[2]]
                                amplicon = Seq(sequence_str[r[1]:f[1]+f[2]])
                                amplicon_r = SeqRecord(amplicon,
                                                       id = barcode,
                                                       name = sequence.id)
                                amplicons.append(amplicon_r)
                                i += 1
        return amplicons
    
    def find_amplicons_v2(self, sequence, primer_dict, maximum_length):
        """Find amplicons given potentially degenerate primers."""
        sequence_str = str(sequence.seq)
        primer_lengths = []
        for primers in primer_dict.values():
            for seq in primers:
                if len(str(seq)) not in primer_lengths:
                    primer_lengths.append(len(seq))
        
        all_kmer_dicts = {}
        for i in primer_lengths:
            all_kmer_dicts[i] = self.build_seq_kmer_dictionary(sequence_str, i)
                    

        amplicons = []
        primer_items = iter(primer_dict.items())
        for primer_name1, primer_seq1 in primer_items:
            primer_name2, primer_seq2 = next(primer_items)
            primer_combo = primer_name1 + "+" + primer_name2
            primer_seq1_rc = [str(Seq(seq).reverse_complement()) for seq in primer_seq1]
            primer_seq2_rc = [str(Seq(seq).reverse_complement()) for seq in primer_seq2]
            match_F_array = [self.find_kmer_positions(F, all_kmer_dicts[len(F)]) for F in primer_seq1]
            match_R_array = [self.find_kmer_positions(R, all_kmer_dicts[len(R)]) for R in primer_seq2]
            match_F_rc_array = [self.find_kmer_positions(F, all_kmer_dicts[len(F)]) for F in primer_seq1_rc]
            match_R_rc_array = [self.find_kmer_positions(R, all_kmer_dicts[len(R)]) for R in primer_seq2_rc]
            match_F_array = self.clean_match_list(match_F_array)
            match_R_array = self.clean_match_list(match_R_array)
            match_F_rc_array = self.clean_match_list(match_F_rc_array)
            match_R_rc_array = self.clean_match_list(match_R_rc_array)
            i=0
            for match_F in match_F_array:
                primer_seq1_match = match_F[0]
                match_F = match_F[1]
                for f in match_F:
                    for match_R_rc in match_R_rc_array:
                        primer_seq2_match = match_R_rc[0]
                        match_R_rc = match_R_rc[1]
                        for r in match_R_rc:
                            length = r[0] + r[1] - f[0]
                            if length > 0 and length < maximum_length:
                                barcode = sequence.id + "_" + primer_combo + "_LEN"  + str(length) + "_" + str(i)
                                desc = sequence.id + "_" + primer_combo + \
                                        "_LEN" + str(length) + \
                                        "_START" + str(f[0]) + \
                                        "_END" + str(r[0]+r[1]) + \
                                        "_FPRIMER" + primer_seq1_match + \
                                        "_RPRIMER" + primer_seq2_match + \
                                        "_" + str(i)
                                amplicon = Seq(sequence_str[f[0]:r[0]+r[1]])
                                amplicon_r = SeqRecord(amplicon,
                                                       id = barcode,
                                                       description = desc,
                                                       name = sequence.id)
                                amplicons.append(amplicon_r)
                                i += 1
            i=0            
            for match_F_rc in match_F_rc_array:
                primer_seq1_match = match_F_rc[0]
                match_F_rc = match_F_rc[1]
                for f in match_F_rc:
                    for match_R in match_R_array:
                        primer_seq2_match = match_R[0]
                        match_R = match_R[1]
                        for r in match_R:
                            length = f[0] + f[1] - r[0]
                            if length > 0 and length < maximum_length:
                                barcode = sequence.id + "_" + primer_combo + "_rev_LEN" + str(length) + "_" + str(i)
                                desc = sequence.id + "_" + primer_combo + \
                                        "_rev_LEN" + str(length) + \
                                        "_START" + str(f[0]+f[1]) + \
                                        "_END" + str(r[0]) + \
                                        "_FPRIMER" + primer_seq1_match + \
                                        "_RPRIMER" + primer_seq2_match + \
                                        "_" + str(i)
                                amplicon = Seq(sequence_str[r[0]:f[0]+f[1]])
                                amplicon_r = SeqRecord(amplicon,
                                                       id = barcode,
                                                       description = desc,
                                                       name = sequence.id)
                                amplicons.append(amplicon_r)
                                i += 1
        return amplicons

    def clean_match_list(self, match_list):
        # Filter out sublists with empty nested lists
        non_empty_lists = [item for item in match_list if item[1]]
        # Initialize a list to store unique items
        unique_lists = []
        for item in non_empty_lists:
            # Convert the nested list to a tuple for comparison
            nested_tuple = tuple(map(tuple, item[1]))
            # Check if this tuple is not in any of the unique lists
            if not any(nested_tuple == tuple(map(tuple, existing_item[1])) for existing_item in unique_lists):
                unique_lists.append(item)
        return unique_lists

    def build_seq_kmer_dictionary(self, sequence, k):
        """Builds a dictionary of k-mers from the genome sequence."""
        kmer_dict = {}
        for i in range(len(sequence) - k + 1):
            kmer = sequence[i:i+k]
            if kmer in kmer_dict:
                kmer_dict[kmer].append(i)
            else:
                kmer_dict[kmer] = [i]
        return kmer_dict

    def find_kmer_positions(self, primer, kmer_dict):
        """Finds positions of a k-mer in the genome."""
        positions = kmer_dict.get(primer, [])
        len_primer = len(primer)
        position_len = []
        for position in positions:
            position_len.append([position, len_primer])
        return [primer, position_len]


    