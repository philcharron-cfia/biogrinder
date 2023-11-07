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


    def find_amplicons(self, sequence, primer_dict):
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
                            length = r[1] - f[1]
                            if length > 0:
                                barcode = sequence.id + "_" + primer_combo + "_" + str(i)
                                amplicon = Seq(sequence_str[f[1]+f[2]:r[1]])
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
                            length = f[1] - r[1]
                            if length > 0:
                                barcode = sequence.id + "_" + primer_combo + "_rev_" + str(i)
                                #amplicons[combo_name] = sequence_str[r[1]:f[1]+f[2]]
                                amplicon = Seq(sequence_str[r[1]+r[2]:f[1]])
                                amplicon_r = SeqRecord(amplicon,
                                                       id = barcode,
                                                       name = sequence.id)
                                amplicons.append(amplicon_r)
                                i += 1
        return amplicons



    