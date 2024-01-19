from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import re

class SimulatedRead(SeqRecord):
    def __init__(self, id, reference, start, end, strand, mid, alphabet, track, coord_style, qual_levels):
        sequence = Seq(reference[start:end])
        if strand == -1:
            if alphabet == "dna":
                sequence = str(Seq(sequence).reverse_complement())
            elif alphabet == "rna":
                sequence = str(Seq(sequence).reverse_complement_rna())

        sequence = Seq(mid) + sequence
        super().__init__(sequence)
        self.id = id
        self.reference = reference
        self.start = start
        self.end = end
        self.strand = strand
        self.mid = mid
        self.track = track
        self.coord_style = coord_style
        self.qual_levels = qual_levels
        self.letter_annotations["phred_quality"] = self.generate_quality_scores(self.seq, self.qual_levels)
        
    def generate_quality_scores(self, seq, qual_levels, error_specs = None):
        """
        Generate quality scores for a sequence based on the provided quality levels.
        """
           
        if not qual_levels or len(qual_levels) != 2:
            raise ValueError("qual_levels should be a list of two quality scores: [good, bad]")

        good_qual, bad_qual = qual_levels
        qual = [good_qual] * len(seq)
        if error_specs is None:
            return qual
        for position, mutations in sorted(error_specs.items(), key=lambda x: x[0]):
            if 0 <= position < len(qual):  # Ensure position is within a valid range
                qual[position] = bad_qual
        return qual
    
    def errors(self, error_specs=None):
        """
        Get or set the sequencing errors and update the read.
        """
        if error_specs is None:
            return self._errors if hasattr(self, '_errors') else {}

        # Modify the sequence based on error specifications
        seq_list = list(self.seq)

        desc = self.description 
        desc = desc + " errors="
        self.description = desc

        for position, mutations in sorted(error_specs.items(), key=lambda x: x[0]):
            # 0-based position adjustment
            #position -= 1

            for mut_type, value in mutations.items():
                if isinstance(value, list):
                    for v in value:
                        self._apply_error(seq_list, position, mut_type, v)
                else:
                    self._apply_error(seq_list, position, mut_type, value)
        self.letter_annotations = {}
        self.seq = Seq("".join(seq_list))
        # Store the error specs
        self._errors = error_specs
        self.description = re.sub(r'(,$)','',self.description)
        self.letter_annotations["phred_quality"] = self.generate_quality_scores(self.seq, self.qual_levels, error_specs)

        return self._errors

    def _apply_error(self, seq_list, position, mut_type, value):
        """Apply a single error to the sequence list."""
        desc = self.description 
            
        if 0 <= position < len(seq_list):  # Ensure position is within valid range
            desc = desc + str(position+1) + mut_type + seq_list[position] + ","
            if mut_type == '%':  # Substitution
                seq_list[position] = value   
            elif mut_type == '+':  # Insertion
                seq_list.insert(position + 1, value)
            elif mut_type == '-':  # Deletion
                del seq_list[position]
            
        self.description = desc

def new_subseq(fragnum, seq_feat, unidirectional, orientation, start, end, mid, alphabet,
               mate_number=None, lib_number=None, tracking=None, qual_levels=None):
    # Adjust start and end if out of bounds
    start = max(0, start)
    end = min(len(seq_feat), end) + 1
    # Build the sequence ID
    name_sep = '_'
    mate_sep = '/'  # mate pair indicator, by convention
    newid = str(fragnum)
    if lib_number is not None:
        newid = str(lib_number) + name_sep + newid
    if mate_number is not None:
        newid += mate_sep + str(mate_number)

    # Create a new simulated read object
    newseq = SimulatedRead(
        id=newid,
        reference=seq_feat.seq,
        start=start,
        end=end,
        strand=orientation,
        mid=mid,
        alphabet=alphabet,
        track=tracking,
        coord_style='genbank',
        qual_levels=qual_levels
    )
    newseq.reference_id = seq_feat.id
    if tracking == 1:
        desc = gen_subseq_desc(newseq, seq_feat, newseq.strand, unidirectional)
        newseq.description = desc


        if hasattr(seq_feat, '_chimera'):
            amplicon_desc = gen_subseq_desc(newseq, seq_feat, newseq.strand, unidirectional)
            #desc = amplicon_desc
            #desc = desc.replace("reference=", "reference=" + amplicon_desc)
            newseq.description = amplicon_desc
    newseq.name = seq_feat.description

    return newseq

def gen_subseq_desc(newseq, seq_feat, strand, unidirectional):
    locations = []
    if hasattr(seq_feat, '_chimera'):
        locations = seq_feat._chimera
    else:
        locations = [(newseq.start, newseq.end)]

    loc_strings = []
    for location in locations:
        strand = strand or 1
        if strand == 1:
            loc_str = f"{location[0]}..{location[1]}"
        elif strand == -1:
            loc_str = f"complement({location[0]}..{location[1]})"
        else:
            raise ValueError(f"Error: Strand should be -1 or 1, but got '{location}'")
        loc_strings.append(loc_str)

    desc = "reference=" + newseq.reference_id + " position=(" + ",".join(loc_strings) + ")"
    if unidirectional == -1:
        newseq.strand *= -1
        if newseq.strand == -1:
            desc = desc.replace("position=(", "position=complement(")

    return desc
