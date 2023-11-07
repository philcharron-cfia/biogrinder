from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from misc import *

class KmerCollection:
    """
    KmerCollection - A collection of kmers from sequences
    
    SYNOPSIS:
      col = KmerCollection(k=10, file='seqs.fa')
    
    DESCRIPTION:
    Manage a collection of kmers found in various sequences. Store information
    about what sequence a kmer was found in and its starting position on the
    sequence.
    
    APPENDIX:
    The rest of the documentation details each of the object methods.
    """
    
    def __init__(self, k=10, revcom=0, seqs=None, ids=None, file=None, weights=None):
        """
        Build a new kmer collection
        
        Args:
        - k:        set the kmer length (default: 10 bp)
        - revcom:   count kmers before and after reverse-complementing sequences
                    (default: 0)
        - seqs:     count kmers in the provided list of sequences
        - ids:      if specified, index the sequences provided to seqs using the
                    IDs in this list instead of using the sequences seq.id()
                    method
        - file:     count kmers in the provided file of sequences
        - weights:  if specified, assign the abundance of each sequence from the
                    values in this list
        
        Returns:
        - KmerCollection object
        """

        self._k = k if k is not None else 10
        self.revcom = revcom
        self.seqs = seqs if seqs else []
        self.ids = ids if ids else []
        self._weights = weights if weights else {}
        self.file = file
        self._collection_by_kmer = {}
        self._collection_by_seq = {}
        
        if seqs:
            self.add_seqs(seqs, ids)
        if file:
            self.add_file(file)

    @property
    def k(self):
        """
        Get the length of the kmers

        Returns:
            Positive integer
        """
        return self._k
    
    @k.setter
    def k(self, value):
        if not isinstance(value, int) or value <= 0:
            raise ValueError("k should be a positive integer.")
        self._k = value

    @property
    def weights(self):
        """
        Get or set countsthe weight of each sequence. Each sequence is given a
        weight of 1 by default.

        Returns:
            dict: A dictionary where the keys are sequence IDs and the values
                  are the weight of the corresponding (e.g. their relative abundance)
        """
        return self._weights
    
    @weights.setter
    def weights(self, value):
        if not isinstance(value, dict):
            raise ValueError("weights should be a dictionary.")
        self._weights = value

    @property
    def collection_by_kmer(self):
        """
        Get the collection of kmers, indexed by kmer.

        Returns:
            dict: A dictionary structure:
                  {kmer: {sequence_ID: [starts of kmer on sequence]}}
        """
        return self._collection_by_kmer

    @collection_by_kmer.setter
    def collection_by_kmer(self, value):
        if not isinstance(value, dict):
            raise ValueError(f"Expected a dictionary for collection_by_kmer, but got {type(value).__name__}")
        self._collection_by_kmer = value

    @property
    def collection_by_seq(self):
        """
        Get the collection of kmers, indexed by sequence ID.

        Returns:
            dict: A dictionary structure:
                  {sequence_ID: {kmer: [starts of kmer on sequence]}}
        """
        return self._collection_by_seq

    @collection_by_seq.setter
    def collection_by_seq(self, value):
        if not isinstance(value, dict):
            raise ValueError(f"Expected a dictionary for collection_by_seq, but got {type(value).__name__}")
        self._collection_by_seq = value

    def add_seqs(self, seqs, ids=None):
        """
        Process the kmers in the given sequences.

        Args:
            seqs (list): List of sequences to count kmers in.
            ids (list, optional): List of IDs to index the sequences.

        Returns:
            KmerCollection: Self for chaining.
        """
        col_by_kmer = self._collection_by_kmer
        col_by_seq = self._collection_by_seq
        for i, seq in enumerate(seqs):
            kmer_counts = self.find_kmers(seq)
            for kmer, positions in kmer_counts.items():
                if ids:
                    seq_id = ids[i]
                else:
                    seq_id = seq.id
                col_by_kmer.setdefault(kmer, {})[seq_id] = positions
                col_by_seq.setdefault(seq_id, {})[kmer] = positions
        return self

    def add_file(self, file):
        """
        Process the kmers in the given file of sequences.

        Args:
            file (str): Filename of sequences to count kmers in.

        Returns:
            KmerCollection: Self for chaining.
        """
        with open(file, "r") as f:
            for record in SeqIO.parse(f, "fasta"):
                self.add_seqs([record])
        return self
    
    def filter_rare(self, min_num):
        """
        Remove kmers occurring at less than the (weighted) abundance specified.

        Args:
            min_num (int): Minimum abundance for a kmer to be retained.

        Returns:
            KmerCollection: Self for chaining.
        """
        changed = False
        col_by_kmer = self._collection_by_kmer.copy()
        col_by_seq = self._collection_by_seq.copy()
        for kmer, sources in col_by_kmer.items():
            count = self.sum_from_sources(sources)
            if count < min_num:
                # Remove this kmer
                changed = True
                del col_by_kmer[kmer]
                for seq, seq_kmers in col_by_seq.items():
                    seq_kmers.pop(kmer, None)
                    if not seq_kmers:
                        del col_by_seq[seq]
        if changed:
            self._collection_by_kmer = col_by_kmer
            self._collection_by_seq = col_by_seq
        return self

    def filter_shared(self, min_shared):
        """
        Remove kmers occurring in less than the number of sequences specified.

        Args:
            min_shared (int): Minimum number of sequences a kmer must be in to be retained.

        Returns:
            KmerCollection: Self for chaining.
        """
        changed = False
        col_by_kmer = self._collection_by_kmer.copy()
        col_by_seq = self._collection_by_seq.copy()
        sequences_to_remove = []
        for kmer in list(col_by_kmer.keys()):
            sources = col_by_kmer[kmer]
            num_shared = len(sources)
            if num_shared < min_shared:
                # Remove this kmer
                changed = True
                del col_by_kmer[kmer]
                
                for seq in sources:
                    if seq in col_by_seq:
                        seq_kmers = col_by_seq[seq]
                        seq_kmers.pop(kmer, None)
                        if not seq_kmers:
                            sequences_to_remove.append(seq)
        for seq in sequences_to_remove:
            del col_by_seq[seq]
                #for seq, seq_kmers in col_by_seq.items():
                #    seq_kmers.pop(kmer, None)
                #    if not seq_kmers:
                #        del col_by_seq[seq]
                    
        #for kmer, sources in col_by_kmer.items():
        #    num_shared = len(sources)
        #    if num_shared < min_shared:
        #        # Remove this kmer
        #        changed = True
        #        del col_by_kmer[kmer]
        #        for seq, seq_kmers in col_by_seq.items():
        #            seq_kmers.pop(kmer, None)
        #            if not seq_kmers:
        #                del col_by_seq[seq]
        if changed:
            self._collection_by_kmer = col_by_kmer
            self._collection_by_seq = col_by_seq
        return self
    
    def counts(self, id=None, start=1, freq=0):
        """
        Calculate the total count of each kmer. Counts are affected by the
        weights given to the sequences.

        Args:
            id (str, optional): Restrict sequences to search to the specified sequence ID.
            start (int, optional): Starting position from which counting should start.
            freq (int): 0 to report counts (default), 1 to report frequencies (normalize to 1).

        Returns:
            tuple: 
                - list of the different kmers.
                - list of the corresponding total counts.
        """
        kmers = []
        counts = []
        total = 0

        for kmer, sources in self._collection_by_kmer.items():
            count = self.sum_from_sources(sources, id, start)
            if count > 0:
                kmers.append(kmer)
                counts.append(count)
                total += count
        
        if freq and total:
            counts = normalize(counts, total)
        
        return kmers, counts

    def sources(self, kmer, excl=None, freq=0):
        """
        Return the sources of a kmer and their (weighted) abundance.

        Args:
            kmer (str): kmer to get the sources of.
            excl (str, optional): Sources to exclude from the results.
            freq (int): 0 to report counts (default), 1 to report frequencies (normalize to 1).

        Returns:
            tuple:
                - list of the different sources.
                - list of the corresponding total counts.
            If the kmer requested does not exist, the list will be empty.
        """
        if not kmer:
            raise ValueError("Error: Need to provide a kmer to sources().")
        
        sources_list = []
        counts = []
        total = 0
        kmer_sources = self._collection_by_kmer.get(kmer, {})

        for source, positions in kmer_sources.items():
            if excl and source == excl:
                continue
            sources_list.append(source)
            weight = self._weights.get(source, 1) if self._weights else 1
            count = weight * len(positions)
            counts.append(count)
            total += count
        
        if freq and total:
            counts = normalize(counts, total)
        
        return sources_list, counts
    
    def kmers(self, seq_id, freq=0):
        """
        This is the inverse of sources(). Return the kmers found in a sequence
        (given its ID) and their (weighted) abundance.

        Args:
            seq_id (str): Sequence ID to get the kmers of.
            freq (int): 0 to report counts (default), 1 to report frequencies (normalize to 1).

        Returns:
            tuple: 
                - list of kmers.
                - list of the corresponding total counts.
            If the sequence ID requested does not exist, the list will be empty.
        """
        kmers_list = []
        counts = []
        total = 0
        seq_kmers = self._collection_by_seq.get(seq_id, {})

        for kmer, positions in seq_kmers.items():
            kmers_list.append(kmer)
            weight = self._weights.get(seq_id, 1) if self._weights else 1
            count = weight * len(positions)
            counts.append(count)
            total += count

        if freq:
            counts = normalize(counts, total)

        return kmers_list, counts

    def positions(self, kmer, source):
        """
        Return the positions of the given kmer on a given sequence. An error
        is reported if the kmer requested does not exist.

        Args:
            kmer (str): Desired kmer.
            source (str): Desired sequence with this kmer.

        Returns:
            list: List of the different positions. The list will be empty if the
            desired combination of kmer and sequence was not found.
        """
        kmer_sources = self._collection_by_kmer.get(kmer, {})
        return kmer_sources.get(source, [])
    
    def find_kmers(self, seq):
        """
        Find all kmers of size k in a sequence (either a Bio.Seq or Bio.SeqFeature)
        and return a dictionary where the keys are the kmers and the values are the
        positions of the kmers in the sequences.
        """
        k = self._k
        if isinstance(seq, Seq):
            seq_str = str(seq)
        elif isinstance(seq, SeqRecord):
            seq_str = str(seq.seq)
        else:
            seq_str = seq
        
        seq_str = seq_str.upper()  # case-insensitive
        seq_len = len(seq_str)
        kmer_positions = {}

        for i in range(seq_len - k + 1):
            kmer = seq_str[i:i+k]
            if kmer not in kmer_positions:
                kmer_positions[kmer] = []
            kmer_positions[kmer].append(i + 1)
        
        return kmer_positions

    def sum_from_sources(self, sources, id=None, start=1):
        """
        Calculate the number of (weighted) occurrences of a kmer. An optional
        sequence ID and start position to restrict the kmers can be specified.
        """
        count = 0
        if id is not None:
            sources = {id: sources.get(id, [])}

        for source, positions in sources.items():
            for position in positions:

                if position >= start:
                    weight = self._weights.get(source, 1) if self._weights else 1
                    count += weight
        
        return count


