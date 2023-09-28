import argparse
import textwrap

# Create an ArgumentParser object
def create_parser():
    parser = argparse.ArgumentParser(
        description='''Grinder: Simulate metagenomic and metatranscriptomic
        reads.''',
        formatter_class=argparse.RawTextHelpFormatter
    )

    # Required parameters
    required_params = parser.add_argument_group(
        'Required parameters'
    )
    required_params.add_argument(
        '-rf', '--reference_file', '-gf', '--genome_file', 
        type=str,
        default='-',
        help=textwrap.dedent('''\
            FASTA file that contains the input reference sequences (full
            genomes, 16S rRNA genes, transcripts, proteins...) or '-' to
            read them from the standard input. See the README file for
            examples of databases you can use and where to get them from.
            Default: %(default)s''')
    )

    # Basic parameters
    basic_params = parser.add_argument_group(
        'Basic parameters'
    )
    basic_params = basic_params.add_mutually_exclusive_group()
    basic_params.add_argument(
        '-tr', '--total_reads', 
        type=check_int_positive, 
        help=textwrap.dedent('''\
            Number of shotgun or amplicon reads to generate for each library.
            Do not specify this if you specify the fold coverage.''')
    )
    basic_params.add_argument(
        '-cf', '--coverage_fold', 
        type=check_float_positive, 
        help=textwrap.dedent('''\
            Desired fold coverage of the input reference sequences (the output
            FASTA length divided by the input FASTA length). Do not specify
            this if you specify the number of reads directly.''')
    ) 

    # Advanced shotgun and amplicon parameters
    shot_amp_params = parser.add_argument_group(
        '''Advanced shotgun and amplicon parameters'''
    )
    shot_amp_params.add_argument(
        '-rd', '--read_dist', 
        nargs='+',
        default=['100'], 
        help=textwrap.dedent('''\
            Desired shotgun or amplicon read length distribution specified as:
            average length, distribution ('uniform' or 'normal') and standard
            deviation. Only the first element is required.
            Examples:
                - All reads exactly 101 bp long: 101
                - Uniform read distribution around 100+-10 bp: 100 uniform 10
            Default: %(default)s''')
    )
    shot_amp_params.add_argument(
        '-id', '--insert_dist', 
        nargs='+',
        default=['0'], 
        help=textwrap.dedent('''\
            Create paired-end or mate-pair reads spanning the given insert
            length. Important: the insert is defined in the biological sense,
            i.e. its length includes the length of both reads and of the
            stretch of DNA between them.
            Options:
                - 0 : off
                - Insert size distribution in bp, in the same format as the 
                  read length distribution (a typical value is 2,500 bp for
                  mate pairs)
            Default: %(default)s''')
    )
    shot_amp_params.add_argument(
        '-mo', '--mate_orientation', 
        type=str,
        choices=['FF','FR','RF','RR'],
        default=['FR'], 
        help=textwrap.dedent('''\
            When generating paired-end or mate-pair reads, specify the 
            orientation of the reads (F: forward, R: reverse).
            Options:
                - FR: --> <--  e.g. Sanger, Illumina paired-end, IonTorrent mate-pair
                - FF: --> -->  e.g. 454
                - RF: <-- -->  e.g. Illumina mate-pair
                - RR: <-- <--
            Default: %(default)s''')
    )
    shot_amp_params.add_argument(
        '-ec', '--exclude_chars', 
        type=str,
        default=None, 
        help=textwrap.dedent('''\
            Do not create reads containing any of the specified characters
            (case insensitive). For example, use 'NX' to prevent reads with
            ambiguities (N or X). Grinder will error if it fails to find a
            suitable read (or pair of reads) after 10 attempts. Consider using 
            <delete_chars>, which may be more appropriate for your case.
            Default: %(default)s''')
    )
    shot_amp_params.add_argument(
        '-dc', '--delete_chars', 
        type=str,
        default=None, 
        help=textwrap.dedent('''\
            Remove the specified characters from the reference sequences
            (case-insensitive), e.g. '-~*' to remove gaps (- or ~) or 
            terminator (*). Removing these characters is done once, when
            reading the reference sequences, prior to taking reads. Hence
            it is more efficient than <exclude_chars>.
            Default: %(default)s''')
    )
    shot_amp_params.add_argument(
        '-fr', '--forward_reverse', 
        type=str,
        help=textwrap.dedent('''\
            Use DNA amplicon sequencing using a forward and reverse PCR primer
            sequence provided in a FASTA file. The reference sequences and 
            their reverse complement will be searched for PCR primer matches.
            The primer sequences should use the IUPAC convention for degenerate
            residues and the reference sequences that do not match the specified
            primers are excluded. If your reference sequences are full genomes,
            it is recommended to use <copy_bias> = 1 and <length_bias> = 0 to 
            generate amplicon reads.''')
    )
    shot_amp_params.add_argument(
        '-un', '--unidirectional', 
        type=int,
        choices=[0,1,-1],
        default=0,
        help=textwrap.dedent('''\
            Instead of producing reads bidirectionally, from the reference
            strand and its reverse complement, proceed unidirectionally, from
            one strand only (forward or reverse).
            Options:
                - 0 (off, i.e. bidirectional)
                - 1 (forward)
                - -1 (reverse)
            Use <unidirectional> = 1 for amplicon and strand-specific
            transcriptomic or proteomic datasets.
            Default: %(default)i''')
    )
    shot_amp_params.add_argument(
        '-lb', '--length_bias', 
        type=int,
        choices=[0,1],
        default=1,
        help=textwrap.dedent('''\
            In shotgun libraries, sample reference sequences proportionally to
            their length. For example, in simulated microbial datasets, this
            means that at the same relative abundance, larger genomes
            contribute more reads than smaller genomes (and all genomes have
            the same fold coverage). If your reference sequences are full
            genomes, it is recommended to use <length_bias> = 0 to generate
            amplicon reads.
            Options:
                - 0 (no)
                - 1 (yes)                 
            Default: %(default)i''')
    )
    shot_amp_params.add_argument(
        '-cb', '--copy_bias', 
        type=int,
        choices=[0,1],
        default=1,
        help=textwrap.dedent('''\
            In amplicon libraries where full genomes are used as input, sample
            species proportionally to the number of copies of the target gene: 
            at equal relative abundance, genomes that have multiple copies of
            the target gene contribute more amplicon reads than genomes that
            have a single copy. If your reference sequences are full genomes,
            it is recommended to use <copy_bias> = 1 to generate amplicon reads.
            Options:
                - 0 (no)
                - 1 (yes)  
            Default: %(default)i''')
    )

    # Aberrations and sequencing errors
    seq_error_params = parser.add_argument_group(
        '''Aberrations and sequencing errors'''
    )
    seq_error_params.add_argument(
        '-md', '--mutation_dist', 
        nargs='+',
        default=['uniform', 0, 0],
        help=textwrap.dedent('''\
            Introduce sequencing errors in the reads, under the form of
            mutations (substitutions, insertions and deletions) at positions
            that follow a specified distribution (with replacement):
            model (uniform, linear, poly4), model parameters.
            For example:
                - To simulate Sanger errors, use a linear model where the error 
                  rate is 1%% at the 5' end of reads and 2%% at the 3' end:
                  linear,1,2.
                - To model Illumina errors using the 4th degree polynome
                  3e-3 + 3.3e-8 * i^4 (Korbel et al 2009), use:
                  poly4,3e-3,3.3e-8.
            Use the <mutation_ratio> option to alter how many of these
            mutations are substitutions or indels. 
            Default: %(default)s''')
    )
    seq_error_params.add_argument(
        '-mr', '--mutation_ratio', 
        nargs='+',
        type=check_float_positive,
        default=[80, 20],
        help=textwrap.dedent('''\
            Indicate the percentage of substitutions and the number of indels
            (insertions and deletions). For example, use '80 20' (4
            substitutions for each indel) for Sanger reads. Note that this
            parameter has no effect unless you specify the <mutation_dist> 
            option. 
            Default: %(default)s''')
    )
    seq_error_params.add_argument(
        '-hd', '--homopolymer_dist', 
        choices=[0,'margulies','richter','balzer'],
        default=0,
        help=textwrap.dedent('''\
            Introduce sequencing errors in the reads under the form of
            homopolymeric stretches (e.g. AAA, CCCCC) using a specified model
            where the homopolymer length follows a normal distribution
            N(mean, standard deviation) that is function of the homopolymer
            length n.
            Options:
                - margulies: N(n, 0.15 * n)                 Margulies et al. 2005.
                - richter  : N(n, 0.15 * sqrt(n))           Richter et al. 2008.
                - balzer   : N(n, 0.03494 + n * 0.06856)    Balzer et al. 2010.
            Default: %(default)s''')
    )
    seq_error_params.add_argument(
        '-cp', '--chimera_perc', 
        type=check_float_percent,
        default=0,
        help=textwrap.dedent('''\
            Specify the percent of reads in amplicon libraries that should be
            chimeric sequences. The 'reference' field in the description of
            chimeric reads will contain the ID of all the reference sequences
            forming the chimeric template. A typical value is 10%% for
            amplicons. This option can be used to generate chimeric shotgun
            reads as well.
            Default: %(default).0f''')
    )
    seq_error_params.add_argument(
        '-cd', '--chimera_dist', 
        type=float,
        nargs='+',
        default=[314, 38, 1],
        help=textwrap.dedent('''\
            Specify the distribution of chimeras: bimeras, trimeras,
            quadrameras and multimeras of higher order. The default is the
            average values from Quince et al. 2011: '314 38 1', which
            corresponds to 89%% of bimeras, 11%% of trimeras and 0.3%% of
            quadrameras. Note that this option only takes effect when you
            request the generation of chimeras with the <chimera_perc> option.
            Default: %(default)s''')
    )
    seq_error_params.add_argument(
        '-ck', '--chimera_kmer', 
        type=check_chimera_kmer,
        default=10,
        help=textwrap.dedent('''\
            Activate a method to form chimeras by picking breakpoints at places 
            where k-mers are shared between sequences. <chimera_kmer> represents
            k, the length of the k-mers (in bp). The longer the kmer, the more
            similar the sequences have to be to be eligible to form chimeras.
            The more frequent a k-mer is in the pool of reference sequences
            (taking into account their relative abundance), the more often this
            k-mer will be chosen. Note that this option only takes effect when
            you request the generation of chimeras with the <chimera_perc>
            option.
            Default: %(default)i''')
    )

    # Community structure and diversity
    diversity_params = parser.add_argument_group(
        '''Community structure and diversity'''
    )
    diversity_params.add_argument(
        '-af', '--abundance_file', 
        type=str,
        help=textwrap.dedent('''\
            Specify the relative abundance of the reference sequences manually
            in an input file. Each line of the file should contain a sequence
            name and its relative abundance (%%), e.g. 'seqABC 82.1' or
            'seqABC 82.1 10.2' if you are specifying two different
            libraries.''')
    )
    diversity_params.add_argument(
        '-am', '--abundance_model', 
        nargs='+',
        default=['uniform', 1],
        help=textwrap.dedent('''\
            Relative abundance model for the input reference sequences:
            uniform, linear, powerlaw, logarithmic or exponential. The uniform
            and linear models do not require a parameter, but the other models
            take a parameter in the range [0, infinity). If this parameter is
            not specified, then it is randomly chosen.
            Options:
                - uniform distribution: uniform
                - powerlaw distribution with parameter 0.1: powerlaw 0.1
                - exponential distribution with automatically chosen parameter: exponential
            Default: %(default)s''')
    )
    diversity_params.add_argument(
        '-nl', '--num_libraries', 
        type=check_int_positive,
        default=1,
        help=textwrap.dedent('''\
            Number of independent libraries to create. Specify how diverse and
            similar they should be with <diversity>, <shared_perc> and
            <permuted_perc>. Assign them different MID tags with
            <multiplex_ids>.
            Default: %(default)i''')
    )
    diversity_params.add_argument(
        '-mi', '--multiplex_ids', 
        type=str,
        help=textwrap.dedent('''\
            Specify an optional FASTA file that contains multiplex sequence
            identifiers (a.k.a MIDs or barcodes) to add to the sequences (one
            sequence per library, in the order given). The MIDs are included in
            the length specified with the -read_dist option and can be altered
            by sequencing errors. See the MIDesigner or BarCrawl programs to
            generate MID sequences.''')
    )
    diversity_params.add_argument(
        '-di', '--diversity', 
        type=check_int_positive_with_0,
        nargs='+',
        default=['0'],
        help=textwrap.dedent('''\
            This option specifies alpha diversity, specifically the richness,
            i.e. number of reference sequences to take randomly and include in
            each library. Use 0 for the maximum richness possible (based on the
            number of reference sequences available). Provide one value to make
            all libraries have the same diversity, or one richness value per
            library otherwise.
            Default: %(default)i''')    
    )
    diversity_params.add_argument(
        '-sp', '--shared_perc', 
        type=check_float_percent,
        default=0,
        help=textwrap.dedent('''\
            This option controls an aspect of beta-diversity. When creating 
            multiple libraries, specify the percent of reference sequences they
            should have in common (relative to the diversity of the least
            diverse library). 
            Default: %(default).0f''')    
    )
    diversity_params.add_argument(
        '-pp', '--permuted_perc', 
        type=check_float_percent,
        default=100,
        help=textwrap.dedent('''\
            This option controls another aspect of beta-diversity. For multiple
            libraries, choose the percent of the most-abundant reference
            sequences to permute (randomly shuffle) the rank-abundance of.
            Default: %(default).0f''')    
    )

    # Miscellaneous options
    misc_params = parser.add_argument_group(
        '''Miscellaneous options'''
    )
    misc_params.add_argument(
        '-rs', '--random_seed', 
        type=check_int_positive,
        help=textwrap.dedent('''\
            Seed number to use for the pseudo-random number generator.''')    
    )
    misc_params.add_argument(
        '-dt', '--desc_track', 
        type=int,
        choices=[0,1],
        default=1,
        help=textwrap.dedent('''\
            Track read information (reference sequence, position, errors, ...)
            by writing it in the read description.
            Default: %(default)i''')    
    )
    misc_params.add_argument(
        '-ql', '--qual_levels', 
        type=check_int_positive_with_0,
        nargs='+',
        default=[None, None],
        help=textwrap.dedent('''\
            Generate basic quality scores for the simulated reads. Good
            residues are given a specified good score (e.g. 30) and residues
            that are the result of an insertion or substitution are given a
            specified bad score (e.g. 10). Specify first the good score and
            then the bad score on the command-line, e.g.: 30 10.''')    
    )
    misc_params.add_argument(
        '-fq', '--fastq_output', 
        type=int,
        choices=[0,1],
        default=0,
        help=textwrap.dedent('''\
            Whether to write the generated reads in FASTQ format (with
            Sanger-encoded quality scores) instead of FASTA and QUAL or not
            (1: yes, 0: no). <qual_levels> need to be specified for this option
            to be effective.
            Default: %(default)i''')    
    )
    misc_params.add_argument(
        '-bn', '--base_name', 
        type=str,
        default='grinder',
        help=textwrap.dedent('''\
            Prefix of the output files.
            Default: %(default)s''')    
    )
    misc_params.add_argument(
        '-od', '--output_dir', 
        type=str,
        default='.',
        help=textwrap.dedent('''\
            Directory where the results should be written. This folder will be
            created if needed.
            Default: %(default)s''')    
    )
    misc_params.add_argument(
        '-pf', '--profile_file', 
        type=str,
        help=textwrap.dedent('''\
            A file that contains Grinder arguments. This is useful if you use
            many options or often use the same options. Lines with comments (#)
            are ignored. Arguments must use the full length name. Consider the
            profile file, 'simple_profile.txt':
                # A simple Grinder profile
                read_dist 105 normal 12
                total_reads 1000
            Running: grinder -reference_file viral_genomes.fa -profile_file 
                simple_profile.txt
            Translates into: grinder -reference_file viral_genomes.fa
                --read_dist 105 normal 12 --total_reads 1000
            Note that the arguments specified in the profile should not be
            specified again on the command line.''')    
    )
    return parser

def check_int_positive(value):
    ivalue = int(value)
    if ivalue <= 0:
        raise argparse.ArgumentTypeError("%s is an invalid positive int value" % value)
    return ivalue

def check_float_positive(value):
    fvalue = float(value)
    if fvalue <= 0:
        raise argparse.ArgumentTypeError("%s is an invalid positive float value" % value)
    return fvalue

def check_float_percent(value):
    fvalue = float(value)
    if fvalue < 0 or fvalue > 100:
        raise argparse.ArgumentTypeError("%s is an invalid percent value" % value)
    return fvalue

def check_chimera_kmer(value):
    ivalue = int(value)
    if ivalue != 0 and ivalue < 2:
        raise argparse.ArgumentTypeError("%s must be 0 or more than 1" % value)
    return ivalue

def check_int_positive_with_0(value):
    ivalue = int(value)
    if ivalue < 0:
        raise argparse.ArgumentTypeError("%s is an invalid positive int value" % value)
    return ivalue
