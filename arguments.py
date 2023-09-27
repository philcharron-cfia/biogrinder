import argparse

# Create an ArgumentParser object
def create_parser():
    parser = argparse.ArgumentParser(description='''Grinder: Simulate metagenomic 
                                    and metatranscriptomic reads.''')

    # Required parameters
    required_params = parser.add_argument_group('Required parameters')
    required_params.add_argument('-rf', '--reference_file', 
                                '-gf', '--genome_file', 
                                type=str,
                                default='-',
                                help='''FASTA file that contains the input
                                reference sequences (full genomes, 16S rRNA genes,
                                transcripts, proteins...) or '-' to read them from
                                the standard input. See the README file for 
                                examples of databases you can use and where to get
                                them from. Default: %(default)s''')

    # Basic parameters
    basic_params = parser.add_argument_group('Basic parameters')
    mut_exc_basic_params = basic_params.add_mutually_exclusive_group()
    mut_exc_basic_params.add_argument('-tr', '--total_reads', 
                              type=int, 
                              help='''Number of shotgun or amplicon reads to 
                              generate for each library. Do not specify this if
                              you specify the fold coverage.''')
    mut_exc_basic_params.add_argument('-cf', '--coverage_fold', 
                              type=float, 
                              help='''Desired fold coverage of the input
                              reference sequences (the output FASTA length 
                              divided by the input FASTA length). Do not
                              specify this if you specify the number of reads 
                              directly.''') 
    
    # Advanced shotgun and amplicon parameters
    shotgun_amplicon_params = parser.add_argument_group('Advanced shotgun and amplicon parameters')

















    return parser

