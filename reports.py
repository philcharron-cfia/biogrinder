
def diversity_report(num_libraries, perc_shared, perc_permuted, overall_diversity):
    """
    Print a diversity report detailing the overall genome diversity, the
    percentage of shared genomes, and the percentage of permuted top genomes.
    
    Args:
    - num_libraries (int): Number of libraries.
    - perc_shared (float): Percentage of shared genomes.
    - perc_permuted (float): Percentage of permuted top genomes.
    - overall_diversity (float): Overall genome diversity.
    
    Returns:
    - None
    
    Note:
    If the number of libraries is greater than 1, additional calculations 
    and print statements for shared and permuted percentages are executed.
    """
    print(f"Overall diversity = {overall_diversity} contigs")
    
    if num_libraries > 1:
        nof_shared = int(perc_shared / 100 * overall_diversity)
        print(f"Percent shared   = {perc_shared:.1f} % ({nof_shared} genomes)")
        
        nof_permuted = int(perc_permuted / 100 * overall_diversity)
        print(f"Percent permuted = {perc_permuted:.1f} % ({nof_permuted} top genomes)")
    
def library_report(cur_lib, alphabet, forward_reverse, ranks_file, fastq_file,
                   fasta_file, qual_file, coverage, nof_seqs, diversity):
    coverage = "{:.1f}".format(coverage)
    lib_alphabet = alphabet.upper()
    
    if "protein" in lib_alphabet.lower():
        lib_alphabet = lib_alphabet.replace("protein", "Proteic", 1)

    lib_type = 'amplicon' if forward_reverse else 'shotgun'
    
    print(f"{lib_alphabet} {lib_type} library {cur_lib}:")
    print(f"   Community structure  = {ranks_file}")
    if fastq_file:
        print(f"   FASTQ file           = {fastq_file}")
    if fasta_file:
        print(f"   FASTA file           = {fasta_file}")
    if qual_file:
        print(f"   QUAL file            = {qual_file}")
    print(f"   Library coverage     = {coverage} x")
    print(f"   Number of reads      = {nof_seqs}")
    print(f"   Diversity (richness) = {diversity}")
