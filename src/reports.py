import os

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
                   fasta_file, qual_file, coverage, nof_seqs, diversity, outdir,
                   basename, lib_num):
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


    filename = os.path.join(outdir, f"{basename}{lib_num}-library_report.txt")
    with open(filename, 'w') as out_file:
        out_file.write(f"{lib_alphabet} {lib_type} library {cur_lib}:\n")
        out_file.write(f"   Community structure  = {ranks_file}\n")
        if fastq_file:
            out_file.write(f"   FASTQ file           = {fastq_file}\n")
        if fasta_file:
            out_file.write(f"   FASTA file           = {fasta_file}\n")
        if qual_file:
            out_file.write(f"   QUAL file            = {qual_file}\n")
        out_file.write(f"   Library coverage     = {coverage} x\n")
        out_file.write(f"   Number of reads      = {nof_seqs}\n")
        out_file.write(f"   Diversity (richness) = {diversity}\n")

def amplicon_report(cur_lib, alphabet, amplicon_dict, outdir, basename, lib_num):
    lib_alphabet = alphabet.upper()
    
    if "protein" in lib_alphabet.lower():
        lib_alphabet = lib_alphabet.replace("protein", "Proteic", 1)

    lib_type = 'amplicon'
    
    sorted_dict = {key: amplicon_dict[key] for key in sorted(amplicon_dict)}

    print(f"{lib_alphabet} {lib_type} library {cur_lib} - amplicon summary:")
    print("   Amplicon Name\tMatch Number\tStart\tEnd\tAmplicon Size (bp)\tPrimer 1\tPrimer 2\tNumber of Reads")
    for amp, value in sorted_dict.items():
        amp_name = amp.split('_LEN')[0]
        size = amp.split('_LEN')[1].split('_')[0]
        start = amp.split('_START')[1].split('_')[0]
        end = amp.split('_END')[1].split('_')[0]
        fprimer = amp.split('_FPRIMER')[1].split('_')[0]
        rprimer = amp.split('_RPRIMER')[1].split('_')[0]
        amp_number = amp.split('_RPRIMER')[1].split('_')[1]
        print(f"   {amp_name}\t{amp_number}\t{start}\t{end}\t{size}\t{fprimer}\t{rprimer}\t{value}")

    filename = os.path.join(outdir, f"{basename}{lib_num}-amplicon_report.txt")
    with open(filename, 'w') as out_file:
        out_file.write(f"{lib_alphabet} {lib_type} library {cur_lib} - amplicon summary:\n")
        out_file.write("   Amplicon Name\tMatch Number\tStart\tEnd\tAmplicon Size (bp)\tPrimer 1\tPrimer 2\tNumber of Reads\n")
        for amp, value in sorted_dict.items():
            amp_name = amp.split('_LEN')[0]
            size = amp.split('_LEN')[1].split('_')[0]
            start = amp.split('_START')[1].split('_')[0]
            end = amp.split('_END')[1].split('_')[0]
            fprimer = amp.split('_FPRIMER')[1].split('_')[0]
            rprimer = amp.split('_RPRIMER')[1].split('_')[0]
            amp_number = amp.split('_RPRIMER')[1].split('_')[1]
            out_file.write(f"   {amp_name}\t{amp_number}\t{start}\t{end}\t{size}\t{fprimer}\t{rprimer}\t{value}\n")
      