
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
    

