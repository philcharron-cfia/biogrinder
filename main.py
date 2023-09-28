import os
import sys
import Bio.SeqIO
import Bio.SeqIO.FastaIO
import Bio.SeqIO.QualityIO
from Grinder import Grinder




def grinder(*args):
    # Create Grinder object (assuming a Grinder class exists)
    factory = Grinder(*args)

def temp():
    # Print diversity and percent shared and permuted
    diversity_report(
        factory.num_libraries,
        factory.shared_perc,
        factory.permuted_perc,
        factory.overall_diversity,
    )

    # Create the output directory if needed
    if not os.path.isdir(factory.output_dir):
        os.makedirs(factory.output_dir, exist_ok=True)

    # Generate sequences
    while True:
        c_struct = factory.next_lib()
        if not c_struct:
            break

        cur_lib = factory.cur_lib

        # Output filenames
        lib_str = ''
        if factory.num_libraries > 1:
            lib_str = f'-{cur_lib:0{len(str(factory.num_libraries))}d}'
        out_reads_basename = os.path.join(
            factory.output_dir, f'{factory.base_name}{lib_str}-reads.'
        )
        out_fasta_file = None
        out_qual_file = None
        out_fastq_file = None
        if factory.fastq_output:
            out_fastq_file = f'{out_reads_basename}fastq'
        else:
            out_fasta_file = f'{out_reads_basename}fa'
            if factory.qual_levels:
                out_qual_file = f'{out_reads_basename}qual'
        out_ranks_file = os.path.join(
            factory.output_dir, f'{factory.base_name}{lib_str}-ranks.txt'
        )

        # Write community structure file
        factory.write_community_structure(c_struct, out_ranks_file)

        # Prepare output FASTA file
        out_fastq = None
        if out_fastq_file is not None:
            out_fastq = Bio.SeqIO.FastaIO.FastaWriter(
                out_fastq_file, 'w', variant='sanger'
            )

        out_fasta = None
        if out_fasta_file is not None:
            out_fasta = Bio.SeqIO.FastaIO.FastaWriter(out_fasta_file, 'w')

        out_qual = None
        if out_qual_file is not None:
            out_qual = Bio.SeqIO.QualIO.QualWriter(out_qual_file)

        # Library report
        diversity = factory.diversity[cur_lib - 1]
        library_report(
            cur_lib,
            factory.alphabet,
            factory.forward_reverse,
            out_ranks_file,
            out_fastq_file,
            out_fasta_file,
            out_qual_file,
            factory.cur_coverage_fold,
            factory.cur_total_reads,
            diversity,
        )

        # Generate shotgun or amplicon reads and write them to a file
        while True:
            read = factory.next_read()
            if not read:
                break

            if out_fastq is not None:
                out_fastq.write_seq(read)
            if out_fasta is not None:
                out_fasta.write_seq(read)
            if out_qual is not None:
                out_qual.write_seq(read)

        if out_fastq is not None:
            out_fastq.close()
        if out_fasta is not None:
            out_fasta.close()
        if out_qual is not None:
            out_qual.close()

    return 1

def diversity_report(num_libraries, perc_shared, perc_permuted, overall_diversity):
    format_str = '{:.1f}'
    print(f"Overall diversity = {overall_diversity} genomes")

    if num_libraries > 1:
        nof_shared = perc_shared / 100 * overall_diversity
        perc_shared = format_str.format(perc_shared)
        print(f"Percent shared   = {perc_shared} % ({nof_shared} genomes)")
        nof_permuted = perc_permuted / 100 * overall_diversity
        perc_permuted = format_str.format(perc_permuted)
        print(f"Percent permuted = {perc_permuted} % ({nof_permuted} top genomes)")

    return 1

def main():
    grinder(*sys.argv[1:])

if __name__ == '__main__':
    main()
