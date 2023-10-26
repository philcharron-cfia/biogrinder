import os
import sys
from Bio import SeqIO
import Bio.SeqIO.FastaIO
import Bio.SeqIO.QualityIO
from Biogrinder import Biogrinder
from contextlib import nullcontext
import reports



def biogrinder(*args):
    # Create Biogrinder object
    factory = Biogrinder(*args)
    # Print diversity and percent shared and permuted
    reports.diversity_report(factory.num_libraries,
                             factory.shared_perc,
                             factory.permuted_perc,
                             factory.overall_diversity)
    # Create the output directory if needed
    if not os.path.isdir(factory.output_dir):
        try:
            os.mkdir(factory.output_dir)
        except OSError as e:
            raise OSError("Error: Could not create output folder "
                          f"{factory.output_dir}\n{e}")
        
    # Generate sequences
    while True:
        c_struct = factory.next_lib()
        if not c_struct:
            break
        cur_lib = factory.cur_lib
        
        # Output filenames
        if factory.num_libraries > 1:
            lib_str = '-' + f"{cur_lib:0{len(str(factory.num_libraries))}d}"
        else:
            lib_str = ''
        out_reads_basename = os.path.join(factory.output_dir, f"{factory.base_name}{lib_str}-reads.")

        if factory.fastq_output:
            out_fastq_file = out_reads_basename + 'fastq'
            out_fasta_file = None
            out_qual_file = None
        else:
            out_fasta_file = out_reads_basename + 'fa'
            if len(factory.qual_levels) > 0:
                out_qual_file = out_reads_basename + 'qual'
            out_fastq_file = None
        out_ranks_file = os.path.join(factory.output_dir, f"{factory.base_name}{lib_str}-ranks.txt")

        # Write community structure file
        factory.write_community_structure(c_struct, out_ranks_file)
        


        # Library report
        reports.library_report(cur_lib, factory.alphabet, factory.forward_reverse,
                    out_ranks_file, out_fastq_file, out_fasta_file, out_qual_file,
                    factory.cur_coverage_fold, factory.cur_total_reads,
                    factory.diversity[cur_lib-1])

        # Generate shotgun or amplicon reads and write them to a file
        with (open(out_fastq_file, 'a') if out_fastq_file else nullcontext()) as fastq_file, \
            (open(out_fasta_file, 'a') if out_fasta_file else nullcontext()) as fasta_file, \
            (open(out_qual_file, 'a') if out_qual_file else nullcontext()) as qual_file:


            while True:
                read = factory.next_read()
                if not read:
                    break
                if out_fastq_file:
                    SeqIO.write(read, fastq_file, "fastq")
                if out_fasta_file:
                    SeqIO.write(read, fasta_file, "fasta")
                if out_qual_file:
                    SeqIO.write(read, qual_file, "qual")

        if out_fastq_file: fastq_file.close()
        if out_fasta_file: fasta_file.close()
        if out_qual_file: qual_file.close()

def main():
    biogrinder(*sys.argv[1:])

if __name__ == '__main__':
    main()