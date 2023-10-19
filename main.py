import os
import sys
import Bio.SeqIO
import Bio.SeqIO.FastaIO
import Bio.SeqIO.QualityIO
from Biogrinder import Biogrinder
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
        

def main():
    biogrinder(*sys.argv[1:])

if __name__ == '__main__':
    main()