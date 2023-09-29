import os
import sys
import Bio.SeqIO
import Bio.SeqIO.FastaIO
import Bio.SeqIO.QualityIO
from Grinder import Grinder




def grinder(*args):
    factory = Grinder(*args)

def main():
    grinder(*sys.argv[1:])

if __name__ == '__main__':
    main()