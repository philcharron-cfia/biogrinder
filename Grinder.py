import arguments
from misc import *

class Grinder:
    def __init__(self, *args):
        # Initialize the Grinder object with any necessary parameters
        self.args = args
        self.parsed_args = self.argparse()
        self.initialized_args = self.initialize()
        print(self.parsed_args)  

    def argparse(self):
        """
        Process arguments
        """
        parser = arguments.create_parser()
        args = parser.parse_args(self.args)    
        if args.profile_file:
            self.process_profile_file(args)
        return args

    def process_profile_file(self, args):
        """
        Find profile file in arguments and read the profiles. The profile file
        only contains Grinder arguments, and lines starting with a '#' are
        comments.
        """
        try:
            with open(args.profile_file, 'r') as file:
                for line in file.readlines():
                    line = line.strip()
                    if not line or line.startswith('#'):
                        continue # Skip empty lines and comments
                    arg_name, arg_value_str = line.split(' ', 1)
                    arg_values = arg_value_str.split()
                    if hasattr(args, arg_name):
                        setattr(args, arg_name, arg_values[0] if
                                len(arg_values) == 1 else arg_values)
                    else:
                        print(f"Warning: {arg_name} is not a recognized argument.")
        except FileNotFoundError:
            print(f"Error: Could not read file '{args.profile_file}'")   

    def initialize(self):
        # Parameter processing: read length distribution
        '''
        args = self.parsed_args
        read_length = args.read_dist[0] if args.read_dist[0] else 100
        setattr(args, 'read_length', int(read_length)) if is_int(read_length) else None
        read_model = args.read_dist[1] if len(args.read_dist) > 1 else 'uniform'
        setattr(args, 'read_model', read_model) if is_option(read_model, ['uniform', 'normal']) else None
        read_delta = args.read_dist[2] if len(args.read_dist) > 2 else 0
        setattr(args, 'read_delta', int(read_delta)) if is_int(read_delta) else None
        delattr(args, 'read_dist')
        '''
        
        args = self.split_distr('read_dist',
                         ['read_length', 'read_model', 'read_delta'],
                         [100, 'uniform', 0])
        args = self.split_distr('insert_dist',
                         ['mate_length', 'mate_model', 'mate_delta'],
                         [0, 'uniform', 0])        
        return args
    

    def split_distr(self, arg, new_args, values):
        args = self.parsed_args
        arg_value = getattr(args, arg, [])
        length = arg_value[0] if len(arg_value) > 0 else values[0]
        if is_int(length):
            setattr(args, new_args[0], int(length)) 
        model = arg_value[1] if len(arg_value) > 1 else values[1]
        if is_option(model, ['uniform', 'normal']):
            setattr(args, new_args[1], model)
        delta = arg_value[2] if len(arg_value) > 2 else values[2]
        if is_int(delta):
            setattr(args, new_args[2], int(delta)) if is_int(delta) else None
        delattr(args, arg)

        return args
        























    def next_lib(self):
        # Implement the logic to get the next library's community structure data
        # and return it
        pass

    def next_read(self):
        # Implement the logic to generate and return the next sequence read
        pass

    def write_community_structure(self, c_struct, out_ranks_file):
        # Implement the logic to write community structure data to a file
        pass

    @property
    def num_libraries(self):
        # Implement a property to get the number of libraries
        pass

    @property
    def shared_perc(self):
        # Implement a property to get the shared percentage
        pass

    @property
    def permuted_perc(self):
        # Implement a property to get the permuted percentage
        pass

    @property
    def overall_diversity(self):
        # Implement a property to get the overall diversity
        pass

    @property
    def output_dir(self):
        # Implement a property to get the output directory
        pass

    @property
    def base_name(self):
        # Implement a property to get the base name
        pass

    @property
    def fastq_output(self):
        # Implement a property to check if fastq output is enabled
        pass

    @property
    def alphabet(self):
        # Implement a property to get the alphabet
        pass

    @property
    def forward_reverse(self):
        # Implement a property to get the forward/reverse information
        pass

    @property
    def qual_levels(self):
        # Implement a property to get quality levels
        pass

    @property
    def cur_lib(self):
        # Implement a property to get the current library
        pass

    @property
    def diversity(self):
        # Implement a property to get diversity data
        pass

    @property
    def cur_coverage_fold(self):
        # Implement a property to get the current coverage fold
        pass

    @property
    def cur_total_reads(self):
        # Implement a property to get the current total reads
        pass