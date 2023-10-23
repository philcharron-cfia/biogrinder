import random
import math

def is_int(value):
    try:
        value = int(value)
        return value
    except ValueError:
        raise ValueError(f"{value} is not an integer.")
    
def is_float(value):
    try:
        value = float(value)
        return value
    except ValueError:
        raise ValueError(f"{value} is not a float.")

class OptionNotFoundError(Exception):
    pass

def is_option(value, options):
    if value not in options:
        raise OptionNotFoundError(f"{value} is not in {options}.")
    else:
        return value
     

def normalize(arr, total):
    # Normalize an arrayref to 1.
    if not total:
        raise ValueError("Error: Need to provide a valid total")
    arr = [x / total for x in arr]
    return arr

def identify_sequence_type(seq):
    dna_chars = set("ATCGN-")
    rna_chars = set("AUCGN-")
    protein_chars = set("ACDEFGHIKLMNPQRSTVWY")  # 20 standard amino acids

    # Checking if the sequence only contains characters specific to DNA, RNA, or proteins
    if set(seq.upper()).issubset(dna_chars):
        return "dna"
    elif set(seq.upper()).issubset(rna_chars):
        return "rna"
    elif set(seq.upper()).issubset(protein_chars):
        return "protein"
    else:
        return "Unknown"

def randig(mu, lambda_):
    """
    Random value sampled from the inverse gaussian (a.k.a. Wald) distribution,
    using the method at http://en.wikipedia.org/wiki/Inverse_Gaussian_distribution
    """
    y = random.gauss(0, 1)**2  
    x = mu + (mu**2 * y) / (2 * lambda_) - mu / (2 * lambda_) * math.sqrt(4 * mu * lambda_ * y + mu**2 * y**2)
    if random.random() <= mu / (mu + x):
        y = x
    else:
        y = mu**2 / x
    return y 

