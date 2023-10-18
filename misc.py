def is_int( value):
    try:
        int(value)
        return True
    except ValueError:
        print(f"Error: {value} is not an integer.")
        return False
    
def is_float( value):
    try:
        float(value)
        return True
    except ValueError:
        print(f"Error: {value} is not a float.")
        return False

class OptionNotFoundError(Exception):
    pass

def is_option(value, options):
    try:
        if value not in options:
            raise OptionNotFoundError(f"Error: {value} is not in {options}.")
        else:
            return ValueError
    except OptionNotFoundError as e:
        print(e)

def normalize(arr, total):
    # Normalize an arrayref to 1.
    if not total:
        raise ValueError("Error: Need to provide a valid total")
    arr = [x / total for x in arr]
    return arr

def proba_cumul(probas):
    sum_val = 0
    cumul_probas = [0]
    for prob in probas:
        sum_val += prob
        cumul_probas.append(sum_val)
    return cumul_probas

def identify_sequence_type(seq):
    dna_chars = set("ATCGN")
    rna_chars = set("AUCGN")
    protein_chars = set("ACDEFGHIKLMNPQRSTVWY")  # 20 standard amino acids

    # Checking if the sequence only contains characters specific to DNA, RNA, or proteins
    if set(seq.upper()).issubset(dna_chars):
        return "DNA"
    elif set(seq.upper()).issubset(rna_chars):
        return "RNA"
    elif set(seq.upper()).issubset(protein_chars):
        return "Protein"
    else:
        return "Unknown"


