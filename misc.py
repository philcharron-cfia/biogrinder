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
