import numpy as np
import os
import re
import sys
import math
from scipy.stats import norm

def hist(data, min_val=None, max_val=None):
    if not data:
        raise ValueError("Error: no data provided to hist()")
    # If min or max are not provided, use the min and max of the data
    min_val = min(data) if min_val is None else min_val
    max_val = max(data) if max_val is None else max_val
    # Count occurrences of each integer in the range [min_val, max_val]
    histogram = [data.count(x) for x in range(min_val, max_val+1)]
    return histogram

def normal(x_min, x_max, mean, variance, num):
    # Evaluate the normal distribution across the range x_min to x_max
    distribution = norm.pdf(np.arange(x_min, x_max+1), mean, math.sqrt(variance)) * num
    return distribution.tolist()

def uniform(x_min, x_max, min_range, max_range, num):
    # Evaluate the uniform distribution across the range x_min to x_max
    width = max_range - min_range + 1
    uniform_distribution = [num / width if min_range <= x <= max_range else 0 for x in range(x_min, x_max+1)]
    return uniform_distribution

def corr_coeff(y, f, mean):
    SSerr = sum((yi - fi) ** 2 for yi, fi in zip(y, f))
    SStot = sum((yi - mean) ** 2 for yi in y)
    R2 = 1 - (SSerr / SStot)
    return R2

def blockPrint():
    sys.stdout = open(os.devnull, 'w')

def error_positions(read):
    error_positions = []
    match = re.search(r'errors=(\S+)', read.description)
    
    if match:
        err_str = match.group(1)
        errors = err_str.split(',')
        for error in errors:
            match = re.match(r'(\d+)([%+-])([a-z]*)', error, re.I)
            if match:
                pos = int(match.group(1))
                error_positions.append(pos)

    return error_positions

def get_references(read):
    """
    Get the number of references that a read comes from.
    """
    desc = read.description
    refs = desc.split('reference=')[1].split(' ')[0]
    return refs.split(',')