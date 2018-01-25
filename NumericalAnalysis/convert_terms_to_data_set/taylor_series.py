import os
import numpy as np

def get_expansion(terms=[], point0=0, values=0):
    result = np.zeros(np.size(values))
    values = np.array(values)
    for n, term in enumerate(terms):
        result += term*np.power(values-point0, n)
    return result
    