from taylor_series import *
import numpy as np

FILE = '../../LamTermsK7-M3-C3-lam10.0-carRate0.2-W1.1-stop1e-04.dat'

def main(args):
    terms = np.loadtxt(FILE)