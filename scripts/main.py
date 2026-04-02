# Main Python file for computing sails and indecomposables

from sage.all import (
    Set, ZZ, RR, pi, euler_phi, CyclotomicField, gap, RealField, sqrt, prod,
    QQ, NumberField, PolynomialRing, latex, pari, cached_function, Permutation)

class NumberFieldData:
    """
     Class for storing number field data with sails and indecomposables
    """
    def __init__(self, label):
        self.label = label


