# Code to compute indecomposables for real quadratic fields
# Uses the continued fraction convergents, by Dress-Scharlau.
#
# IMPORTANT: The main implementation is now in real_quadratic.py (Python)
# See real_quadratic.py for the RealQuadraticField class.
#
# This .sage file is kept for reference/historical purposes.
# New code should use the Python real_quadratic.py module instead.
#
# Usage:
#   from real_quadratic import RealQuadraticField
#   rq = RealQuadraticField(13)  # Q(sqrt(13))
#   indecomps = rq.compute_indecomposables_dress_scharlau()
#
# Or use with NumberFieldData for automatic optimization:
#   from main import NumberFieldData
#   from sage.all import QuadraticField
#   K = QuadraticField(13)
#   nfd = NumberFieldData(field=K)
#   indecomps = nfd.compute_indecomposables_optimized()
#
# References:
#   Dress, A.; Scharlau, W. (1975).
#   "Indecomposable integral representations of finite groups over real quadratic number fields"
#   Advances in Mathematics, 17(3), 231-273.
#
#   Kala, Vítězslav. "An effective tool for determining indecomposable elements of orders"
#   arXiv:2108.15387

# See real_quadratic.py for complete, modern implementation

