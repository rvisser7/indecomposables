# Scripts for Computing Indecomposables

This folder contains Sage/Python code for computing additively indecomposable elements in totally real number fields.

## 📖 Documentation

For a comprehensive reference of the `NumberFieldData` class including all methods, properties, and usage examples, see [DEVELOPMENT.md](DEVELOPMENT.md).

## Files

### `main.py`

Core implementation with the `NumberFieldData` class for computing indecomposables.

**Main Class: `NumberFieldData`**

Represents a totally real number field and computes its indecomposable elements using the ideal bound brute force method.

**Key Properties:**
- `label`: LMFDB label of the field
- `K`: The Sage NumberField object
- `discriminant`: Discriminant of K
- `regulator`: Regulator of K
- `class_number`: Class number of K
- `degree`: Degree of K
- `indecomposables`: List of computed indecomposables (up to units)
- `fundamental_units`: Fundamental units of K
- `totally_positive_units`: Basis for totally positive units
- `unit_representatives`: All unit classes modulo squares

**Main Methods:**

- `compute_indecomposables(verbose=True)`: Main algorithm to find all indecomposables
  - Enumerates ideals of bounded norm
  - For each principal ideal with totally positive generator, checks decomposability
  - Returns list sorted by norm then trace

- `to_data_row()`: Format as pipe-delimited output row for database

**Configuration:**
```python
nfd = NumberFieldData(field=K, label="3.3.49.1")
nfd.exit_for_nonunit = False    # Exit early if non-unit indecomp found
nfd.check_asserts = False       # Run expensive assertions
nfd.max_norm_bound = None       # Override default norm bound
nfd.precision = 100             # RealField precision
```

**Usage Example:**

```python
from sage.all import QQ, PolynomialRing

# Define a cubic field via polynomial
R = PolynomialRing(QQ, 'x')
x = R.gen()
K = QQ.extension(x**3 - 2, names='a')

# Create NumberFieldData and compute
nfd = NumberFieldData(label="3.3.49.1", field=K)
indecomps = nfd.compute_indecomposables(verbose=True)

# Access results
print(f"Found {len(indecomps)} indecomposables")
for ind in indecomps:
    print(f"  {ind} (norm {K(ind).norm()})")

# Output to database format
row = nfd.to_data_row()
print(row)
```

**Utility Functions:**

- `create_field_from_polynomial(poly, var='a')`: Create NumberFieldData from polynomial string
- `process_degree_file(filename)`: Parse data file and create NumberFieldData objects
- `batch_compute_indecomposables(fields, degree=3, verbose=True, output_file=None)`: Process multiple fields

### `real_quadratic.py`

Efficient implementation of Dress-Scharlau algorithm for real quadratic fields.

**Class: `RealQuadraticField`**

Computes indecomposables in Q(√D) using continued fraction convergents, much faster than the general brute force.

**Key Methods:**
- `compute_indecomposables_dress_scharlau(verbose=True)` - Main computation
- Properties: `fund_unit`, `tp_unit`, `cf_data`, discriminant, regulator, etc.

**Key Advantages:**
- ⚡ 10-100x faster than brute force for quadratic fields
- Uses mathematical structure (continued fractions) rather than enumeration
- Satisfies Dress-Scharlau bound: max norm ≤ disc(K)/4

**Usage:**
```python
from real_quadratic import RealQuadraticField

# Q(sqrt(13))
rq = RealQuadraticField(13)
indecomps = rq.compute_indecomposables_dress_scharlau()
print(f"Found {len(indecomps)} indecomposables")

# Or use with NumberFieldData for automatic optimization
from main import NumberFieldData
from sage.all import QuadraticField

K = QuadraticField(13)
nfd = NumberFieldData(field=K)
indecomps = nfd.compute_indecomposables_optimized()  # Uses Dress-Scharlau automatically
```

**Algorithm Overview:**
1. Compute continued fraction expansion of (√D - 1)/2 or √D depending on D mod 4
2. Generate convergents p_i/q_i using standard CF algorithm
3. Compute alpha elements: α_i = p_i + q_i·δ (where δ = (√D ± 1)/2)
4. Extract candidates: α_{i,t} = α_i + t·α_{i+1} for all valid i, t
5. Normalize up to multiplication by totally positive units
6. Return unique representatives sorted by norm

**References:**
- Dress, A.; Scharlau, W. (1975). "Indecomposable integral representations of finite groups over real quadratic number fields". Advances in Mathematics, 17(3), 231-273.
- Kala, Vítězslav. "An effective tool for determining indecomposable elements of orders". arXiv:2108.15387

### `simplest_cubic.py`

Efficient implementation of Kala-Tinková classification for simplest cubic fields.

**Class: `SimplestCubicField`**

Computes indecomposables in simplest cubic fields Q[x]/(x³ - n*x² - (n+3)*x - 1) where Z[ρ] has index 1 or 3 in the maximal order, using explicit formulas from Kala-Tinková and Gil-Muñoz-Tinková.

**Key Methods:**
- `compute_indecomposables_kala_tinkova(verbose=True)` - Main computation using explicit classification
- Properties: `discriminant`, `regulator`, `index`, `can_use_classification`

**Key Advantages:**
- ⚡ Extremely fast: direct enumeration from formulas (no search required)
- Uses mathematical structure rather than brute force enumeration
- Handles both index 1 (O_K = Z[ρ]) and index 3 ([O_K : Z[ρ]] = 3) cases
- Provides exact theoretical bounds on counts and maximal norms

**Usage:**
```python
from simplest_cubic import SimplestCubicField

# Simplest cubic field with n=1: Q[x]/(x^3 - x^2 - 4*x - 1)
scf = SimplestCubicField(1)
indecomps = scf.compute_indecomposables_kala_tinkova()
print(f"Found {len(indecomps)} indecomposables")

# Or use with NumberFieldData for automatic optimization
from main import NumberFieldData
from sage.all import QQ, PolynomialRing

R = PolynomialRing(QQ, 'x')
K = QQ.extension(x**3 - x**2 - 4*x - 1, names='rho')
nfd = NumberFieldData(field=K)
indecomps = nfd.compute_indecomposables_optimized()  # Uses Kala-Tinková automatically
```

**Algorithm Overview:**
1. **Field Classification**: Check if Z[ρ] has index 1 or 3 in O_K
2. **Index 1 Case** (O_K = Z[ρ]): Enumerate from Theorem 1.2 formulas with parameters v, w
3. **Index 3 Case** ([O_K : Z[ρ]] = 3): Enumerate from Theorem 1.1 with 8 different cases and parameters v, r
4. **Normalization**: Remove duplicates up to multiplication by totally positive units
5. **Verification**: Check against theoretical predictions for counts and maximal norms

**References:**
- Kala, Vítězslav; Tinková, Markéta. "Universal quadratic forms, small norms and traces in families of number fields". Theorem 1.2 (index 1 case)
- Gil-Muñoz, Daniel; Tinková, Markéta. "Additive structure of non-monogenic simplest cubic fields". Theorem 1.1 (index 3 case)

## Algorithm Details

### Dress-Scharlau Method (Real Quadratic Fields)

For Q(√D):

1. **Continued Fraction**: Compute CF of (√D - 1)/2 or √D with period s
2. **Generate Alphas**: Use convergents to build α_i via recurrence relation
3. **Extract Indecomposables**: From α_{i,t} = α_i + t·α_{i+1}
4. **Normalize**: Adjust by powers of totally positive units to get canonical representatives
5. **Filter**: Keep only unique representatives

**Complexity:** O(s² × num_candidates) where s is CF period length (typically small).

**Bounds:** Max indecomposable norm ≤ |disc(K)|/4

### Ideal Bound Brute Force Method

For a totally real number field K of degree d:

1. **Compute Unit Group**: Find fundamental units and basis for totally positive units
2. **Set Norm Bound**: Use Kala-Yatsyna bound: max norm ≤ disc(K)
3. **Enumerate Ideals**: For each norm n from 1 to bound:
   - Get all principal ideals of norm n
   - Find totally positive generator (if exists)
   - Normalize generator using log-space matrices
4. **Check Decomposability**: For each candidate α, check if α = β + γ where β is previous indecomp and γ is totally positive
   - Use Cramer's rule on pairs of embeddings to bound coefficients in terms of totally positive units
   - Test all bounded vectors

### Complexity

- Time: O(max_norm × num_ideals × |previous_indecomposables| × d³)
- Space: O(|indecomposables| × d)

Practical bounds make this tractable for degree ≤ 6 and discriminant ≤ 10^6.

## Data Format

Output follows pipe-delimited format compatible with PostgreSQL import via psycodict:

```
lmfdb_label|discriminant|regulator|class_number|num_indecomposables|min_norm|max_norm|indecomposables
```

Example:
```
3.3.49.1|49|0.525454682122572|4|1|2|7|[1, -a^2 + 4]
```

## References

- Kala, Vítězslav; Yatsyna, Pavlo (2023). "On Kitaoka's conjecture and lifting problem for universal quadratic forms". Bull. Lond. Math. Soc. 55(2), 854–864.
- Brunotte, Horst (1983). "Zur Zerlegung totalpositiver Zahlen in Ordnungen totalreeller algebraischer Zahlkörper". Arch. Math. 41(6), 502–503.
- Brunotte, Horst (1982). "The computation of a certain metric invariant of an algebraic number field". Math. Comp. 38(158), 627–632.

