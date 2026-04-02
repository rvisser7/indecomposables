# Scripts for Computing Indecomposables

This folder contains Sage/Python code for computing additively indecomposable elements in totally real number fields.

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

### `real_quadratic.sage`

*TODO:* Implement algorithm for real quadratic fields using Dress-Scharlau's method with continued fraction convergents.

Reference: Dress, A.; Scharlau, W. (1975). "Indecomposable integral representations of finite groups over real quadratic number fields"

### `simplest_cubic.sage`

*TODO:* Implement algorithm for simplest cubic fields using classification by Kala-Tinková and Gil-Muñoz-Tinková.

References:
- Kala, Vítězslav; Tinková, Markéta. "Indecomposable elements in number fields"
- Gil-Muñoz, Daniel; Tinková, Markéta. "On indecomposables in simplest cubic fields"

## Algorithm Details

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

