# Development Guide: NumberFieldData Class

This document provides a comprehensive reference for the `NumberFieldData` class in `main.py`, which is the core API for computing indecomposables in totally real number fields.

## Overview

The `NumberFieldData` class represents a totally real number field and provides methods to:
- Compute indecomposable elements
- Analyze units and sails
- Export results in database format

**Key Properties:**
- `K`: The Sage `NumberField` object
- `discriminant`: Discriminant of the field
- `regulator`, `class_number`, `degree`: Field invariants
- `indecomposables`: List of computed indecomposables (up to units)
- `fundamental_units`, `totally_positive_units`: Unit groups
- `unit_representatives`: All unit classes modulo squares

## Initialization

```python
from main import NumberFieldData
from sage.all import PolynomialRing, QQ

# Create a number field via defining polynomial
R = PolynomialRing(QQ, 'x')
poly = [1, -2, 0, 1]  # x^3 - 2
K = QQ.extension(R(poly), names='a')

# Initialize NumberFieldData
nfd = NumberFieldData(
    coeffs=poly,
    lmfdb_label="3.3.49.1",    # Optional: LMFDB label
    lmfdb_index="1",           # Optional: last part of label
    metadata={                  # Optional: precomputed field data
        'degree': 3,
        'discriminant': 49,
        'regulator': 2.5,
        'class_number': 1
    }
)
```

## Properties (Read-Only)

### Field Invariants

**`degree`**
- Returns the degree of the number field
- Computed from `K.degree()` or metadata
- Type: `int`

**`discriminant`**
- Returns the discriminant of the number field
- Computed on first access if not in metadata
- Type: `Integer`

**`regulator`**
- Returns the regulator of the unit group
- Computed on first access using Sage's `K.regulator()`
- Type: `float`

**`class_number`**
- Returns the class number of the field
- Computed on first access using Sage's `K.class_number()`
- Type: `int`

**`indecomposables`**
- Returns the list of indecomposables computed by `compute_indecomposables()`
- Returns `None` if not yet computed
- Elements are sorted by norm, then trace
- Type: `list` or `None`

### Unit Groups

**`fundamental_units`**
- Returns a set of fundamental units for the field
- Auto-computed on first access via `_compute_fundamental_units()`
- Units normalized to be positive under the first real embedding
- Type: `list` of field elements

**`totally_positive_units`**
- Returns a basis for the totally positive unit group
- Auto-computed on first access via `_compute_totally_positive_units()`
- Type: `list` of field elements
- Vector space dimension: `degree - 1`

**`unit_representatives`**
- Returns one representative from each unit class modulo squares
- Auto-computed on first access via `_compute_unit_representatives()`
- Type: `list` of field elements
- Number of representatives: $2^{\text{degree}}$

## Public Methods

### Main Computation

**`compute_indecomposables(verbose=True)`**

Compute all indecomposables up to multiplication by totally positive units using the optimal algorithm for the field type.

**Behavior:**
- For degree 2 (real quadratic fields): Uses the fast Dress-Scharlau continued fraction algorithm
- For degree 3 (simplest cubics): Uses the Kala-Tinková classification if applicable
- For degree 4 (certain biquadratics): Uses Siu Hang Man's classification if applicable
- For other fields: Falls back to brute force enumeration

**Parameters:**
- `verbose` (bool): Print progress information (default: `True`)

**Returns:**
- List of indecomposables sorted by norm then trace
- Elements are Sage number field elements in `K`

**Example:**
```python
nfd = NumberFieldData(coeffs=[1, -2, 0, 1], lmfdb_label="3.3.49.1")
indecomps = nfd.compute_indecomposables(verbose=True)
print(f"Found {len(indecomps)} indecomposables")
for elem in indecomps[:5]:
    print(f"  {elem} (norm: {K(elem).norm()})")
```

---

**`compute_indecomposables_brute_force(verbose=True, algorithm="new", GRH=False)`**

Compute indecomposables using brute force enumeration over ideals of bounded norm.

**Description:**
Uses the ideal bound brute force method based on the Kala-Yatsyna bound. Enumerates all ideals of norm up to the bound, and for each principal ideal with a totally positive generator, checks if it can be decomposed as a sum of already-found indecomposables.

**Parameters:**
- `verbose` (bool): Print progress information (default: `True`)
- `algorithm` (str): Version of algorithm to use; "new" recommended (default: `"new"`)
- `GRH` (bool): Whether to assume the Generalized Riemann Hypothesis for faster computation (default: `False`)

**Returns:**
- List of indecomposables sorted by norm then trace

**Example:**
```python
# Force brute force even if faster method available
nfd.compute_indecomposables_brute_force(verbose=True, GRH=False)

# With GRH assumption (faster but less rigorous)
nfd.compute_indecomposables_brute_force(GRH=True)
```

---

**`is_indecomposable_milp(x, eps=None, precision=None, require_integral=True, return_witness=False)`**

Test whether a given element is indecomposable using mixed-integer linear programming (MILP).

**Description:**
Tests if there exists a totally positive element $y = \sum_i k_i w_i$ (with $k_i \in \mathbb{Z}$ and $w_i$ an integral basis) such that:
$$0 < \sigma_j(y) < \sigma_j(x)$$
for all real embeddings $\sigma_j$ of $K$.

If such $y$ exists, then $x$ is decomposable (since $x = y + (x - y)$ with both summands totally positive). Otherwise, $x$ is indecomposable.

**Parameters:**
- `x`: Element of `K` to test
- `eps` (float): Tolerance for strict inequalities (default: auto-computed based on precision)
- `precision` (int): Real field precision in bits (default: `self.precision`)
- `require_integral` (bool): Enforce that `x ∈ O_K` (default: `True`)
- `return_witness` (bool): If `True`, return a decomposition witness when decomposable (default: `False`)

**Returns:**
- If `return_witness=False`: `bool` — `True` if indecomposable, `False` if decomposable
- If `return_witness=True`: `tuple` — `(is_indecomposable, witness_or_None)`
  - `witness_or_None` is an element `y` such that `x - y` is totally positive, or `None`

**Example:**
```python
# Simple check
is_indecomp = nfd.is_indecomposable_milp(nfd.K.gen())

# Get witness for decomposable element
is_indecomp, witness = nfd.is_indecomposable_milp(
    x=some_element,
    return_witness=True
)
if not is_indecomp:
    print(f"Decomposable as {witness} + ({some_element} - {witness})")

# High precision check
is_indecomp = nfd.is_indecomposable_milp(x, precision=200)
```

---

**`compute_sails(signature, verbose=True, GRH=False)`**

*Under development* — Will compute the "sails" of the field for a given signature.

---

**`to_data_row()`**

Format the computed results as a pipe-delimited data row for database output.

**Description:**
Returns a formatted string suitable for appending to a database file. Format:
```
LMFDB_LABEL|DISCRIMINANT|REGULATOR|CLASS_NUMBER|NUM_INDECOMPS|MIN_NORM|MAX_NORM|[other fields]
```

**Requires:**
- `compute_indecomposables()` must have been called first

**Returns:**
- String with pipe-delimited columns

**Example:**
```python
nfd.compute_indecomposables()
row = nfd.to_data_row()
print(row)
# Output: 3.3.49.1|49|2.5|1|5|2|13|...

# Write to file
with open('results_deg3.txt', 'a') as f:
    f.write(row + '\n')
```

## Configuration Options

These attributes control algorithm behavior:

**`self.precision`** (int)
- Real field precision in bits for numerical embeddings
- Default: `100`
- Used in unit computation and MILP checks
- Example: `nfd.precision = 150`

---

**`self.max_norm_bound`** (int or None)
- Manual override for the bound on indecomposable norms
- If `None`, computed automatically based on degree and discriminant
- Default: `None`
- Example: `nfd.max_norm_bound = 1000`

---

**`self.check_asserts`** (bool)
- Enable expensive assertion checks (slows computation)
- Default: `False`
- Example: `nfd.check_asserts = True`

---

**`self.exit_for_nonunit`** (bool)
- Exit brute force early upon finding first non-unit indecomposable
- Default: `False`
- Example: `nfd.exit_for_nonunit = True`

---

**Example: Configuration**
```python
nfd = NumberFieldData(coeffs=[1, -2, 0, 1])
nfd.precision = 150        # Higher precision
nfd.max_norm_bound = 500   # Limit search
nfd.check_asserts = True   # Enable validation
indecomps = nfd.compute_indecomposables()
```

## Internal Methods (Helpers)

These methods are called automatically but documented for understanding:

**`_compute_fundamental_units()`**
- Computes fundamental units via Sage's `unit_group()`
- Caches result in `_fundamental_units`

**`_compute_totally_positive_units()`**
- Computes totally positive unit basis using GF(2) kernel computation
- Also sets `_unit_signature_rank`

**`_compute_unit_representatives()`**
- Builds representatives of all unit classes modulo squares
- Result has size $2^{\text{degree}}$

**`_compute_max_norm_bound()`**
- Computes the Kala-Yatsyna bound for maximal indecomposable norm
- Includes optimizations for specific degree/discriminant ranges

**`_signature(x)`**
- Returns the signature of element `x` as tuple of ±1 values
- One ±1 per real embedding

## Usage Examples

### Example 1: Basic Computation

```python
from main import NumberFieldData
from sage.all import PolynomialRing, QQ

R = PolynomialRing(QQ, 'x')
x = R.gen()

# Define Q(cbrt(2))
K = QQ.extension(x**3 - 2, names='a')

# Compute indecomposables
nfd = NumberFieldData(coeffs=[1, -2, 0, 1], lmfdb_label="3.3.49.1")
indecomps = nfd.compute_indecomposables(verbose=True)

print(f"Degree: {nfd.degree}")
print(f"Discriminant: {nfd.discriminant}")
print(f"Found {len(indecomps)} indecomposables")
```

### Example 2: Batch Processing

```python
from main import NumberFieldData

def process_fields_from_file(filename, degree, disc_min, disc_max):
    """Process all fields in a data file within discriminant range."""
    results = []
    with open(filename) as f:
        for line in f:
            parts = line.strip().split(':')
            coeffs = [int(c) for c in parts[0].split(',')]
            disc = int(parts[1])
            label = parts[2]
            
            if disc_min <= disc <= disc_max:
                nfd = NumberFieldData(
                    coeffs=coeffs,
                    lmfdb_label=label,
                    metadata={'discriminant': disc}
                )
                indecomps = nfd.compute_indecomposables(verbose=False)
                results.append((label, len(indecomps)))
    
    return results
```

### Example 3: High-Precision Verification

```python
# Verify indecomposability with multiple precisions
nfd = NumberFieldData(coeffs=[...])
nfd.compute_indecomposables()

for elem in nfd.indecomposables[:3]:
    for prec in [100, 150, 200]:
        is_indecomp = nfd.is_indecomposable_milp(elem, precision=prec)
        print(f"  {elem}: precision={prec}, indecomposable={is_indecomp}")
```

### Example 4: Unit Group Analysis

```python
nfd = NumberFieldData(coeffs=[...])

# Access unit groups
fund_units = nfd.fundamental_units
print(f"Fundamental units: {fund_units}")

tp_units = nfd.totally_positive_units
print(f"Totally positive units (rank {len(tp_units)}): {tp_units}")

# Unit representatives
unit_reps = nfd.unit_representatives
print(f"Unit classes modulo squares: {len(unit_reps)} representatives")
```

## Performance Tips

1. **Use degree-specific algorithms**: `compute_indecomposables()` automatically selects the fastest algorithm for each degree.

2. **Parallelize**: Use `--num-threads` in `compute_indecomposables.sh` to split work across processes.

3. **Precision tuning**: 
   - Higher precision (e.g., `200` bits) for rigorous verification
   - Lower precision (e.g., `100` bits) for rapid screening

4. **GRH assumption**: Use `GRH=True` in brute force to speed computation (less rigorous).

5. **Early exit**: Set `exit_for_nonunit=True` to stop upon finding first non-unit indecomposable.

## References

- Kala, Vítězslav; Yatsyna, Pavlo. "On Kitaoka's conjecture and lifting problem for universal quadratic forms." *Mathematika* 65.3 (2019): 612-630.
- Dress-Scharlau algorithm: Fast continued fraction method for quadratic fields
- Kala-Tinková classification: Efficient algorithm for simplest cubic fields
