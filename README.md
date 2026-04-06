# Database of indecomposable elements over number fields

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://github.com/codespaces/new/rvisser7/indecomposables)

This is a (work-in-progress) database of additively indecomposable totally positive integral elements over number fields.

For a given a number field $K$, an element $\alpha \in K$ is totally positive if $\sigma(\alpha) > 0$ for all real embeddings $\sigma$ of $K$.  We say that a totally positive integral element $\alpha \in \mathcal{O}_K$ is **(additively) indecomposable** if it cannot be written as the sum of two other totally positive integral elements in $\mathcal{O}_K$.

Given a number field $K$, it's known that there are only finitely many indecomposable elements up to multiplicatiion by totally positive units.  There is however no known general classification of all such indecomposable elements (other than some particular families of number fields, e.g. real quadratic fields. simplest cubic fields, or certain families of biquadratic fields).

Sage code to enumerate all additively indecomposable elements (modulo totally positive units) is provided in the `scripts` folder, tables of data is given in the `data` folder, and various plots of the data is given in the `plots` folder.

Any contributions to either the code or data are very welcome!


## Usage

After installing Sage (see Development section below), you can compute indecomposables for totally real number fields using the provided shell script.

### Basic Usage

Compute indecomposables for degree 3 fields with discriminants between 1 and 100:

```bash
./compute_indecomposables.sh 3 1 100
```

This will:
- Process all degree 3 fields from `totally_real_fields/deg3.txt`
- Filter fields with discriminants in the range [1, 100]
- Output results to `results_deg3.txt`
- Show basic progress information

### With Verbose Output

For detailed computation progress, enable verbose mode:

```bash
./compute_indecomposables.sh 3 1 100 1
```

Verbose mode shows:
- Discriminant values for each field
- Progress through norm checking (percentage complete)
- When indecomposables are found
- Total count of indecomposables per field

### Parameters

```
./compute_indecomposables.sh DEGREE DISC_MIN DISC_MAX [VERBOSE]
```

- `DEGREE`: Degree of number fields to process (3, 4, 5, etc.)
- `DISC_MIN`: Minimum discriminant value (default: 1)
- `DISC_MAX`: Maximum discriminant value (required)
- `VERBOSE`: Optional, set to `1` to enable verbose output

### Examples

```bash
# Degree 4 fields, discriminants 1-500, quiet mode
./compute_indecomposables.sh 4 1 500

# Degree 5 fields, discriminants 100-1000, verbose mode
./compute_indecomposables.sh 5 100 1000 1

# Large discriminant range for degree 3
./compute_indecomposables.sh 3 1 10000 1
```


## Development

### Setup

This project requires [Sage](https://www.sagemath.org/) for number field computations.

### Project Structure

```
indecomposables/
├── scripts/           # Core computation code
│   ├── main.py       # NumberFieldData class with optimized algorithm selection
│   ├── real_quadratic.py      # Dress-Scharlau algorithm for quadratic fields
│   ├── simplest_cubic.py      # Kala-Tinková classification for cubic fields
│   ├── demo_real_quadratic.py # Performance demo for quadratic fields
│   └── demo_simplest_cubic.py # Performance demo for cubic fields
├── data/             # Database files (pipe-delimited, PostgreSQL-compatible)
├── plots/            # Visualization and analysis
├── tests/            # Test suite
│   ├── test_main.py          # Main NumberFieldData tests
│   ├── test_real_quadratic.py # RealQuadraticField tests
│   └── test_simplest_cubic.py # SimplestCubicField tests
└── Makefile          # Development commands
```

### Main API

See [`scripts/README.md`](scripts/README.md) for detailed API documentation of the `NumberFieldData` class.

Quick example:
```python
from sage.all import QQ, PolynomialRing
from main import NumberFieldData

# Create a number field
R = PolynomialRing(QQ, 'x')
K = QQ.extension(R.gen()**3 - 2, names='a')

# Compute indecomposables
nfd = NumberFieldData(label="3.3.49.1", field=K)
indecomps = nfd.compute_indecomposables()

print(f"Found {len(indecomps)} indecomposables")
```
