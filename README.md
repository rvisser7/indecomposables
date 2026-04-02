# Database of indecomposable elements over number fields

This is a (work-in-progress) database of additively indecomposable totally positive integral elements over number fields.

For a given a number field $K$, an element $\alpha \in K$ is totally positive if $\sigma(\alpha) > 0$ for all real embeddings $\sigma$ of $K$.  We say that a totally positive integral element $\alpha \in \mathcal{O}_K$ is **(additively) indecomposable** if it cannot be written as the sum of two other totally positive integral elements in $\mathcal{O}_K$.

Given a number field $K$, it's known that there are only finitely many indecomposable elements up to multiplicatiion by totally positive units.  There is however no known general classification of all such indecomposable elements (other than some particular families of number fields, e.g. real quadratic fields. simplest cubic fields, or certain families of biquadratic fields).

Sage code to enumerate all additively indecomposable elements (modulo totally positive units) is provided in the `scripts` folder, tables of data is given in the `data` folder, and various plots of the data is given in the `plots` folder.

Any contributions to either the code or data are very welcome!


## Development

### Setup

This project requires [Sage](https://www.sagemath.org/) for number field computations.

**Installation:**
- **Ubuntu/Debian:** `sudo apt install sagemath`
- **macOS (Homebrew):** `brew install sage`
- **Conda:** `conda install -c conda-forge sage`
- **Docker:** Use the official Sage Docker image

After installing Sage, install development dependencies:
```bash
make dev
```

### Running Tests

Tests are organized into fast and comprehensive suites:

**Fast tests (no Sage required):**
```bash
make test-minimal
```

**Full tests with Sage:**
```bash
make test-sage
```

**All tests:**
```bash
make test
```

**With coverage report:**
```bash
make coverage
```

### Code Quality

**Lint with flake8:**
```bash
make lint
```

**Auto-format with black:**
```bash
make format
```

### Continuous Integration

Tests are automatically run on GitHub Actions for:
- Linting and syntax checks
- Fast minimal tests (on every push/PR)
- Full test suite with Sage (for comprehensive validation)
- Code coverage reporting

See [`.github/workflows/python-app.yml`](.github/workflows/python-app.yml) for workflow details.

### Project Structure

```
indecomposables/
├── scripts/           # Core computation code
│   ├── main.py       # NumberFieldData class
│   ├── real_quadratic.sage      # (TODO) Quadratic field optimization
│   └── simplest_cubic.sage      # (TODO) Cubic field optimization
├── data/             # Database files (pipe-delimited, PostgreSQL-compatible)
├── plots/            # Visualization and analysis
├── tests/            # Test suite
│   ├── test_main.py  # Main tests
│   └── conftest.py   # Pytest configuration
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
