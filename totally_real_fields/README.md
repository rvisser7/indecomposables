# Data for totally real fields 

This folder includes data giving totally real fields up to some discriminant - one file per degree, `degree<n>.txt`. These are the *inputs* to the computation; the results go in `data/`.  All fields are taken from the LMFDB.

These files are deliberately kept small to keep the total file size low. The basic schema is as follows:

| Column | Type | Description |
| --- | --- | --- |
| lmfdb_label | string | The label of this field on the LMFDB (if it exists) |
| lmfdb_index | smallint | The index $i$ of this field on the LMFDB (if it exists) |
| coeffs | numeric[] | A list of coefficients of a polredabs defining polynomials, starting with the constant term.  (the leading coefficient 1 is omitted) |
| monogenic | smallint | Whether the field is monogenic: 1 for yes, -1 for no, 0 for not computed |

Rows are sorted by `|disc|` then `lmfdb_index`.
