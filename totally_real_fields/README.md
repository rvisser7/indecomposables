# Data for totally real fields 

This folder includes data giving totally real fields up to some discriminant.  All fields are taken from the LMFDB.

The schema is as follows:

| Column | Type | Description |
| --- | --- | --- |
| lmfdb_label | string | The label of this field on the LMFDB (if it exists) |
| lmfdb_index | smallint | The index $i$ of this field on the LMFDB (if it exists) |
| coeffs | numeric[] | A list of coefficents of a polredabs defining polynomials, starting with the constant term.  (the leading coefficient 1 is omitted) |
| monogenic | smallint | Whether the field is monogenic: 1 for yes, -1 for no, 0 for not computed |
