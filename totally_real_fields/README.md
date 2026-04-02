# Data for totally real fields 

This folder includes data giving totally real fields up to some discriminant.  All fields are taken from the LMFDB.

The schema is as follows:

| Column | Type | Description |
| --- | --- | --- |
| lmfdb_index | smallint | The index $i$ of this field on the LMFDB (if it exists) |
| coeffs | numeric[] | A list of ciefficents of a polredabs definnig polynomials, starting with the constant term.  (the leading coefficinet 1 is omitted) |
| mongenic | smallint | Whether the field monogenic: 1 for yes, -1 for no, 0 for not computed |
