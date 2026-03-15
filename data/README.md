# Tables of additively indecomposable elements

Data on the indecomposables for each number field is given in this folder.  Fields are seperated by their degree - currently only totally real number fields are computed in the database.  Currently, all fields given in our data are also in the [LMFDB](https://www.lmfdb.org/NumberField/).

Each file is formatted to be easily uploadable into a PostgreSQL database using [psycodict](https://github.com/roed314/psycodict). The first line always consists of the column headers, the second line is the postgres data type for each corresponding column, and the third line is empty.  The fourth line onwards contains the indecomposables data, given in the same order as the column headers (separated by `|`).

The schema for each file is as follows.  Some of the columns below are already present in the LMFDB, though we give them here again for convenience (to easily compute stats and plots).  We note that not all columns are guaranteed to be present in each data file.

| Column | Type | Description |
| --- | --- | --- |
| lmfdb_label | string | the LMFDB label of the number field K |
| discriminant | numeric | the discriminant of the number field K |
| regulator | numeric | the regulator of the number field K |
| class_number | numeric | the class number of the number field K |
| is_monogenic | smallint | whether the field K is monogenic |
| num_subfields | smallint | the number of subfields of K, excluding Q and K itself |
| unit_signature_rank | smallint | the unit signature rank of the number field |
| positive_unit_basis | jsonb | a free basis for the group of totally positive units |
| num_indecomposables | numeric | the number of indecomposables, up to multiplication by totally positive units |
| min_norm_indecomposable | numeric | the smallest norm of a non-unit indecomposable in K, if it exists.  Otherwise, this is -1. |
| max_norm_indecomposable | numeric | the maximal norm of an indecomposable in K |
| min_norm_all_signatures | numeric | the smallest positive integer k such that every signature has an element of norm at most k |
| min_norm_fprimea_signature | numeric | the smallest norm of an element in the signature containing f'(a), where f is a defining polynomial for K |
| indecomposables | numeric[] | a list of all indecomposables of K, up to multiplication by totally positive units |

