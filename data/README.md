# Tables of indecomposables

Data on the indecomposables for each number field is given in this folder.  Fields are seperated by their degree - currently only totally real number fields are computed in the database.  Currently, all fields given in our data are also in the [LMFDB](https://www.lmfdb.org/NumberField/).

Each file is formatted to be easily uploadable into a PostgreSQL database using [psycodict](https://github.com/roed314/psycodict). The first line always consists of the column headers, the second line is the postgres data type for each corresponding column, and the third line is empty.  The fourth line onwards contains the indecomposables data, given in the same order as the column headers (separated by `|`).

The schema for each file is as follows:

| Column | Type | Description |
| --- | --- | --- |
| lmfdb_label | string | the LMFDB label of the number field K |
| unit_signature_rank | smallint | the unit signature rank of the number field |
| num_indecomposables | int | the number of indecomposables, up to multiplication by totally positive units |
| min_norm_indecomposable | int | the smallest norm of a non-unit indecomposable in K, if it exists.  Otherwise, this is -1. |
| max_norm_indecomposable | int | the maximal norm of an indecomposable in K |
| indecomposables | int[] | a list of all indecomposables of K, up to multiplication by totally positive units |

