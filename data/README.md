# Tables of indecomposables

Data on the indecomposables for each number field is given in this folder.  Fields are seperated by their degree - currently only totally real number fields are computed in the database.

The schema for each file is as follows:

| Column | Type | Description |
| --- | --- | --- |
| lmfdb_label | string | the LMFDB label of the number field K |
| unit_signature_rank | smallint | the unit signature rank of the number field |
| num_indecomposables | int | the number of indecomposables, up to multiplication by totally positive units |
| min_norm_indecomposable | int | the smallest norm of a non-unit indecomposable in K, if it exists.  Otherwise, this is 0 |
| max_norm_indecomposable | int | the maximal norm of an indecomposable in K |
| indecomposables | int[] | a list of all indecomposables of K, up to multiplication by totally positive units |

