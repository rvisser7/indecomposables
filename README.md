# Database of indecomposable elements over number fields

This is a (work-in-progress) database of additively indecomposable totally positive integral elements over number fields.

For a given a number field $K$, an element $\alpha \in K$ is totally positive if $\sigma(\alpha) > 0$ for all real embeddings $\sigma$ of $K$.  We say that a totally positive integral element $\alpha \in \mathcal{O}_K$ is **(additively) indecomposable** if it cannot be written as the sum of two other totally positive integral elements in $\mathcal{O}_K$.

Given a number field $K$, it's known that there are only finitely many indecomposable elements up to multiplicatiion by totally positive units.  There is however no known general classification of all such indecomposable elements (other than some particular families of number fields, e.g. real quadratic fields. simplest cubic fields, or certain families of biquadratic fields).

Sage code to enumerate all additively indecomposable elements (modulo totally positive units) is provided in the `scripts` folder, tables of data is given in the `data` folder, and various plots of the data is given in the `plots` folder.

Contributions to either the code or data are very welcome!
