r"""The Schur–Weyl measure.

Measure $k$ i.i.d. copies of a state $\rho$ in the Schur basis and you get a
random partition $\lambda \vdash k$, with

$$
\Pr[\lambda] \;=\; f^\lambda \, s_\lambda(\operatorname{spec} \rho),
$$

where $s_\lambda$ is the Schur polynomial. This is the distribution behind
spectrum estimation and quantum tomography. Partitions with more rows than
$\rho$ has eigenvalues get probability zero.

Examples:
    >>> from schur_weyl.sw_measure import schur_weyl_measure
    >>> {lam: round(float(p), 4) + 0.0 for lam, p in schur_weyl_measure([0.7, 0.3], k=4).items()}
    {(4,): 0.4141, (3, 1): 0.4977, (2, 2): 0.0882, (2, 1, 1): 0.0, (1, 1, 1, 1): 0.0}
"""

import numpy as np
from .young_diagrams import partitions
from .dimensions import dim_specht
from .symmetric_functions import schur_polynomial


def schur_weyl_measure(spectrum: list[float], k: int) -> dict[tuple[int, ...], float]:
    """The probabilities $tr(\\Pi^\\lambda rho^{\\otimes k})$ where rho has the defined spectrum

    Args:
        spectrum (list[float]): the spectrum of rho
        k (int): the number of copies of rho

    Returns:
        dict[tuple[int, ...], float]: The probabilities and their associated partition.
    """
    ret: dict[tuple[int, ...], float] = dict()
    for lam in partitions(k):
        ret[lam] = dim_specht(lam) * schur_polynomial(lam, spectrum)
    return ret
