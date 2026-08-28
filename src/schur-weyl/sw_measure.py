import numpy as np
from young_diagrams import partitions
from dimensions import dim_specht
from symmetric_functions import schur_polynomial


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
