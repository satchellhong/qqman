import numbers

import numpy as np


def ppoints(n, a=None):
    """numpy analogue or `R`'s `ppoints` function
    see details at https://www.rdocumentation.org/packages/stats/versions/3.6.2/topics/ppoints
    https://docs.tibco.com/pub/enterprise-runtime-for-R/5.0.0/doc/html/Language_Reference/stats/ppoints.html
    :param n: array type or number"""

    n = np.float64(n) if isinstance(n, numbers.Number) else np.float64(len(n))
    if a is None:
        a = 0.375 if n <= 10 else 0.5

    return (np.arange(n) + 1 - a) / (n + 1 - 2 * a)
