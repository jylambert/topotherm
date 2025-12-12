"""Provides helper functions common to all models."""


def annuity(c_i, n):
    """
    Calculate the annuity factor.

    Parameters
    ----------
    c_i : float
        Interest rate.
    n : float
        Number of years.

    Returns
    -------
    float
        Annuity factor.

    Examples
    --------
    Calculate the annuity factor for a 5% interest rate over 10 years:

    >>> annuity(0.05, 10)
    0.129504...

    With a higher interest rate (10%) and the same duration:

    >>> annuity(0.10, 10)
    0.162745...
    """
    a = ((1 + c_i) ** n * c_i) / ((1 + c_i) ** n - 1)
    return a
