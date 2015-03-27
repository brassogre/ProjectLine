'''
Created on Feb 3, 2015

@author: 10034888
'''
import decimal
decimal.getcontext().prec = 50
decimal.getcontext().rounding = decimal.ROUND_HALF_UP
import math

def pi():
    """Compute Pi to the current precision.

    >>> print(pi())
    3.141592653589793238462643383

    """
    decimal.getcontext().prec += 2  #  extra digits for intermediate steps
    three = decimal.Decimal(3)      #  substitute "three=3.0" for regular floats
    lasts, t, s, n, na, d, da = 0, three, 3, 1, 0, 0, 24
    while s != lasts:
        lasts = s
        n, na = n + na, na + 8
        d, da = d + da, da + 32
        t = (t * n) / d
        s += t
    decimal.getcontext().prec -= 2
    return +s               #  unary plus applies the new precision

def exp(x):
    """Return e raised to the power of x.  Result type matches input type.

    >>> print(exp(Decimal(1)))
    2.718281828459045235360287471
    >>> print(exp(Decimal(2)))
    7.389056098930650227230427461
    >>> print(exp(2.0))
    7.38905609893
    >>> print(exp(2+0j))
    (7.38905609893+0j)

    """
    decimal.getcontext().prec += 2
    i, lasts, s, fact, num = 0, 0, 1, 1, 1
    while s != lasts:
        lasts = s
        i += 1
        fact *= i
        num *= x
        s += num / fact
    decimal.getcontext().prec -= 2
    return +s               #  unary plus applies the new precision

def cos(x):
    """Return the cosine of x as measured in radians.

    The Taylor series approximation works best for a small value of x.
    For larger values, first compute x = x % (2 * pi).

    >>> print(cos(Decimal('0.5')))
    0.8775825618903727161162815826
    >>> print(cos(0.5))
    0.87758256189
    >>> print(cos(0.5+0j))
    (0.87758256189+0j)

    """
    decimal.getcontext().prec += 2
    i, lasts, s, fact, num, sign = 0, 0, 1, 1, 1, 1
    while s != lasts:
        lasts = s
        i += 2
        fact *= i * (i - 1)
        num *= x * x
        sign *= -1
        s += num / fact * sign
    decimal.getcontext().prec -= 2
    return +s               #  unary plus applies the new precision

def sin(x):
    """Return the sine of x as measured in radians.

    The Taylor series approximation works best for a small value of x.
    For larger values, first compute x = x % (2 * pi).

    >>> print(sin(Decimal('0.5')))
    0.4794255386042030002732879352
    >>> print(sin(0.5))
    0.479425538604
    >>> print(sin(0.5+0j))
    (0.479425538604+0j)

    """
    decimal.getcontext().prec += 2
    i, lasts, s, fact, num, sign = 1, 0, x, 1, x, 1
    while s != lasts:
        lasts = s
        i += 2
        fact *= i * (i - 1)
        num *= x * x
        sign *= -1
        s += num / fact * sign
    decimal.getcontext().prec -= 2
    return +s               #  unary plus applies the new precision

def atan2(x):
    pass

def acos(x):
    return 2.*atan2((decimal.Decimal(1. - x)).sqrt(), (decimal.Decimal(1. + x)).sqrt())



if __name__ == '__main__':
    Pi = pi()
    print(sin(Pi / 2))
    print(sin(Pi / 6))
    print(math.asin(decimal.Decimal(1)) * 2)
    print(cos(Pi / 6))
    print(decimal.Decimal(3).sqrt() / 2)

    import fractions
    ratpi = fractions.Fraction(Pi).limit_denominator(1000)
    print(ratpi)
    Ratnum = decimal.Decimal(ratpi.numerator)
    Ratden = decimal.Decimal(ratpi.denominator)
    Ratpi = Ratnum / Ratden
    print(Ratpi)
    print(Pi)
