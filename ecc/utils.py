"""
ecc/utils.py — Number-Theory Primitives
========================================
Pure-Python implementations of the core math tools needed by the rest
of the library.  No third-party dependencies are required.

Public API
----------
is_prime(n, rounds)       – Miller-Rabin probabilistic primality test
generate_prime(bits)      – Random prime of exact bit-length
modular_inverse(a, m)     – Modular inverse via Extended Euclidean Algorithm
legendre_symbol(a, p)     – Euler's criterion     a^((p-1)/2) mod p
tonelli_shanks(n, p)      – Square root of n modulo prime p
"""

import random


# ---------------------------------------------------------------------------
# Primality
# ---------------------------------------------------------------------------

def is_prime(n: int, rounds: int = 20) -> bool:
    """
    Miller-Rabin probabilistic primality test.

    Parameters
    ----------
    n      : Integer to test.
    rounds : Number of random bases to try.  20 gives a false-positive
             probability < 4^{-20} ≈ 10^{-12}.

    Returns
    -------
    True  – n is *probably* prime.
    False – n is definitely composite.
    """
    if n < 2:
        return False
    # Small prime short-circuit
    small_primes = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37]
    if n in small_primes:
        return True
    if any(n % p == 0 for p in small_primes):
        return False

    # Write n-1 as 2^r * d with d odd
    r, d = 0, n - 1
    while d % 2 == 0:
        r += 1
        d //= 2

    # Witness loop
    for _ in range(rounds):
        a = random.randrange(2, n - 1)
        x = pow(a, d, n)          # a^d mod n  (Python built-in is fast)

        if x == 1 or x == n - 1:
            continue

        for _ in range(r - 1):
            x = pow(x, 2, n)
            if x == n - 1:
                break
        else:
            return False           # Composite witness found

    return True


def generate_prime(bits: int) -> int:
    """
    Generate a random prime whose bit-length is exactly `bits`.

    Strategy: repeatedly sample a random odd integer of the correct
    bit-length and run Miller-Rabin until one passes.

    Parameters
    ----------
    bits : Desired bit-length (must be ≥ 2).

    Returns
    -------
    A prime p with 2^(bits-1) ≤ p < 2^bits.
    """
    if bits < 2:
        raise ValueError("bits must be >= 2")

    low  = 1 << (bits - 1)          # 2^(bits-1)  — smallest bits-bit number
    high = (1 << bits) - 1          # 2^bits - 1  — largest  bits-bit number

    while True:
        candidate = random.randrange(low | 1, high, 2)  # always odd
        if is_prime(candidate):
            return candidate


# ---------------------------------------------------------------------------
# Modular Arithmetic
# ---------------------------------------------------------------------------

def modular_inverse(a: int, m: int) -> int:
    """
    Compute the modular inverse of `a` modulo `m` using the
    Extended Euclidean Algorithm.

    Returns x such that  a * x ≡ 1 (mod m).
    Raises  ValueError if gcd(a, m) ≠ 1 (inverse does not exist).
    """
    g, x, _ = _extended_gcd(a % m, m)
    if g != 1:
        raise ValueError(f"Modular inverse does not exist: gcd({a}, {m}) = {g}")
    return x % m


def _extended_gcd(a: int, b: int):
    """
    Extended Euclidean Algorithm.
    Returns (gcd, x, y) such that a*x + b*y = gcd(a, b).
    """
    if a == 0:
        return b, 0, 1
    g, x1, y1 = _extended_gcd(b % a, a)
    return g, y1 - (b // a) * x1, x1


# ---------------------------------------------------------------------------
# Quadratic Residues
# ---------------------------------------------------------------------------

def legendre_symbol(a: int, p: int) -> int:
    """
    Compute the Legendre symbol (a | p) using Euler's criterion:

        (a | p) = a^((p-1)/2)  mod p

    Returns
    -------
     1  – a is a non-zero quadratic residue mod p
     0  – p divides a
    -1  – a is a quadratic non-residue mod p
    """
    ls = pow(a, (p - 1) // 2, p)
    return -1 if ls == p - 1 else ls


def tonelli_shanks(n: int, p: int) -> int | None:
    """
    Tonelli-Shanks algorithm: find x such that x² ≡ n (mod p).

    Parameters
    ----------
    n : The number whose square root we want.
    p : A prime modulus.

    Returns
    -------
    x  (the smaller of the two square roots) if n is a QR mod p,
    or None if no square root exists.
    """
    if legendre_symbol(n, p) != 1:
        return None          # n has no square root mod p

    # Trivial cases
    if n == 0:
        return 0
    if p == 2:
        return n % p
    if p % 4 == 3:
        return pow(n, (p + 1) // 4, p)

    # Factor out powers of 2 from p-1: p-1 = Q * 2^S
    Q, S = p - 1, 0
    while Q % 2 == 0:
        Q //= 2
        S += 1

    # Find a quadratic non-residue z
    z = 2
    while legendre_symbol(z, p) != -1:
        z += 1

    M = S
    c = pow(z, Q, p)
    t = pow(n, Q, p)
    R = pow(n, (Q + 1) // 2, p)

    while True:
        if t == 0:
            return 0
        if t == 1:
            return min(R, p - R)

        # Find the least i such that t^(2^i) ≡ 1 (mod p)
        i, tmp = 1, pow(t, 2, p)
        while tmp != 1:
            tmp = pow(tmp, 2, p)
            i += 1

        b  = pow(c, 1 << (M - i - 1), p)
        M  = i
        c  = pow(b, 2, p)
        t  = (t * c) % p
        R  = (R * b) % p
