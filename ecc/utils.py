"""
ecc/utils.py — Number-Theory Primitives (gmpy2-accelerated)
============================================================
Core math tools needed by the rest of the library.

Uses gmpy2 (GNU Multiple Precision) for fast modular arithmetic
on large integers (64–128 bit), with pure-Python fallbacks.

Public API
----------
is_prime(n)               – Primality test (gmpy2.is_prime or Miller-Rabin)
generate_prime(bits)      – Random prime of exact bit-length
modular_inverse(a, m)     – Modular inverse (gmpy2.invert or extended GCD)
legendre_symbol(a, p)     – Euler's criterion     a^((p-1)/2) mod p
tonelli_shanks(n, p)      – Square root of n modulo prime p
"""

import random

# ---------------------------------------------------------------------------
# Try to import gmpy2 for fast modular arithmetic (C-level speed).
# Fall back to pure Python if unavailable.
# ---------------------------------------------------------------------------
try:
    import gmpy2
    _HAS_GMPY2 = True
except ImportError:
    _HAS_GMPY2 = False


# ---------------------------------------------------------------------------
# Primality
# ---------------------------------------------------------------------------

def is_prime(n: int, rounds: int = 25) -> bool:
    """
    Primality test.

    Uses gmpy2.is_prime() if available (deterministic for small n,
    BPSW + extra Miller-Rabin rounds for large n).
    Falls back to Miller-Rabin with `rounds` random bases.

    Parameters
    ----------
    n      : Integer to test.
    rounds : Number of Miller-Rabin rounds (only used without gmpy2).

    Returns
    -------
    True if n is (almost certainly) prime.
    """
    if n < 2:
        return False

    if _HAS_GMPY2:
        # gmpy2.is_prime returns 0 (not prime), 1 (probably prime), 2 (definitely prime)
        return gmpy2.is_prime(int(n), rounds) > 0

    # ── Pure-Python Miller-Rabin fallback ──
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
        x = pow(a, d, n)

        if x == 1 or x == n - 1:
            continue

        for _ in range(r - 1):
            x = pow(x, 2, n)
            if x == n - 1:
                break
        else:
            return False

    return True


def generate_prime(bits: int) -> int:
    """
    Generate a random prime whose bit-length is exactly `bits`.

    Strategy: repeatedly sample a random odd integer of the correct
    bit-length and run primality tests until one passes.

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
    Compute the modular inverse of `a` modulo `m`.

    Uses gmpy2.invert() if available (C-speed), otherwise the
    Extended Euclidean Algorithm.

    Returns x such that  a * x ≡ 1 (mod m).
    Raises  ValueError if gcd(a, m) ≠ 1 (inverse does not exist).
    """
    a = a % m
    if a == 0:
        raise ValueError(f"Modular inverse does not exist: gcd(0, {m}) = {m}")

    if _HAS_GMPY2:
        try:
            result = gmpy2.invert(gmpy2.mpz(a), gmpy2.mpz(m))
            return int(result)
        except ZeroDivisionError:
            raise ValueError(f"Modular inverse does not exist: gcd({a}, {m}) ≠ 1")

    # ── Pure-Python Extended Euclidean Algorithm fallback ──
    g, x, _ = _extended_gcd(a, m)
    if g != 1:
        raise ValueError(f"Modular inverse does not exist: gcd({a}, {m}) = {g}")
    return x % m


def _extended_gcd(a: int, b: int):
    """
    Extended Euclidean Algorithm (iterative, avoids recursion limit).
    Returns (gcd, x, y) such that a*x + b*y = gcd(a, b).
    """
    old_r, r = a, b
    old_s, s = 1, 0
    while r != 0:
        quotient = old_r // r
        old_r, r = r, old_r - quotient * r
        old_s, s = s, old_s - quotient * s
    # old_r = gcd, old_s = x, compute y
    if b != 0:
        y = (old_r - old_s * a) // b
    else:
        y = 0
    return old_r, old_s, y


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
    if _HAS_GMPY2:
        result = gmpy2.legendre(gmpy2.mpz(a), gmpy2.mpz(p))
        return int(result)

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
        r = pow(n, (p + 1) // 4, p)
        return min(r, p - r)

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
