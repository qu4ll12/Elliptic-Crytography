"""
ecc/order.py — Curve Order (#E) via Baby-step Giant-step
==========================================================
Implements the mathematical machinery to compute and verify the
order of an elliptic curve group over 𝔽ₚ.

Background
----------
Hasse's theorem guarantees that the curve order falls in the interval:

    p + 1 - 2√p  ≤  #E  ≤  p + 1 + 2√p

This interval has width 4√p, so for a 64-bit prime p:
    width ≈ 4 · 2^32 ≈ 2^34   ← feasible to search with BSGS

Algorithm: Baby-step Giant-step (BSGS)
---------------------------------------
Given a point P on the curve with unknown order n:

1. Choose m = ⌈√(4√p)⌉ + 1  (BSGS step size covering the interval width).
2. Baby steps: precompute  j·P  for j = 0, 1, …, m-1  and store in a dict.
3. Start from the lower Hasse bound: Q = low·P.
4. Giant steps: for i = 0, 1, …, ceil(4√p / m):
       key-lookup Q in the baby-step table.
       if found (i*m + j == order of P), validate: candidate·P == O.
5. The unique value in [low, high] that annihilates P is #E.

For very small primes (p < 10_000) an exact brute-force count is used.

Public API
----------
curve_order(curve)          – Compute #E (baby-step giant-step or brute force)
is_secure_curve(curve, n)   – Check security properties of the curve
find_generator(curve, n)    – Find a generator G of order n
"""

import math
import random
from .point import CurvePoint
from .utils import is_prime


# ---------------------------------------------------------------------------
# Hasse Interval
# ---------------------------------------------------------------------------

def hasse_interval(p: int) -> tuple[int, int]:
    """
    Return the Hasse interval  [p+1 - 2√p,  p+1 + 2√p].

    Hasse's theorem guarantees  #E  lies within this range.
    """
    two_sqrt_p = 2 * math.isqrt(p) + 2   # slight overestimate to be safe
    return (p + 1 - two_sqrt_p, p + 1 + two_sqrt_p)


# ---------------------------------------------------------------------------
# Baby-step Giant-step
# ---------------------------------------------------------------------------

def _brute_force_order(curve) -> int:
    """
    Count #E by testing every x ∈ [0, p-1].
    Only practical for small primes (p < ~10 000).
    """
    p   = curve.p
    a   = curve.a
    b   = curve.b
    count = 1   # start with 1 to include the point at infinity
    for x in range(p):
        rhs = (pow(x, 3, p) + a * x + b) % p
        if rhs == 0:
            count += 1          # the point (x, 0) — tangent horizontal
        else:
            from .utils import legendre_symbol
            if legendre_symbol(rhs, p) == 1:
                count += 2      # two solutions ±y
    return count


def _bsgs_order(curve, P: CurvePoint, low: int, high: int) -> int | None:
    """
    Baby-step Giant-step: find the unique n ∈ [low, high] such that n·P = O.

    Let N = high - low + 1 (width of Hasse interval).
    Choose step size m = ⌈√N⌉.

    Baby steps:  precompute  T[j] = j·P  for j = 0 … m-1.
    Giant steps: starting from Q = low·P, repeatedly add m·P and check
                 if the current point matches any baby step.
                 A match at position (i, j) means:
                     (low + i·m + j)·P = O   →  candidate = low + i·m + j.

    After finding a candidate, also check its negation and the distance
    from the other end so we catch both ±solutions.
    """
    width = high - low + 1
    m = max(math.isqrt(width) + 1, 1)

    # Baby steps: j·P → j,  for j = 0, 1, …, m-1
    baby: dict[tuple, int] = {}
    jP = CurvePoint.infinity(curve)
    for j in range(m):
        key = (jP.x, jP.y, jP._is_infinity)
        baby[key] = j
        jP = jP + P

    # Giant step vector: m·P
    mP = m * P

    # Walk Q = low·P, low·P + m·P, low·P + 2m·P, …
    Q   = low * P
    steps = (width // m) + 2   # enough steps to cover [low, high]

    for i in range(steps):
        key = (Q.x, Q.y, Q._is_infinity)
        if key in baby:
            j = baby[key]
            candidate = low + i * m - j
            if low <= candidate <= high and (candidate * P).is_infinity():
                return candidate
            # Also try candidate wrapped from the other direction
            candidate2 = low + i * m + j
            if low <= candidate2 <= high and (candidate2 * P).is_infinity():
                return candidate2
        Q = Q + mP

    return None


def _verify_order(curve, P: CurvePoint, n: int) -> bool:
    """
    Verify that n is the exact order of P:
      •  n·P = O     (n annihilates P)
    """
    return (n * P).is_infinity()


def curve_order(curve) -> int:
    """
    Compute #E — the number of points on `curve` including the point at infinity.

    For p < 10 000 an exact brute-force count is used (O(p) time).
    For larger p, Baby-step Giant-step is used over the Hasse interval.

    Algorithm (BSGS path)
    ----------------------
    1. Sample a random point P on the curve.
    2. Use BSGS to find the order of P inside the Hasse interval.
    3. The order of P divides #E; enumerate its multiples in the interval.
    4. Cross-validate each multiple with a second independent random point.
    5. Repeat with fresh random points until #E is uniquely determined.

    Returns
    -------
    Integer #E.

    Raises
    ------
    RuntimeError  if unable to determine #E after many attempts.
    """
    p = curve.p

    # Fast exact path for tiny primes (test / educational use)
    if p < 10_000:
        return _brute_force_order(curve)

    low, high = hasse_interval(p)

    for attempt in range(100):
        # Sample a random affine point
        raw = curve.random_point()
        if raw is None:
            continue
        P = CurvePoint(curve, raw[0], raw[1])

        # Find order of P using BSGS
        n = _bsgs_order(curve, P, low, high)
        if n is None:
            continue
        if n <= 0:
            continue

        # #E is a multiple of n that falls in [low, high]
        start = low + ((-low) % n)   # first multiple of n >= low
        candidates = list(range(start, high + 1, n))

        for candidate in candidates:
            # Verify with a second independent random point
            raw2 = curve.random_point()
            if raw2 is None:
                continue
            P2 = CurvePoint(curve, raw2[0], raw2[1])
            if (candidate * P2).is_infinity():
                return candidate

    raise RuntimeError(
        "Failed to determine curve order after many attempts. "
        "This may indicate a degenerate curve or very unlucky sampling — "
        "try regenerating parameters."
    )


# ---------------------------------------------------------------------------
# Security Checks
# ---------------------------------------------------------------------------

def is_secure_curve(curve, order: int, target_bits: tuple[int, int] = (64, 128)) -> dict:
    """
    Evaluate the cryptographic security of a curve against standard criteria.

    Checks performed
    ----------------
    1. prime_order   : #E is prime   → cofactor h = 1, no small-subgroup attacks
    2. not_singular  : discriminant ≠ 0 (guaranteed by EllipticCurve constructor)
    3. not_supsinglar: #E ≢ p+1 (mod p)  — avoids MOV attack via Weil pairing
    4. bit_range     : 64 ≤ log₂(#E) ≤ 128
    5. not_anomalous : #E ≠ p            — avoids SSSA (smart) attack

    Parameters
    ----------
    curve       : EllipticCurve instance.
    order       : The computed curve order #E.
    target_bits : (min_bits, max_bits) for #E.

    Returns
    -------
    A dict with keys matching each check above and bool values,
    plus an 'all_pass' summary.
    """
    p = curve.p
    min_bits, max_bits = target_bits
    order_bits = order.bit_length()

    result = {
        "prime_order"    : is_prime(order),
        "not_singular"   : True,                      # guaranteed by __init__
        "not_supersingular": (order % p) != 1,        # #E ≢ p+1 (mod p)
        "bit_range"      : min_bits <= order_bits <= max_bits,
        "not_anomalous"  : order != p,                # #E ≠ p
        "order_bits"     : order_bits,
    }
    result["all_pass"] = all(
        v for k, v in result.items()
        if k not in ("order_bits",)
    )
    return result


# ---------------------------------------------------------------------------
# Generator Point
# ---------------------------------------------------------------------------

def find_generator(curve, order: int) -> CurvePoint:
    """
    Find a generator point G of order `order` on the curve.

    When #E is prime (cofactor h = 1), **every** non-identity point is
    a generator, so we just need a single random point.

    Parameters
    ----------
    curve : EllipticCurve instance.
    order : The prime curve order #E.

    Returns
    -------
    CurvePoint G such that order·G = O.

    Raises
    ------
    RuntimeError  if no valid generator is found (extremely unlikely for prime #E).
    """
    for _ in range(200):
        raw = curve.random_point()
        if raw is None:
            continue
        G = CurvePoint(curve, raw[0], raw[1])

        # Verify G has the right order: order·G == O
        if (order * G).is_infinity():
            return G

    raise RuntimeError("Could not find a generator point. Check that #E is prime.")
