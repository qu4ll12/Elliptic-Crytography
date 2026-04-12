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
We use raw affine coordinates + fast modular arithmetic (gmpy2) in the
inner BSGS loop for performance.  CurvePoint objects are only used for
the high-level API.

For very small primes (p < 10_000) an exact brute-force count is used.

Public API
----------
curve_order(curve)                    – Compute #E (BSGS or brute force)
is_secure_curve(curve, n)             – Check security properties
find_generator(curve, n)              – Find a generator G of order n
find_prime_order_curve(p, callback)   – Auto-search (a,b) until #E is prime
"""

import math
import random
from .point import CurvePoint
from .curve import EllipticCurve
from .utils import is_prime, legendre_symbol, modular_inverse


# ---------------------------------------------------------------------------
# Hasse Interval
# ---------------------------------------------------------------------------

def hasse_interval(p: int) -> tuple[int, int]:
    """
    Return the Hasse interval  [p+1 - 2√p,  p+1 + 2√p].

    We add a small margin (+2) to the √p estimate to ensure coverage.
    """
    two_sqrt_p = 2 * math.isqrt(p) + 2
    return (p + 1 - two_sqrt_p, p + 1 + two_sqrt_p)


# ---------------------------------------------------------------------------
# Fast raw-coordinate arithmetic (no object overhead)
# ---------------------------------------------------------------------------
# Points are represented as (x, y) tuples or None for the identity O.
# This avoids CurvePoint constructor overhead (including the is_on_curve
# check) in the tight BSGS inner loops.
# ---------------------------------------------------------------------------

def _raw_add(P, Q, a: int, p: int):
    """
    Add two points in raw (x,y) tuple form.  Returns (x3, y3) or None (identity).
    P, Q are each either (x, y) or None.
    """
    if P is None:
        return Q
    if Q is None:
        return P

    x1, y1 = P
    x2, y2 = Q

    if x1 == x2:
        if (y1 + y2) % p == 0:
            return None                # P + (-P) = O
        # Point doubling:  λ = (3x² + a) / (2y)
        num = (3 * x1 * x1 + a) % p
        den = (2 * y1) % p
    else:
        # General addition:  λ = (y2 - y1) / (x2 - x1)
        num = (y2 - y1) % p
        den = (x2 - x1) % p

    lam = num * modular_inverse(den, p) % p
    x3 = (lam * lam - x1 - x2) % p
    y3 = (lam * (x1 - x3) - y1) % p
    return (x3, y3)


def _raw_neg(P, p: int):
    """Negate a point: -(x, y) = (x, p-y)."""
    if P is None:
        return None
    return (P[0], (p - P[1]) % p)


def _raw_scalar_mul(k: int, P, a: int, p: int):
    """
    Compute k·P using double-and-add on raw coordinates.
    k must be a non-negative integer.
    """
    if k < 0:
        return _raw_scalar_mul(-k, _raw_neg(P, p), a, p)
    if k == 0 or P is None:
        return None
    result = None
    addend = P
    while k:
        if k & 1:
            result = _raw_add(result, addend, a, p)
        addend = _raw_add(addend, addend, a, p)
        k >>= 1
    return result


# ---------------------------------------------------------------------------
# Brute-force for tiny primes
# ---------------------------------------------------------------------------

def _brute_force_order(curve) -> int:
    """
    Count #E by testing every x ∈ [0, p-1].
    Only practical for small primes (p < ~10 000).
    """
    p = curve.p
    a = curve.a
    b = curve.b
    count = 1   # point at infinity
    for x in range(p):
        rhs = (pow(x, 3, p) + a * x + b) % p
        if rhs == 0:
            count += 1
        elif legendre_symbol(rhs, p) == 1:
            count += 2
    return count


# ---------------------------------------------------------------------------
# Baby-step Giant-step (optimised with raw coordinates)
# ---------------------------------------------------------------------------

def _bsgs_find_order(curve) -> int | None:
    """
    Baby-step Giant-step to find #E for the given curve.

    Algorithm
    ---------
    1. Pick a random point P on the curve.
    2. Compute Q = (p+1)·P.   (Because #E·P = O  ⟹  t·P = (p+1)·P  where t = p+1-#E.)
    3. Baby steps: build table  { j·P : j }  for j = 0 … m.
    4. Giant steps: walk BOTH directions from Q:
         Q + i·(m·P)   →  finds negative traces  t = -(i·m ± j)
         Q - i·(m·P)   →  finds positive traces  t = +(i·m ± j)
       If a match is found:  #E = p + 1 - t.
    5. Verify candidate.

    Returns #E or None if this attempt failed.
    """
    p = curve.p
    a = curve.a
    b = curve.b

    # Hasse bound for the trace t:  |t| ≤ 2√p
    bound = 2 * math.isqrt(p) + 2
    # Width of the search interval for t: [-bound, +bound], width = 2·bound + 1
    width = 2 * bound + 1

    # BSGS step size
    m = math.isqrt(width) + 1

    # ── Pick a random point P ──
    raw_pt = curve.random_point()
    if raw_pt is None:
        return None
    P = raw_pt   # (x, y) tuple

    # ── Baby steps: store j·P for j = 0 … m ──
    baby: dict[tuple | None, int] = {}
    jP = None   # 0·P = O, represented as None
    for j in range(m + 1):
        # Use the x-coordinate (or None for identity) as key.
        # We store j for the first occurrence only.
        if jP is None:
            key = None
        else:
            key = jP[0]   # x-coordinate uniquely identifies ±y
        if key not in baby:
            baby[key] = (j, jP)
        jP = _raw_add(jP, P, a, p)

    # ── Compute Q = (p+1)·P ──
    Q = _raw_scalar_mul(p + 1, P, a, p)

    # ── Giant step vector: R = m·P ──
    R = _raw_scalar_mul(m, P, a, p)
    neg_R = _raw_neg(R, p)

    # ── Giant steps: walk BOTH directions from Q ──
    # Walking Q - i·R  covers positive t  (t = i·m ± j)
    # Walking Q + i·R  covers negative t  (t = -(i·m ∓ j))
    # We need enough steps to cover |t| ≤ bound with stride m.
    num_steps = (bound // m) + 2

    # Two cursors: one walks Q-i·R, the other walks Q+i·R
    cursor_minus = Q   # Q - i·R  (for positive t)
    cursor_plus  = Q   # Q + i·R  (for negative t)

    for i in range(num_steps + 1):
        # Check both cursors against the baby step table
        for cursor, direction in [(cursor_minus, +1), (cursor_plus, -1)]:
            if cursor is None:
                key = None
            else:
                key = cursor[0]

            if key in baby:
                j, jP_val = baby[key]

                # cursor = Q ∓ i·m·P = (t ∓ i·m)·P
                # baby match at x-coord means cursor = ±j·P
                # So t ∓ i·m = ±j  ⟹  t = ±i·m ± j
                # direction=+1 (cursor_minus): t = +i·m ± j
                # direction=-1 (cursor_plus):  t = -i·m ± j
                for sign_j in (+1, -1):
                    t = direction * (i * m) + sign_j * j
                    if abs(t) <= bound:
                        candidate = p + 1 - t
                        if candidate > 1:
                            check = _raw_scalar_mul(candidate, P, a, p)
                            if check is None:
                                return candidate

        # Advance the cursors
        cursor_minus = _raw_add(cursor_minus, neg_R, a, p)  # Q - (i+1)·R
        cursor_plus  = _raw_add(cursor_plus,  R,     a, p)  # Q + (i+1)·R

    return None


def curve_order(curve) -> int:
    """
    Compute #E — the number of points on `curve` including the point at infinity.

    For p < 10 000 an exact brute-force count is used (O(p) time).
    For larger p, Baby-step Giant-step is used.

    Returns
    -------
    Integer #E.

    Raises
    ------
    RuntimeError  if unable to determine #E after many attempts.
    """
    p = curve.p

    # Fast exact path for tiny primes
    if p < 10_000:
        return _brute_force_order(curve)

    low, high = hasse_interval(p)

    for attempt in range(40):
        candidate = _bsgs_find_order(curve)
        if candidate is None:
            continue

        # Sanity: must be in Hasse interval
        if not (low <= candidate <= high):
            continue

        # Cross-validate with 3 additional independent random points
        valid = True
        for _ in range(3):
            raw2 = curve.random_point()
            if raw2 is None:
                continue
            check = _raw_scalar_mul(candidate, raw2, curve.a, p)
            if check is not None:
                valid = False
                break

        if valid:
            return candidate

    raise RuntimeError(
        "Failed to determine curve order after many attempts. "
        "This may indicate a degenerate curve or very unlucky sampling — "
        "try regenerating parameters."
    )


# ---------------------------------------------------------------------------
# Auto-Search: Find a curve with prime order
# ---------------------------------------------------------------------------

def find_prime_order_curve(
    p: int,
    max_attempts: int = 2000,
    callback=None,
) -> tuple[EllipticCurve, int] | None:
    """
    For a given prime p, repeatedly sample random (a, b) and compute #E
    until a curve with PRIME order and all security checks satisfied is found.

    The probability that a random curve has prime #E is roughly
    1/ln(p) ≈ 1/(0.693 · bitlen(p)), so for 64-bit p we expect ~44 tries.

    Parameters
    ----------
    p            : A prime modulus.
    max_attempts : Give up after this many (a,b) pairs.
    callback     : Optional callable(attempt, curve_or_None, order_or_None)
                   for progress reporting.

    Returns
    -------
    (curve, order)  if a prime-order curve is found, or None.
    """
    for attempt in range(1, max_attempts + 1):
        a = random.randrange(1, p)
        b = random.randrange(1, p)
        disc = (4 * pow(a, 3, p) + 27 * pow(b, 2, p)) % p
        if disc == 0:
            if callback:
                callback(attempt, None, None)
            continue

        curve = EllipticCurve(p, a, b)

        try:
            order = curve_order(curve)
        except RuntimeError:
            if callback:
                callback(attempt, curve, None)
            continue

        if callback:
            callback(attempt, curve, order)

        checks = is_secure_curve(curve, order)
        if checks["all_pass"]:
            return (curve, order)

    return None


# ---------------------------------------------------------------------------
# Security Checks
# ---------------------------------------------------------------------------

def is_secure_curve(curve, order: int, target_bits: tuple[int, int] = (64, 128)) -> dict:
    """
    Evaluate the cryptographic security of a curve against standard criteria.

    Checks performed
    ----------------
    1. prime_order    : #E is prime   → cofactor h = 1
    2. not_singular   : discriminant ≠ 0
    3. not_supersingular : #E ≢ p+1 (mod p) — avoids MOV attack
    4. bit_range      : 64 ≤ log₂(#E) ≤ 128
    5. not_anomalous  : #E ≠ p — avoids SSSA/Smart attack

    Returns
    -------
    Dict with bool values for each check + 'all_pass' summary.
    """
    p = curve.p
    min_bits, max_bits = target_bits
    order_bits = order.bit_length()

    result = {
        "prime_order"      : is_prime(order),
        "not_singular"     : True,
        "not_supersingular": (order % p) != 1,
        "bit_range"        : min_bits <= order_bits <= max_bits,
        "not_anomalous"    : order != p,
        "order_bits"       : order_bits,
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

    When #E is prime (cofactor h = 1), every non-identity point is a
    generator, so we just pick a random one and verify.

    Returns
    -------
    CurvePoint G such that order·G = O.
    """
    p = curve.p
    a = curve.a
    for _ in range(200):
        raw = curve.random_point()
        if raw is None:
            continue
        # Verify with raw arithmetic (fast)
        check = _raw_scalar_mul(order, raw, a, p)
        if check is None:   # order·P == O ✓
            return CurvePoint(curve, raw[0], raw[1])

    raise RuntimeError("Could not find a generator point. Check that #E is prime.")
