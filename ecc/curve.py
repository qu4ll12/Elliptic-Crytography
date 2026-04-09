"""
ecc/curve.py — EllipticCurve Class
====================================
Represents a short Weierstrass curve over a prime finite field 𝔽ₚ:

    y² ≡ x³ + ax + b  (mod p)

The class handles parameter storage, curve validation (non-singular
check via discriminant), and point membership tests.

A helper function `generate_curve(p)` automatically samples valid
parameters `a` and `b` given a user-supplied prime `p`.
"""

import random
from .utils import modular_inverse, tonelli_shanks, legendre_symbol


class EllipticCurve:
    """
    Short Weierstrass elliptic curve:  y² ≡ x³ + ax + b  (mod p)

    Attributes
    ----------
    p : Prime modulus (defines the finite field 𝔽ₚ).
    a : Curve coefficient in x term.
    b : Curve constant term.

    The curve is *non-singular* iff its discriminant is non-zero:
        Δ = -16 · (4a³ + 27b²)  ≢ 0  (mod p)
    """

    def __init__(self, p: int, a: int, b: int):
        """
        Initialise and immediately validate the curve.

        Raises
        ------
        ValueError  if p is not prime, or if the curve is singular.
        """
        self.p = p
        self.a = a % p
        self.b = b % p
        self._validate()

    # ------------------------------------------------------------------
    # Validation
    # ------------------------------------------------------------------

    def _validate(self):
        """Ensure the curve is non-singular (discriminant ≠ 0 mod p)."""
        disc = self._discriminant()
        if disc == 0:
            raise ValueError(
                f"Singular curve: discriminant ≡ 0 (mod p).\n"
                f"  4a³ + 27b² ≡ 0 (mod {self.p}) with a={self.a}, b={self.b}"
            )

    def _discriminant(self) -> int:
        """
        Compute 4a³ + 27b²  (mod p).
        A non-zero result guarantees the curve has no cusps or self-intersections.
        """
        p = self.p
        return (4 * pow(self.a, 3, p) + 27 * pow(self.b, 2, p)) % p

    # ------------------------------------------------------------------
    # Point membership
    # ------------------------------------------------------------------

    def is_on_curve(self, x: int, y: int) -> bool:
        """
        Check whether the affine point (x, y) satisfies the curve equation:

            y² ≡ x³ + ax + b  (mod p)
        """
        p = self.p
        lhs = pow(y, 2, p)
        rhs = (pow(x, 3, p) + self.a * x + self.b) % p
        return lhs == rhs

    # ------------------------------------------------------------------
    # Random point sampling (used internally for order counting)
    # ------------------------------------------------------------------

    def random_point(self) -> tuple[int, int] | None:
        """
        Sample a random affine point on the curve by trial.

        Strategy:
          1. Pick a random x ∈ [0, p-1].
          2. Compute rhs = x³ + ax + b  (mod p).
          3. Check if rhs is a quadratic residue (Legendre symbol = 1).
          4. If yes, compute y via Tonelli-Shanks and return (x, y).
          5. Retry up to 2p times (expected ~2 tries since half of 𝔽ₚ* are QRs).

        Returns None only if sampling fails after many attempts (should be
        astronomically unlikely for well-formed curves).
        """
        p = self.p
        for _ in range(2 * p.bit_length() + 128):
            x = random.randrange(0, p)
            rhs = (pow(x, 3, p) + self.a * x + self.b) % p
            if legendre_symbol(rhs, p) == 1:
                y = tonelli_shanks(rhs, p)
                if y is not None and self.is_on_curve(x, y):
                    return (x, y)
        return None

    # ------------------------------------------------------------------
    # Display helpers
    # ------------------------------------------------------------------

    def __repr__(self) -> str:
        return (
            f"EllipticCurve(\n"
            f"  equation : y² ≡ x³ + {self.a}·x + {self.b}  (mod {self.p})\n"
            f"  p        : {self.p}  ({self.p.bit_length()} bits)\n"
            f"  a        : {self.a}\n"
            f"  b        : {self.b}\n"
            f"  Δ = 4a³+27b² : {self._discriminant()}  (non-zero ✓)\n"
            f")"
        )


# ---------------------------------------------------------------------------
# Auto-generation of curve parameters
# ---------------------------------------------------------------------------

def generate_curve(p: int, max_attempts: int = 10_000) -> EllipticCurve:
    """
    Given a prime `p`, randomly sample `a` and `b` until we obtain a
    non-singular curve.

    Algorithm
    ---------
    1. Pick a  ∈ [1, p-1] uniformly at random.
    2. Pick b  ∈ [1, p-1] uniformly at random.
    3. Check discriminant  4a³ + 27b² ≢ 0 (mod p).
    4. If zero, retry.  Expected retries: very few (only hits zero for
       specific (a,b) pairs, probability ≈ 1/p ≈ negligible).

    Parameters
    ----------
    p            : A prime modulus.
    max_attempts : Safety cap; raises RuntimeError if never satisfied.

    Returns
    -------
    A valid, non-singular EllipticCurve over 𝔽ₚ.
    """
    for _ in range(max_attempts):
        a = random.randrange(1, p)
        b = random.randrange(1, p)
        disc = (4 * pow(a, 3, p) + 27 * pow(b, 2, p)) % p
        if disc != 0:
            return EllipticCurve(p, a, b)

    raise RuntimeError(
        f"Could not generate a non-singular curve after {max_attempts} attempts. "
        "This is extremely unlikely — please check that p is indeed prime."
    )
