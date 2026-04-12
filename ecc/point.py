"""
ecc/point.py — Elliptic Curve Point Arithmetic
================================================
Implements the group law for points on a short Weierstrass curve
using *affine coordinates*.

    y² ≡ x³ + ax + b  (mod p)

The group operation rules:
  • Point at infinity O  is the identity element (O + P = P).
  • Negation:  -(x, y) = (x, -y mod p).
  • Addition:  P + Q  via the chord formula (P ≠ Q).
  • Doubling:  P + P  via the tangent formula.
  • Scalar multiplication:  k·P  via Double-and-Add (O(log k) operations).
"""

from __future__ import annotations
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from .curve import EllipticCurve

from .utils import modular_inverse


# ---------------------------------------------------------------------------
# Sentinel: Point at Infinity
# ---------------------------------------------------------------------------

class _Infinity:
    """Singleton representing the point at infinity O (identity element)."""
    _instance = None

    def __new__(cls):
        if cls._instance is None:
            cls._instance = super().__new__(cls)
        return cls._instance

    def __repr__(self):
        return "O  (point at infinity)"

    def __eq__(self, other):
        return isinstance(other, _Infinity)

    def __hash__(self):
        return hash("infinity")


INFINITY = _Infinity()


# ---------------------------------------------------------------------------
# CurvePoint
# ---------------------------------------------------------------------------

class CurvePoint:
    """
    An affine point (x, y) on an EllipticCurve, or the point at infinity O.

    Do NOT construct O directly — use the module-level constant `INFINITY`
    or the class method `CurvePoint.infinity(curve)`.

    Examples
    --------
    >>> G = CurvePoint(curve, gx, gy)
    >>> two_G = G + G
    >>> k_G = 42 * G       # scalar multiplication
    >>> neg_G = -G
    """

    def __init__(self, curve: "EllipticCurve", x: int, y: int):
        """
        Create a point and verify it lies on the curve.

        Parameters
        ----------
        curve : The EllipticCurve this point belongs to.
        x, y  : Affine coordinates (integers in [0, p-1]).

        Raises
        ------
        ValueError  if (x, y) does not satisfy the curve equation.
        """
        self.curve = curve
        self.x = x % curve.p
        self.y = y % curve.p
        self._is_infinity = False

        if not curve.is_on_curve(self.x, self.y):
            raise ValueError(
                f"Point ({x}, {y}) is NOT on the curve:\n  {curve}"
            )

    @classmethod
    def infinity(cls, curve: "EllipticCurve") -> "CurvePoint":
        """Return a CurvePoint instance representing O (identity element)."""
        obj = object.__new__(cls)
        obj.curve = curve
        obj.x = None
        obj.y = None
        obj._is_infinity = True
        return obj

    # ------------------------------------------------------------------
    # Identity check
    # ------------------------------------------------------------------

    def is_infinity(self) -> bool:
        """Return True if this point is the identity element O."""
        return self._is_infinity

    # ------------------------------------------------------------------
    # Negation:  -(x, y) = (x, p - y)
    # ------------------------------------------------------------------

    def __neg__(self) -> "CurvePoint":
        """
        Negate the point:  -(x, y) = (x, -y mod p).
        The negation of O is O.
        """
        if self.is_infinity():
            return self
        return CurvePoint(self.curve, self.x, self.curve.p - self.y)

    # ------------------------------------------------------------------
    # Addition / Doubling
    # ------------------------------------------------------------------

    def __add__(self, other: "CurvePoint") -> "CurvePoint":
        """
        Add two points using the standard chord-and-tangent formulas.

        Cases handled:
          1. P + O = P               (identity)
          2. O + Q = Q               (identity)
          3. P + (-P) = O            (inverse)
          4. P + P   = double(P)     (tangent formula)
          5. P + Q   (general case)  (chord formula)
        """
        if not isinstance(other, CurvePoint):
            return NotImplemented
        if self.curve.p != other.curve.p:
            raise ValueError("Cannot add points from different curves.")

        # Case 1 & 2: identity element
        if self.is_infinity():
            return other
        if other.is_infinity():
            return self

        p = self.curve.p

        # Case 3: P + (-P) = O
        if self.x == other.x and (self.y + other.y) % p == 0:
            return CurvePoint.infinity(self.curve)

        # Case 4: P == Q  →  point doubling
        if self.x == other.x and self.y == other.y:
            return self._double()

        # Case 5: General addition (chord formula)
        #   slope λ = (y₂ - y₁) / (x₂ - x₁)  mod p
        dy = (other.y - self.y) % p
        dx = (other.x - self.x) % p
        lam = (dy * modular_inverse(dx, p)) % p

        x3 = (lam * lam - self.x - other.x) % p
        y3 = (lam * (self.x - x3) - self.y) % p
        return CurvePoint(self.curve, x3, y3)

    def _double(self) -> "CurvePoint":
        """
        Point doubling using the tangent formula.

        Slope of tangent at P = (x, y):
            λ = (3x² + a) / (2y)  mod p
        """
        if self.is_infinity():
            return self

        p = self.curve.p
        a = self.curve.a

        # λ = (3x² + a) * (2y)^{-1}  mod p
        numerator   = (3 * pow(self.x, 2, p) + a) % p
        denominator = (2 * self.y) % p
        lam = (numerator * modular_inverse(denominator, p)) % p

        x3 = (lam * lam - 2 * self.x) % p
        y3 = (lam * (self.x - x3) - self.y) % p
        return CurvePoint(self.curve, x3, y3)

    # ------------------------------------------------------------------
    # Scalar Multiplication — Double-and-Add  O(log k)
    # ------------------------------------------------------------------

    def __rmul__(self, k: int) -> "CurvePoint":
        """
        Compute k·P using the Double-and-Add algorithm.

        Algorithm (left-to-right binary method):
        -----------------------------------------
            Q = O
            for each bit of k from MSB to LSB:
                Q = Q + Q         # double
                if bit == 1:
                    Q = Q + P     # add

        Complexity: O(log k) point additions/doublings.

        Examples
        --------
        >>> 3 * G  ==  G + G + G
        >>> 0 * G  ==  O
        """
        if not isinstance(k, int):
            return NotImplemented

        # NOTE: Do NOT reduce k modulo p here.  The scalar must be the exact
        # integer because the group order #E ≈ p (by Hasse) and reducing
        # mod p would silently corrupt BSGS order-counting which uses k ≈ p.
        if k < 0:
            return (-k) * (-self)   # negate both scalar and point

        if k == 0:
            return CurvePoint.infinity(self.curve)
        if k == 1:
            return CurvePoint(self.curve, self.x, self.y)
        if self.is_infinity():
            return self

        result = CurvePoint.infinity(self.curve)   # Q = O
        addend = CurvePoint(self.curve, self.x, self.y)  # R = P

        while k:
            if k & 1:                 # if current bit is 1
                result = result + addend
            addend = addend + addend  # double
            k >>= 1                   # shift to next bit

        return result

    def __mul__(self, k: int) -> "CurvePoint":
        """Support P * k  (delegates to k * P)."""
        return self.__rmul__(k)

    # ------------------------------------------------------------------
    # Comparison
    # ------------------------------------------------------------------

    def __eq__(self, other) -> bool:
        if not isinstance(other, CurvePoint):
            return False
        if self.is_infinity() and other.is_infinity():
            return True
        if self.is_infinity() or other.is_infinity():
            return False
        return self.x == other.x and self.y == other.y

    def __hash__(self):
        return hash((self.x, self.y, self.curve.p))

    # ------------------------------------------------------------------
    # Display
    # ------------------------------------------------------------------

    def __repr__(self) -> str:
        if self.is_infinity():
            return "CurvePoint(O — identity)"
        return (
            f"CurvePoint(\n"
            f"  x = {self.x}\n"
            f"  y = {self.y}\n"
            f"  on curve mod {self.curve.p}\n"
            f")"
        )
