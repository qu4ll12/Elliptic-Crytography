"""
ecc — A pure-Python Elliptic Curve Cryptography library.

Curve form:  y² ≡ x³ + ax + b  (mod p)
Target:      64–128 bit prime curve order (#E)
"""

from .utils import is_prime, generate_prime
from .curve import EllipticCurve
from .point import CurvePoint
from .order import curve_order, is_secure_curve
from .keygen import generate_keypair, ecdh_shared_secret

__all__ = [
    "EllipticCurve",
    "CurvePoint",
    "generate_prime",
    "is_prime",
    "curve_order",
    "is_secure_curve",
    "generate_keypair",
    "ecdh_shared_secret",
]
