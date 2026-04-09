"""
ecc/keygen.py — Key Generation & ECDH Key Exchange
====================================================
Demonstrates Elliptic Curve Diffie-Hellman (ECDH) key agreement.

Domain parameters (shared publicly):
  curve  — the elliptic curve E over 𝔽ₚ
  G      — a generator point of prime order n
  n      — order of G  (= #E for prime-order curves)

Key generation:
  private key  d  — a random integer in [1, n-1]
  public  key  Q  — the point Q = d·G

ECDH key exchange (Alice and Bob):
  Alice sends  Q_A = d_A · G  to Bob
  Bob   sends  Q_B = d_B · G  to Alice
  Shared secret:  S = d_A · Q_B = d_B · Q_A = d_A · d_B · G

Public API
----------
generate_keypair(curve, G, n)           – Returns (private_key, public_key)
ecdh_shared_secret(private_key, public) – Returns shared CurvePoint
"""

import random
from .point import CurvePoint
from .curve import EllipticCurve


# ---------------------------------------------------------------------------
# Key Pair
# ---------------------------------------------------------------------------

def generate_keypair(
    curve: EllipticCurve,
    G: CurvePoint,
    n: int,
) -> tuple[int, CurvePoint]:
    """
    Generate an ECDH key pair.

    Parameters
    ----------
    curve : The elliptic curve.
    G     : Generator point of prime order n.
    n     : Order of G (should equal #E for a prime-order curve).

    Returns
    -------
    (d, Q) where
      d  – private key: a random integer in [1, n-1]
      Q  – public  key: Q = d·G  (a point on the curve)
    """
    d = random.randint(1, n - 1)   # private key — keep secret!
    Q = d * G                      # public key  — share openly
    return d, Q


# ---------------------------------------------------------------------------
# ECDH Shared Secret
# ---------------------------------------------------------------------------

def ecdh_shared_secret(private_key: int, peer_public_key: CurvePoint) -> CurvePoint:
    """
    Compute the ECDH shared secret.

    Each party multiplies their *own* private key by the *other's* public key:

        S = d_A · Q_B = d_A · (d_B · G) = d_B · (d_A · G) = d_B · Q_A

    The x-coordinate of S is conventionally used as the raw shared secret.

    Parameters
    ----------
    private_key     : Your own private key d  (integer).
    peer_public_key : The other party's public key Q  (CurvePoint).

    Returns
    -------
    The shared CurvePoint S.  Use S.x as the raw shared value.
    """
    return private_key * peer_public_key
