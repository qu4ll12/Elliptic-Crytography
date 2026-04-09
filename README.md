# Elliptic Curve Cryptography — Python Implementation

A **pure-Python** educational implementation of elliptic curve cryptography (ECC)
over a prime finite field 𝔽ₚ, targeting a **prime curve order #E** between **64 and 128 bits**.

---

## Curve Form

```
y² ≡ x³ + a·x + b  (mod p)
```

---

## Project Structure

```
Elliptic Cryptography/
│
├── ecc/
│   ├── __init__.py   — Public API re-exports
│   ├── utils.py      — Miller-Rabin, modular inverse, Tonelli-Shanks
│   ├── curve.py      — EllipticCurve class + generate_curve(p)
│   ├── point.py      — CurvePoint with Double-and-Add scalar multiplication
│   ├── order.py      — Baby-step Giant-step order counting + security checks
│   └── keygen.py     — Key pair generation & ECDH demo
│
├── main.py           — Interactive CLI
└── README.md         — This file
```

---

## Quick Start

```bash
# No dependencies required — pure Python 3.10+
python main.py
```

---