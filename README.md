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

## How It Works

### 1. Choose prime p
The user supplies or auto-generates a prime `p` of 64–128 bits. This defines 𝔽ₚ.

### 2. Generate parameters a, b
Random `a`, `b ∈ 𝔽ₚ` are sampled until the **non-singular condition** holds:
```
4a³ + 27b² ≢ 0  (mod p)
```

### 3. Count #E — Baby-step Giant-step
**Hasse's theorem** bounds the group order:
```
p + 1 - 2√p  ≤  #E  ≤  p + 1 + 2√p
```
The Baby-step Giant-step (BSGS) algorithm finds **#E exactly** within this interval in O(p^(1/4)) time.

### 4. Security checks
| Check | Threat avoided |
|---|---|
| #E is prime | Small-subgroup attacks (cofactor h = 1) |
| Not supersingular | MOV attack via Weil pairing |
| Not anomalous (#E ≠ p) | SSSA / Smart's attack |
| 64 ≤ log₂(#E) ≤ 128 | Brute-force / Pohlig-Hellman |

### 5. Generator G
Since #E is **prime**, every non-identity point is a generator.
We verify `#E · G = O` before using G.

### 6. ECDH key exchange
```
Alice: d_A (private),  Q_A = d_A · G (public)
Bob:   d_B (private),  Q_B = d_B · G (public)

Shared secret:  S = d_A · Q_B = d_B · Q_A = d_A · d_B · G
```

---

## Key Algorithms

### Double-and-Add (scalar multiplication)
```python
Q = O          # identity
R = P          # accumulator
while k > 0:
    if k & 1:
        Q = Q + R
    R = R + R  # double
    k >>= 1
return Q
```
Runs in **O(log k)** point operations (vs O(k) for naïve addition).

### Miller-Rabin (primality)
Tests `n` against `rounds` random bases; false-positive rate < 4^{-rounds}.

### Tonelli-Shanks (square root mod p)
Efficiently finds `y` such that `y² ≡ n (mod p)` — used when sampling random curve points.

---

## Requirements

- Python **3.10+** (for `int | None` type hints)
- No third-party libraries

---

## Example Session

```
╔══════════════════════════════════════════════════════════╗
║    Elliptic Curve Cryptography — Educational Demo        ║
║    Curve: y² ≡ x³ + ax + b  (mod p)                     ║
║    Target: #E is prime, 64–128 bits                      ║
╚══════════════════════════════════════════════════════════╝

──────────────────────────────────────────────────────────
  Step 1 — Choose a prime p
──────────────────────────────────────────────────────────
  Options:
    [1]  Enter your own prime p
    [2]  Auto-generate a random prime of chosen bit-length
  Choice (1 or 2): 2
  bit-length of p: 80
  ✓  Generating a random 80-bit prime — done.
  ✓  p = 1208925819614629174706189
  ✓  Bit-length: 80 bits
```
