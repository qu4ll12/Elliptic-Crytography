"""
main.py — Interactive ECC Demo
================================
A guided command-line tool for exploring elliptic curve cryptography
over a prime field 𝔽ₚ with a prime curve order #E in [64, 128] bits.

Run with:
    python main.py
"""

import sys
import time

from ecc.utils     import is_prime, generate_prime
from ecc.curve     import EllipticCurve, generate_curve
from ecc.point     import CurvePoint
from ecc.order     import curve_order, is_secure_curve, find_generator
from ecc.keygen    import generate_keypair, ecdh_shared_secret


# ---------------------------------------------------------------------------
# Terminal colours (works on Linux/macOS; degrades gracefully on Windows)
# ---------------------------------------------------------------------------

RESET  = "\033[0m"
BOLD   = "\033[1m"
CYAN   = "\033[96m"
GREEN  = "\033[92m"
YELLOW = "\033[93m"
RED    = "\033[91m"
DIM    = "\033[2m"

def c(text: str, colour: str) -> str:
    return f"{colour}{text}{RESET}"

def section(title: str):
    """Print a visual section divider."""
    width = 58
    print()
    print(c("─" * width, CYAN))
    print(c(f"  {title}", BOLD + CYAN))
    print(c("─" * width, CYAN))


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def prompt_int(label: str, lo: int | None = None, hi: int | None = None) -> int:
    """Ask the user for an integer, with optional bounds checking."""
    while True:
        raw = input(f"  {label}: ").strip()
        if not raw:
            continue
        try:
            val = int(raw)
        except ValueError:
            print(c("  ✗  Please enter a valid integer.", RED))
            continue
        if lo is not None and val < lo:
            print(c(f"  ✗  Value must be ≥ {lo}.", RED))
            continue
        if hi is not None and val > hi:
            print(c(f"  ✗  Value must be ≤ {hi}.", RED))
            continue
        return val


def tick(msg: str):
    """Print a success line."""
    print(c(f"  ✓  {msg}", GREEN))


def info(msg: str):
    """Print an informational line."""
    print(f"  {c('›', CYAN)} {msg}")


def warn(msg: str):
    """Print a warning line."""
    print(c(f"  ⚠  {msg}", YELLOW))


def fail(msg: str):
    """Print an error line."""
    print(c(f"  ✗  {msg}", RED))


def spinner(msg: str, func, *args, **kwargs):
    """
    Run `func(*args, **kwargs)` while printing a simple progress indicator,
    then return the result.
    """
    frames = ["⠋", "⠙", "⠹", "⠸", "⠼", "⠴", "⠦", "⠧", "⠇", "⠏"]
    import threading

    result   = [None]
    exc      = [None]
    done     = [False]

    def worker():
        try:
            result[0] = func(*args, **kwargs)
        except Exception as e:
            exc[0] = e
        finally:
            done[0] = True

    t = threading.Thread(target=worker)
    t.start()

    i = 0
    while not done[0]:
        print(f"\r  {c(frames[i % len(frames)], CYAN)}  {msg} …", end="", flush=True)
        time.sleep(0.1)
        i += 1

    print(f"\r  {c('✓', GREEN)}  {msg} — done.{' ' * 20}")
    t.join()

    if exc[0]:
        raise exc[0]
    return result[0]


# ---------------------------------------------------------------------------
# Step 1 — Choose prime p
# ---------------------------------------------------------------------------

def step_choose_p() -> int:
    section("Step 1 — Choose a prime p")
    print()
    print("  The prime p defines the finite field 𝔽ₚ over which the curve lives.")
    print("  The larger p is, the harder the discrete-log problem becomes.")
    print()
    print("  Options:")
    print(c("    [1]", CYAN) + "  Enter your own prime p")
    print(c("    [2]", CYAN) + "  Auto-generate a random prime of chosen bit-length")
    print()

    choice = prompt_int("Choice (1 or 2)", lo=1, hi=2)

    if choice == 1:
        print()
        info("Enter a prime p (suggested: 18–39 decimal digits for 64–128 bit security)")
        p = prompt_int("p")
        if not is_prime(p):
            fail(f"{p} is NOT prime. Please restart and enter a prime.")
            sys.exit(1)
        tick(f"p = {p}  ({p.bit_length()} bits) — prime ✓")
        return p

    else:
        print()
        info("Choose a bit-length for p in the range [64, 128].")
        info("For the curve order #E to also be ≥64 bits, p should be ≥ 64 bits.")
        bits = prompt_int("bit-length of p", lo=64, hi=128)
        p = spinner(f"Generating a random {bits}-bit prime", generate_prime, bits)
        tick(f"p = {p}")
        tick(f"Bit-length: {p.bit_length()} bits")
        return p


# ---------------------------------------------------------------------------
# Step 2 — Generate curve parameters a, b
# ---------------------------------------------------------------------------

def step_generate_curve(p: int) -> EllipticCurve:
    section("Step 2 — Generate curve parameters a, b")
    print()
    print("  The short Weierstrass equation is:")
    print(c("      y² ≡ x³ + a·x + b  (mod p)", BOLD))
    print()
    print("  We need:  4a³ + 27b² ≢ 0 (mod p)  — the non-singular (discriminant) condition.")
    print("  Parameters a and b are chosen uniformly at random and validated.")
    print()

    curve = spinner("Generating non-singular curve parameters", generate_curve, p)

    print()
    print(c("  Curve parameters:", BOLD))
    info(f"a = {curve.a}")
    info(f"b = {curve.b}")
    info(f"Δ = 4a³ + 27b² = {curve._discriminant()}  (mod p)  ≠ 0 ✓")
    print()
    print(c("  Curve equation:", BOLD))
    print(f"    y² ≡ x³  +  {curve.a}·x  +  {curve.b}  (mod {curve.p})")

    return curve


# ---------------------------------------------------------------------------
# Step 3 — Count curve order #E via Baby-step Giant-step
# ---------------------------------------------------------------------------

def step_count_order(curve: EllipticCurve) -> int:
    section("Step 3 — Count curve order  #E  (Baby-step Giant-step)")
    print()
    print("  Hasse's theorem guarantees:")
    print(c("      p + 1 - 2√p  ≤  #E  ≤  p + 1 + 2√p", BOLD))
    print()
    print("  We use the Baby-step Giant-step (BSGS) algorithm to find #E")
    print("  exactly within this interval without enumerating all points.")
    print()

    order = spinner("Computing #E via Baby-step Giant-step", curve_order, curve)

    print()
    info(f"#E = {order}")
    info(f"Bit-length of #E: {order.bit_length()} bits")

    return order


# ---------------------------------------------------------------------------
# Step 4 — Security check
# ---------------------------------------------------------------------------

def step_security_check(curve: EllipticCurve, order: int) -> bool:
    section("Step 4 — Security checks")
    print()

    checks = is_secure_curve(curve, order)

    labels = {
        "prime_order"      : "#E is prime            (cofactor h = 1, no small-subgroup attacks)",
        "not_singular"     : "Curve is non-singular  (discriminant ≠ 0)",
        "not_supersingular": "Not supersingular       (safe from MOV/Weil-pairing attack)",
        "bit_range"        : f"#E bit-length in [64, 128]  (currently {checks['order_bits']} bits)",
        "not_anomalous"    : "Not anomalous          (#E ≠ p, safe from SSSA/Smart attack)",
    }

    all_pass = True
    for key, label in labels.items():
        passed = checks[key]
        if passed:
            tick(label)
        else:
            fail(label)
            all_pass = False

    print()
    if all_pass:
        print(c("  ✓  All security checks passed — this is a cryptographically sound curve.", GREEN + BOLD))
    else:
        warn("One or more security checks FAILED.")
        warn("You should regenerate the curve with different parameters.")

    return all_pass


# ---------------------------------------------------------------------------
# Step 5 — Generator point G
# ---------------------------------------------------------------------------

def step_find_generator(curve: EllipticCurve, order: int) -> CurvePoint:
    section("Step 5 — Generator point G")
    print()
    print("  Since #E is prime, every non-identity point on the curve is a")
    print("  generator of the full group.  We sample one at random and verify:")
    print(c("      #E · G = O  (point at infinity)", BOLD))
    print()

    G = spinner("Finding a generator point G", find_generator, curve, order)

    print()
    info(f"G.x = {G.x}")
    info(f"G.y = {G.y}")
    tick(f"Verified:  #E · G = O ✓")

    return G


# ---------------------------------------------------------------------------
# Step 6 — ECDH key exchange demo
# ---------------------------------------------------------------------------

def step_ecdh(curve: EllipticCurve, G: CurvePoint, order: int):
    section("Step 6 — ECDH Key Exchange Demo")
    print()
    print("  Elliptic Curve Diffie-Hellman (ECDH) lets two parties establish")
    print("  a shared secret over an insecure channel without ever transmitting")
    print("  their private keys.")
    print()
    print("  Protocol:")
    print("    Alice: picks private key d_A,  computes Q_A = d_A · G")
    print("    Bob:   picks private key d_B,  computes Q_B = d_B · G")
    print("    Alice: computes S = d_A · Q_B")
    print("    Bob:   computes S = d_B · Q_A")
    print("    Both arrive at the same S = d_A · d_B · G ✓")
    print()

    # Alice
    d_A, Q_A = generate_keypair(curve, G, order)
    print(c("  Alice:", BOLD))
    info(f"Private key (d_A) = {d_A}")
    info(f"Public  key (Q_A.x) = {Q_A.x}")
    info(f"Public  key (Q_A.y) = {Q_A.y}")
    print()

    # Bob
    d_B, Q_B = generate_keypair(curve, G, order)
    print(c("  Bob:", BOLD))
    info(f"Private key (d_B) = {d_B}")
    info(f"Public  key (Q_B.x) = {Q_B.x}")
    info(f"Public  key (Q_B.y) = {Q_B.y}")
    print()

    # Compute shared secret
    S_A = ecdh_shared_secret(d_A, Q_B)   # Alice computes d_A · Q_B
    S_B = ecdh_shared_secret(d_B, Q_A)   # Bob   computes d_B · Q_A

    print(c("  Shared secret:", BOLD))
    info(f"Alice computes  S_A.x = {S_A.x}")
    info(f"Bob   computes  S_B.x = {S_B.x}")
    print()

    if S_A == S_B:
        tick("S_A == S_B — shared secret matches! ✓")
        info("The x-coordinate of S can now be used as a symmetric key (e.g., for AES).")
    else:
        fail("Shared secrets do NOT match — something went wrong!")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def banner():
    print()
    print(c("╔══════════════════════════════════════════════════════════╗", CYAN))
    print(c("║", CYAN) + c("    Elliptic Curve Cryptography — Educational Demo        ", BOLD) + c("║", CYAN))
    print(c("║", CYAN) + c("    Curve: y² ≡ x³ + ax + b  (mod p)                     ", DIM)  + c(" ║", CYAN))
    print(c("║", CYAN) + c("    Target: #E is prime, 64–128 bits                      ", DIM)  + c("║", CYAN))
    print(c("╚══════════════════════════════════════════════════════════╝", CYAN))


def main():
    banner()

    # Keep looping until we find a curve that passes all security checks
    while True:
        p     = step_choose_p()
        curve = step_generate_curve(p)
        order = step_count_order(curve)
        ok    = step_security_check(curve, order)

        if ok:
            break

        print()
        warn("Curve did not pass all checks. Let's try again with different parameters.")
        again = input(c("\n  Try again? [Y/n]: ", YELLOW)).strip().lower()
        if again == "n":
            print()
            info("Exiting. Goodbye!")
            sys.exit(0)

    G = step_find_generator(curve, order)
    step_ecdh(curve, G, order)

    # Final summary
    section("Summary")
    print()
    print(c("  Domain Parameters (share these publicly):", BOLD))
    print()
    print(f"    Curve equation : y² ≡ x³ + {curve.a}·x + {curve.b}  (mod {curve.p})")
    print(f"    Prime  p       : {curve.p}")
    print(f"    Curve order #E : {order}  ({order.bit_length()} bits, prime ✓)")
    print(f"    Generator G.x  : {G.x}")
    print(f"    Generator G.y  : {G.y}")
    print()
    print(c("  Done. Thank you for exploring ECC!", CYAN + BOLD))
    print()


if __name__ == "__main__":
    main()
