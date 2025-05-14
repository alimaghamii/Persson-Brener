"""
This code computes the effective gamma value for a range of unloading rate (r_u)
values and plots the results for different power-law exponents (n).
The viscoelastic response is modeled via the MPL framework from:

Maghami, A., et al., JMPS, Vol. 193, 2024, 105844.

Author: Ali Maghami
Date: 2025-05-14
License: MIT License
"""

import math
import numpy as np
import mpmath as mp
import matplotlib.pyplot as plt

# ----------------------------------------------------------
# Model parameters
k = 0.10
n_list = [0.2, 0.4, 0.6, 0.8, 1.6]  # MPL exponents
ru_min, ru_max = 1e-2, 1e10         # Unloading rate range
num_points = 200                   # Grid resolution
# ----------------------------------------------------------

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def _gamma_product(b_list):
    prod = 1
    for b in b_list:
        prod *= mp.gamma(b)
    return prod

def hyper_reg(a_list, b_list, z):
    return mp.hyper(a_list, b_list, z) / _gamma_product(b_list)

def I_value(n, nu, Gamma):
    pre = (2 ** (-3 - 2 * n)) * (math.pi ** (-1.5 - n)) / ((n - 1) * nu)
    x_sq = (Gamma / (4 * math.pi * nu)) ** 2
    z = -x_sq

    hyper1 = mp.hyper([-0.5], [0.5 - n / 2, 1 - n / 2], z)
    term1 = -(4 ** (1 + n)) * (math.pi ** (0.5 + n)) * (
        Gamma - 2 * (n - 1) * math.pi * nu * hyper1
    )

    hyper2 = hyper_reg([(n - 1) / 2], [0.5, (2 + n) / 2], z)
    hyper3 = hyper_reg([n / 2], [1.5, (3 + n) / 2], z)

    term2 = (
        2
        * math.pi
        * (Gamma / nu) ** n
        * (
            4 * math.pi * nu * mp.gamma(1 - n / 2) * hyper2
            + Gamma * mp.gamma(1.5 - n / 2) * hyper3
        )
    )

    return pre * (term1 + term2)

def gamma_eff_PB(n, k, nu, tol=1e-10, max_iter=200):
    Gamma = 1.0
    for _ in range(max_iter):
        I = I_value(n, nu, Gamma)
        Gamma_next = 1.0 / (1.0 - (1.0 - k) * I)
        if abs(Gamma_next - Gamma) < tol:
            return float(Gamma_next)
        Gamma = Gamma_next
    raise RuntimeError(f"Convergence failed at ν̂ = {nu:g}")

# ---------------------------------------------------------------------------
# Plot
# ---------------------------------------------------------------------------
ru_vals = np.logspace(np.log10(ru_min), np.log10(ru_max), num_points)

plt.figure(figsize=(8, 5))
for n in n_list:
    # Convert r_u to ν̂
    # Conversion constants (from paper)
    C1 = 2.887
    C2 = 3.24 * np.pi**(2/3)
    nu_vals = C1 * ((ru_vals / C2) ** 1.171)
    gamma_vals = np.array([gamma_eff_PB(n, k, nu) for nu in nu_vals])
    plt.loglog(ru_vals, gamma_vals, label=fr"$n = {n}$")

plt.xlabel(r"Unloading Rate $r_u$", fontsize=12)
plt.ylabel(r"$\hat{\Gamma}_{\mathrm{eff}}$", fontsize=12)
plt.title(r"$\hat{\Gamma}_{\mathrm{eff}}$ vs. $r_u$" + f"   ($k={k}$)", fontsize=13, pad=12)
plt.grid(True, which="both", ls=":")
plt.legend(title="Power-law exponent", fontsize=10)
plt.tight_layout()
plt.savefig("gamma_vs_ru_plot.png", dpi=300)
plt.show()
