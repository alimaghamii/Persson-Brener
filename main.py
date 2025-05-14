"""
This code computes the effective gamma value for a range of nu values
and plots the results for different power-law exponents (n).
This code is based on the MPL viscoelastic model, which is a theoretical
model used to describe the behavior of viscoelastic materials.

Based on the paper:
Maghami, Ali, et al. "Bulk and fracture process zone contribution to the rate-dependent adhesion amplification in viscoelastic broad-band materials." 
Journal of the Mechanics and Physics of Solids, vol. 193, 2024, pp. 105844. Elsevier.

Plot  \hat{Γ}_eff  vs  \hat{ν}  for the MPL visco‑elastic model
using Eq. (B.1) with the closed‑form integral (B.2–B.3), for multiple n.

Author: <Ali Maghami>
Date: 2025-05-14
License: MIT License
"""

import math
import numpy as np
import mpmath as mp
import matplotlib.pyplot as plt

# ----------------------------------------------------------
# model parameters
k = 0.10
n_list = [0.2,0.4,0.6, 0.8, 1.6]  # multiple values of n
nu_min, nu_max = 1e-2, 1e8
num_points = 200
# ----------------------------------------------------------

# ---------------------------------------------------------------------------
# helpers
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

def gamma_eff(n, k, nu, tol=1e-10, max_iter=200):
    Gamma = 1.0
    for _ in range(max_iter):
        I = I_value(n, nu, Gamma)
        Gamma_next = 1.0 / (1.0 - (1.0 - k) * I)
        if abs(Gamma_next - Gamma) < tol:
            return float(Gamma_next)
        Gamma = Gamma_next
    raise RuntimeError(f"Convergence failed at ν̂ = {nu:g}")

# ---------------------------------------------------------------------------
# plot
# ---------------------------------------------------------------------------
nu_vals = np.logspace(math.log10(nu_min), math.log10(nu_max), num=num_points)

plt.figure(figsize=(8, 5))
for n in n_list:
    gamma_vals = np.array([gamma_eff(n, k, nu) for nu in nu_vals])
    plt.loglog(nu_vals, gamma_vals, label=fr"$n = {n}$")

plt.xlabel(r"$\hat{\nu}$", fontsize=12)
plt.ylabel(r"$\hat{\Gamma}_{\mathrm{eff}}$", fontsize=12)
plt.title(r"$\hat{\Gamma}_{\mathrm{eff}}$ vs. $\hat{\nu}$"
          + f"   ($k={k}$)", fontsize=13, pad=12)
plt.grid(True, which="both", ls=":")
plt.legend(title="Power-law exponent", fontsize=10)
plt.tight_layout()
plt.show()
