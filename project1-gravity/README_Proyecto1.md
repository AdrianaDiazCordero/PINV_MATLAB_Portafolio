# Project 1 — Gravitational Inverse Problem

Reconstruction of a subsurface mass density distribution $\lambda(x)$
from surface gravitational field measurements $F_z(t)$.

---

## Problem statement

Given noisy measurements of the vertical gravitational field component at
surface positions $t_1, \ldots, t_n$, recover the mass density function
$\lambda : [0,1] \to \mathbb{R}$ satisfying

$$F_z(t_i) = \int_0^1 \frac{x^{j-1}}{((x - t_i)^2 + 1)^{3/2}} \, dx, \qquad i = 1, \ldots, n$$

The density is approximated by a polynomial of degree $k$:
$\lambda(x) \approx \sum_{j=0}^{k} c_j x^j$,
leading to the linear system $AC = F_z$ where each entry of $A$ is computed
by numerical integration.

This is a **Fredholm integral equation of the first kind** — a classic
ill-posed inverse problem.

---

## Methods implemented

### 1. Polynomial least squares
Solves $(A^T A) c = A^T F_z$ for increasing polynomial degree $k = 1, \ldots, 12$.
Condition number of $A^T A$ is monitored as a function of $k$.

### 2. L-curve method
Plots $\log \|Ac - F_z\|$ vs. $\log \|c\|$ in log-log scale.
The corner of the L-curve identifies the optimal degree $k^*$ that
balances data fidelity against solution norm (regularization effect).

### 3. Cross-validation (80/20 split)
The dataset is split: 80% for training, 20% for testing.
The optimal $k$ is selected by minimizing the prediction error on the
held-out set, providing an independent stability criterion.

---

## Key results

| Dataset | Optimal k (L-curve) | Optimal k (cross-val) |
|---------|--------------------|-----------------------|
| CampoGravitacional1 | 4–5 | 4–5 |
| CampoGravitacional2 | varies | varies |

> **Noise sensitivity:** even small perturbations ($\varepsilon = 0.001$)
> produce significant changes in $\lambda(x)$ for high $k$, confirming the
> ill-posed character of the problem.

---

## Files

```
project1-gravity/
├── Proyecto1V2.m          # Main script
├── CampoGravitacional1.mat  # Dataset 1 (required)
├── CampoGravitacional2.mat  # Dataset 2 (required)
└── README.md
```

---

## How to run

```matlab
% Place .mat files in the same folder as the script, then:
run('Proyecto1V2.m')
```

The script will generate:
- Density reconstructions for both datasets (with and without noise)
- L-curves in log-log scale
- Cross-validation error curves

---

## Mathematical connection to finance

The same least-squares framework with regularization via polynomial degree
selection appears in **yield curve fitting** (e.g., Nelson-Siegel models),
where market bond prices are used to recover the term structure of interest rates —
another Fredholm-type inverse problem over financial data.
