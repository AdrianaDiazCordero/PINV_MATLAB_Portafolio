# Project 1 — Gravitational Inverse Problem

Reconstruction of a subsurface mass density distribution $\lambda(x)$
from noisy surface gravitational field measurements $F_z(t)$.

Implements and compares three strategies for selecting the optimal polynomial
degree $k$: minimum residual norm, automatic L-curve elbow detection, and
repeated random cross-validation.

---

## Problem statement

Given noisy measurements of the vertical gravitational field component at
surface positions $t_1, \ldots, t_n$, recover the mass density function
$\lambda : [0,1] \to \mathbb{R}$ satisfying

$$F_z(t_i) = \int_0^1 \frac{\lambda(x)}{((x - t_i)^2 + 1)^{3/2}} \, dx, \qquad i = 1, \ldots, n$$

The density is approximated by a polynomial of degree $k$:
$\lambda(x) \approx \sum_{j=0}^{k} c_j x^j$,
leading to the **normal equations** $A^T A\, c = A^T F_z$ where each entry of $A$ is

$$A_{ij} = \int_0^1 \frac{x^{j-1}}{((x - t_i)^2 + 1)^{3/2}} \, dx$$

computed by adaptive quadrature (`integral`). This is a **Fredholm integral equation
of the first kind** — a classical ill-posed inverse problem.

---

## Files

```
project1-gravity/
├── Proyecto1V2.m          # Dataset 1: k=1..12, L-curve, cross-validation
├── Proyecto1Ejerc2.m      # Dataset 2: k=1..5, adaptive noise, L-curve
├── calc_codoCurvaL.m      # Utility: automatic L-curve elbow detection
├── CampoGravitacional1.mat  # Gravitational field dataset 1 (required)
├── CampoGravitacional2.mat  # Gravitational field dataset 2 (required)
└── README.md
```

---

## Methods

### 1. Polynomial least squares

Solves $(A^T A)\,c = A^T F_z$ for polynomial degrees $k = 1, \ldots, k_{\max}$.
The condition number of $A^T A$ is monitored: it grows rapidly with $k$,
illustrating the ill-posed character of the problem.

### 2. L-curve method

Plots $\log \|Ac - F_z\|$ vs. $\log \|c\|$ in log-log scale.
The optimal $k^*$ is found by two approaches:
- **Visual inspection** of the corner of the L
- **Automatic elbow detection** via `calc_codoCurvaL` (see below)

### 3. Automatic L-curve elbow — `calc_codoCurvaL.m`

Computes the geometric elbow of the L-curve by finding the point with
maximum perpendicular distance to the straight line connecting the two
endpoints of the curve (in log-log space):

$$d_i = \frac{|(y_2 - y_1)\,x_i - (x_2 - x_1)\,y_i + x_2 y_1 - y_2 x_1|}{\sqrt{(y_2 - y_1)^2 + (x_2 - x_1)^2}}$$

where $x_i = \log \|r_i\|$, $y_i = \log \|c_i\|$.

```matlab
function codo_opt = calc_codoCurvaL(norma_res, norma_sol)
```

Returns the index of the optimal $k$ in the input arrays. Works for any
monotone L-shaped curve — reusable across projects.

### 4. Repeated random cross-validation (Dataset 1)

Runs `rep = 20` random 80/20 train-test splits. For each split and each $k$:
- Builds $A_{\text{train}}$ and $A_{\text{test}}$ by numerical integration
- Fits coefficients on training data: $c = (A_{\text{train}}^T A_{\text{train}})^{-1} A_{\text{train}}^T F_{\text{train}}$
- Evaluates prediction error on test set: $\|A_{\text{test}}\,c - F_{\text{test}}\|$

Reports the **median** optimal $k$ across repetitions (both by minimum residual
and by L-curve elbow), making the selection robust to random split variance.

---

## Noise model

| Dataset | Noise level | Rationale |
|---------|-------------|-----------|
| Dataset 1 | $\varepsilon = 0.001$ (absolute) | Small fixed perturbation |
| Dataset 2 | $\varepsilon = 0.05 \cdot \text{range}(F_z)$ | Proportional to data scale |

Noise is uniform: $F_z^\delta = F_z + \varepsilon_i$, $\varepsilon_i \sim U(-\varepsilon, \varepsilon)$.
Dataset 2 uses `rng(42)` for reproducibility.

---

## Key results

| Dataset | $k_{\max}$ explored | Optimal $k$ (L-curve) | Optimal $k$ (cross-val median) |
|---------|--------------------|-----------------------|-------------------------------|
| CampoGravitacional1 | 1–12 | reported at runtime | reported at runtime |
| CampoGravitacional2 | 1–5 | reported at runtime | — |

> For Dataset 2 the range is deliberately limited to $k \le 5$ because the condition
> number of $A^T A$ grows beyond $10^{10}$ for higher degrees, making the
> least-squares solution numerically unreliable.

---

## How to run

```matlab
% Place all .mat files in the same folder as the scripts, then:

% Dataset 1 (full analysis: L-curve + cross-validation)
run('Proyecto1V2.m')

% Dataset 2
run('Proyecto1Ejerc2.m')
```

`calc_codoCurvaL.m` is called automatically — no separate execution needed.

**Generated output (Dataset 1):**
1. Density reconstructions for $k = 1, \ldots, 12$ (with and without noise)
2. Goodness-of-fit comparison for the last $k$
3. L-curve in log-log scale (noisy data)
4. L-curve in log-log scale (clean data)
5. Console: automatic elbow $k^*$ for both datasets
6. Console: cross-validation median $k^*$ over 20 random splits

---

## Mathematical connection to finance

The least-squares framework with regularization via degree selection appears
directly in quantitative finance in two settings:

- **Yield curve fitting** (Nelson-Siegel, Svensson): market bond prices are
  used to recover the term structure of interest rates — another Fredholm-type
  inverse problem over financial observations.
- **Implied volatility surface calibration**: fitting a smooth parametric
  surface $\sigma(K, T)$ to discrete option quotes requires balancing data
  fidelity against overfitting, exactly as the L-curve balances residual norm
  against solution norm here.

The automatic elbow detection in `calc_codoCurvaL` is directly analogous to
information criteria (AIC, BIC) used in model selection for financial time series.
