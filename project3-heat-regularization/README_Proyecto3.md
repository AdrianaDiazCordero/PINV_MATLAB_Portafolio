# Project 3 — Regularization Methods for the Backward Heat Equation

Spectral formulation of the backward heat inverse problem with a systematic
comparison of four regularization families. This project is the basis of a
scientific article submitted to *Ciencias Matemáticas* (2026).

---

## Problem statement

Recover the initial temperature $f(x) = u(x,0)$ of a rod from observations
$g(x) = u(x,T)$ at a later time $T > 0$. The problem is discretized via
a **spectral method** (truncated Fourier sine series), yielding the linear system

$$AF = G, \qquad \kappa(A) \approx 10^{18}$$

The condition number makes direct inversion computationally equivalent to
division by zero. All four methods below address this instability.

---

## Methods implemented

### 1. Tikhonov regularization (classic)
Replaces $AF = G$ with the penalized minimization

$$F_\alpha = \arg\min_F \left\{ \|AF - G\|^2 + \alpha\|F\|^2 \right\}$$

**Parameter selection:** Morozov discrepancy principle — finds the largest $\alpha$
such that $\|AF_\alpha - G\| \in [\delta, 1.2\delta]$ with $\delta = 0.01\|G\|$.  
**Result:** condition number drops from $10^{18}$ to $\approx 10^2$.

### 2. Iterated Tikhonov
Reduces the smoothing bias of classic Tikhonov by iterating:

$$(A^TA + \alpha I)\,\varphi_n = A^T G + \alpha\,\varphi_{n-1}$$

Three iterations (Tikhonov 1, 2, 3) progressively recover sharper features
with minimal additional computational cost (matrix already factored).

### 3. Landweber iteration
Gradient-descent minimization of the residual:

$$f_{k+1} = f_k - \mu\, A^T(Af_k - G)$$

Two relaxation parameters tested: $\mu = 0.95/\|A^TA\|$ and $\mu = 0.5/\|A^TA\|$.
Demonstrates **semi-convergence**: $K$ acts as regularization parameter;
too many iterations over-fits the noise.

### 4. Spectral regularization: TSVD and WSVD

**TSVD** (Truncated SVD): discards components beyond index $k$:
$$F_k = \sum_{j=1}^{k} \frac{1}{\sigma_j}\langle u_j, G\rangle v_j$$

**WSVD** (Weighted SVD): caps the denominator at $\sigma_k$ instead of truncating:
$$c_j^\alpha = \begin{cases} 1/\sigma_j & \text{if } \sigma_j > \sigma_k \\ 1/\sigma_k & \text{otherwise} \end{cases}$$

Optimal truncation: $k = 5$ or $6$ (for $T = 0.4$, matrix rank $r = 9$).

---

## Key results

| Method | Reconstruction quality | Noise robustness | Requires matrix inverse |
|--------|----------------------|-----------------|------------------------|
| Tikhonov classic | Medium | **High** | Yes |
| Tikhonov iterated | High | **High** | Yes |
| Landweber | High | Low | No |
| TSVD ($k$ optimal) | **Very high** | Medium | No |
| WSVD ($k$ optimal) | High | Low | No |

---

## Files

```
project3-heat-regularization/
├── CopiaProyecto3.m    # Main script (T=0.4 and T=1 experiments)
└── README.md
```

> **Note:** the `.mlx` live script version contains additional formatted
> output and intermediate plots not included in the `.m` file.

---

## How to run

```matlab
% Set T = 0.4 or T = 1 at line 8, then:
run('CopiaProyecto3.m')
```

Generated figures (in order):
1. Direct problem: spectral vs. finite differences comparison
2. Tikhonov classic reconstruction
3. Tikhonov iterated (iterations 0–3 overlaid)
4. Landweber ($\mu = 0.95$), clean data
5. Tikhonov with noisy data
6. Tikhonov iterated with noisy data
7. Landweber ($\mu = 0.95$), noisy data
8. Singular values (semilog scale)
9–17. TSVD and WSVD for $k = 1, \ldots, r$ (clean data)
18–26. TSVD and WSVD for $k = 1, \ldots, r$ (noisy data)

---

## Associated publication

> Díaz Cordero, A. (2026). *Métodos de regularización para el problema inverso
> de la ecuación del calor: un análisis comparativo*.
> Submitted to *Ciencias Matemáticas*.

---

## Mathematical connection to finance

The Black-Scholes PDE becomes the heat equation under the change of variables
$x = \ln S$, $\tau = T - t$, $D = \sigma^2/2$.
**Implied volatility calibration** — recovering $\sigma(S,t)$ from observed
option prices — is an inverse problem of exactly this type.
The Tikhonov and TSVD methods implemented here are used in practice by
quantitative analysts to stabilize volatility surface calibration.
