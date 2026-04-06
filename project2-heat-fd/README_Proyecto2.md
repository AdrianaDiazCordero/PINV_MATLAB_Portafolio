# Project 2 — Backward Heat Equation via Finite Differences

Numerical solution of the **direct** and **inverse** heat equation using
explicit finite differences. Includes a rigorous study of the *inverse crime*
and noise sensitivity.

---

## Problem statement

$$\begin{cases} u_t - D\,u_{xx} = 0, & 0 < x < L,\; t > 0 \\ u(0,t) = u(L,t) = 0 \\ u(x,0) = f(x) \end{cases}$$

**Direct problem:** given $f(x)$, compute $u(x,T)$.  
**Inverse problem:** given $u(x,T)$, recover $f(x) = u(x, T_0)$ for some $T_0 < T$.

---

## Methods implemented

### Direct problem
Explicit forward-in-time finite difference scheme:

$$u_j^{n+1} = d\,u_{j-1}^n + (1-2d)\,u_j^n + d\,u_{j+1}^n, \qquad d = \frac{D\Delta t}{(\Delta x)^2}$$

The stability condition $d \le 1/2$ is enforced throughout.

### Inverse problem — Method 1: direct inversion
Runs the time-stepping backward using $A^{-1}$ (the inverse of the
forward operator). Demonstrates catastrophic instability for large $K$ steps.

### Inverse problem — Method 2: implicit backward scheme
Constructs $\tilde{A} = I + 2dI - d(\text{off-diagonals})$ and steps
backward using $\tilde{A}$ instead of $A^{-1}$. More stable for small $K$.

### Noise study
Adds uniform random perturbation $\varepsilon \in [-10^{-4}, 10^{-4}]$ to the
final state and observes amplification in the backward reconstruction.

---

## The inverse crime

> **Definition:** using the same discretization grid to *generate* the synthetic
> data and to *solve* the inverse problem artificially inflates accuracy.

This project explicitly avoids the inverse crime: the forward solution uses
grid $(M=40, N=400)$ and the inversion uses a different grid, so that any
error in the reconstruction is a genuine mathematical difficulty, not a
numerical artifact.

---

## Parameters explored

| Parameter | Values tested |
|-----------|--------------|
| Spatial nodes $M$ | 20, 40 |
| Time steps $N$ | 25, 50, 100, 400 |
| Final time $T$ | 0.1, 1.0 |
| Backward steps $K$ | 10 |
| Noise level $\varepsilon$ | $10^{-4}$ |

---

## Files

```
project2-heat-fd/
├── Proyecto2.m    # Main script
└── README.md
```

---

## How to run

```matlab
run('Proyecto2.m')
```

Outputs: surface plots of $u(x,t)$, comparison of recovered vs. exact
initial condition, noise sensitivity analysis.

---

## Mathematical connection to finance

In the Black-Scholes framework, the **backward induction** of option prices
from payoff to present value is analogous to this backward heat problem.
The instability observed here explains why naive numerical inversion of
option pricing operators amplifies bid-ask spread noise, motivating the
regularization techniques of Projects 2 and 3.
