# Project 4 — Gappy Data Recovery via Iterative SVD

Matrix completion for functions and images with randomly missing entries,
using low-rank approximation and iterative singular value decomposition.
Includes a comparative spectral analysis of initialization strategies and
a systematic evaluation of reconstruction quality across multiple settings.

---

## Problem statement

Given a matrix $A \in \mathbb{R}^{n \times n}$ where a fraction $p$ of entries
are missing (set to NaN), recover the complete matrix by exploiting the
assumption that $A$ has **low numerical rank**.

Two data types are studied:

- **Synthetic function:** a complex 3D function $f(x,y,z)$ normalized globally
  to $[-1, 1]$, sampled on a $100 \times 100$ grid for the slice $z = 0.5$.
- **Real grayscale images:** two photographs with $p = 0.40$ of pixels removed
  at random.

---

## Files

```
project4-gappy-svd/
├── Proyecto4_AdrianaDiaz.m   # Main script (complete pipeline)
├── mansion.jpg                 # Test image 1 (required)
├── palm.jpg                    # Test image 2 (required)
└── README.md
```

---

## The synthetic function

The function under study is:

$$f(x,y,z) = x^2(\sin(5\pi y + 3\log(x^3+y^2+z+\pi^2))-1)^2 - 4x^2 y^3 (1-z)^{3/2} + (x+z-1)(2y-z)\cos(30(x+z))\log(6+x^2y^2+z^3)$$

evaluated on $[0,1]^3$, then **globally normalized** to $[-1,1]$:

$$F(x,y,z) = 2 \cdot \frac{f(x,y,z) - f_{\min}}{f_{\max} - f_{\min}} - 1$$

where $f_{\min}$ and $f_{\max}$ are computed over the full 3D volume on a
$100 \times 100 \times 100$ grid, then the 3D arrays are cleared from memory.
Experiments use the $z = 0.5$ slice (matrix $A_2$).

---

## Methods

### Initialization strategies

| Strategy | Implementation | Notes |
|----------|---------------|-------|
| Zero-fill | `A2_gappy(missing) = 0` | Simple baseline |
| MAKIMA spline | `fillmissing(..., "makima", 2, "EndValues","nearest")` | Along columns; best starting point |
| White-fill | `fillmissing(..., 'constant', 1)` | For image experiments |

### Spectral analysis of initialization

Before iterating, the script computes and plots the **normalized singular values**
$\sigma_j / \sigma_1$ of each initialized matrix in semilog scale.
This comparison reveals how well each initialization preserves the low-rank
structure of the original data — a key diagnostic for convergence speed and
reconstruction quality.

The number of numerically significant modes is estimated as:

$$m_{op} = \#\{j : \sigma_j > 10^{-10}\}$$

Both zero-fill and MAKIMA initializations report their own $m_{op}$, which
differ in general.

### Single-step low-rank projection (no iteration)

For the synthetic function with zero-fill initialization, the script computes
a one-shot rank-$m$ approximation by projecting the rank-$m$ SVD directly
onto the missing positions:

```matlab
for i = 1:m_op
    Ei = U_rec(:,i) * V_rec(:,i)';
    Vector_vals_gappy = Vector_vals_gappy + singularValues(i)*Ei;
    A2_recons2(missing) = Vector_vals_gappy(missing);
end
```

This serves as a baseline to quantify how much the iterative refinement
actually gains over a single projection step.

### Core algorithm: iterative SVD (Gappy-SVD)

```
Given: A_init (initialized matrix), missing positions, m (modes), s (iterations)

for k = 1 to s:
    [U, S, V] = svd(A_current, "econ")
    Ak = U[:,1:m] * S[1:m,1:m] * V[:,1:m]'    % rank-m approximation
    A_current[missing] = Ak[missing]             % update only missing entries
end
```

Known entries are **never modified**. Missing entries converge toward the
low-rank completion of the matrix.

### Mode selection

For the **synthetic function**: $m = 14$ modes (fixed), based on visual
inspection of the singular value spectrum.

For **images**: $m = \lfloor 0.1 \cdot m_{op} \rfloor$ (10% of the numerical
rank of the initialized matrix), computed separately for each image and
each initialization.

---

## Experiments

### Synthetic function ($p = 0.40$, $z = 0.5$ slice, $m = 14$, $s = 20$)

| Method | Metric |
|--------|--------|
| Zero-fill, no iteration | Max error + RMSE reported |
| Zero-fill + Gappy-SVD ($m=14$, $s=20$) | Max error + RMSE reported |
| MAKIMA only (no iteration) | Max error + RMSE reported |
| MAKIMA + Gappy-SVD ($m=14$, $s=20$) | Max error + RMSE reported |

### Image 1 — `mansion.jpg` ($p = 0.40$, $s = 50$, $m = \lfloor 0.1 \cdot m_{op} \rfloor$)

| Initialization | Metric reported |
|---------------|-----------------|
| MAKIMA (pre-iteration baseline) | Max error + RMSE |
| MAKIMA + Gappy-SVD | Max error + RMSE |
| White-fill (pre-iteration baseline) | Max error + RMSE |
| White-fill + Gappy-SVD | Max error + RMSE |

### Image 2 — `palm.jpg` ($p = 0.40$, $s = 50$, $m = \lfloor 0.1 \cdot m_{op} \rfloor$)

| Initialization | Metric reported |
|---------------|-----------------|
| White-fill (pre-iteration baseline) | Max error + RMSE |
| White-fill + Gappy-SVD | Max error + RMSE |
| MAKIMA (pre-iteration baseline) | Max error + RMSE |
| MAKIMA + Gappy-SVD | Max error + RMSE |

> Both initialization strategies are tested on both images, enabling a
> direct comparison of MAKIMA vs. white-fill across different image content.

---

## Error metrics

Two metrics are reported for every reconstruction:

$$\text{Max error} = \max_{i,j} |A_{ij} - \hat{A}_{ij}|$$

$$\text{RMSE} = \frac{1}{n} \|A - \hat{A}\|_F$$

where $n$ is the number of rows and $\|\cdot\|_F$ is the Frobenius norm.

---

## How to run

```matlab
% Place image files in the same folder as the script, then:
run('Proyecto4_AdrianaDiaz.m')
```

The script prints progress markers to the console (`disp` calls) as it moves
through each experiment, so you can track execution in long runs.

**Generated figures (in order):**
1. Original function (unnormalized), $z = 0.5$
2. Normalized function, $z = 0.5$
3. Gappy function (40% missing)
4. MAKIMA initialization
5. Singular values — zero-fill initialization (normalized, semilog)
6. Singular values — MAKIMA initialization (normalized, semilog)
7. Single-step SVD reconstruction ($m = 14$, no iteration)
8. Iterative SVD from zero-fill ($m = 14$, $s = 20$)
9. Iterative SVD from MAKIMA ($m = 14$, $s = 20$)
10. `mansion.jpg` — original
11. `mansion.jpg` — 40% missing
12. `mansion.jpg` — MAKIMA initialization
13. `mansion.jpg` — MAKIMA + Gappy-SVD
14. `mansion.jpg` — white-fill initialization
15. `mansion.jpg` — white-fill + Gappy-SVD
16. `palm.jpg` — original
17. `palm.jpg` — 40% missing
18. `palm.jpg` — white-fill initialization
19. `palm.jpg` — white-fill + Gappy-SVD
20. `palm.jpg` — MAKIMA initialization
21. `palm.jpg` — MAKIMA + Gappy-SVD

---

## Mathematical connection to finance

Matrix completion underpins several core problems in quantitative finance:

- **Covariance matrix estimation:** financial return matrices are often
  incomplete (assets with different trading histories, missing prices during
  halts). Low-rank SVD imputation regularizes the covariance estimator
  before portfolio optimization, analogously to how the rank-$m$ projection
  here fills missing entries.
- **Yield curve interpolation:** reconstructing interest rate surfaces from
  sparse market quotes at discrete maturities and tenors.
- **Factor models and collaborative filtering:** estimating unobserved asset
  returns using latent factor structure — the singular vectors $u_j$, $v_j$
  play the role of market factors.

The iterative SVD algorithm implemented here is a direct precursor to
**nuclear norm minimization** and **alternating least squares** methods
used in modern compressed sensing and recommender systems (Netflix Prize).
The initialization comparison (MAKIMA vs. zero vs. white-fill) parallels
the "warm-start" strategies used in convex solvers for these problems.
