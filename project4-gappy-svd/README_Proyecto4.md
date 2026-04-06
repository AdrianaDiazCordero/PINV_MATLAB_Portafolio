# Project 4 — Gappy Data Recovery via Iterative SVD

Matrix completion for functions and images with randomly missing entries,
using low-rank approximation and iterative singular value decomposition.

---

## Problem statement

Given a matrix $A \in \mathbb{R}^{n \times n}$ where a fraction $p$ of entries
are missing (set to NaN), recover the complete matrix by exploiting the
assumption that $A$ has **low numerical rank**.

Two data types are studied:
- **Synthetic function:** a complex 3D function $f(x,y,z)$ normalized to $[-1,1]$,
  sampled on a $100 \times 100$ grid for fixed $z$ slices.
- **Real grayscale images:** two photographs with up to 40% of pixels removed at random.

---

## Methods implemented

### Initialization strategies

| Strategy | Description |
|----------|-------------|
| Zero-fill | Replace NaN with 0 |
| White-fill | Replace NaN with 1 (max intensity) |
| MAKIMA interpolation | Replace NaN via modified Akima spline along columns |

### Core algorithm: iterative SVD (matrix completion)

```
Given: A_gappy (with NaN replaced by initialization value)
       m = number of modes, s = number of iterations

for k = 1 to s:
    [U, S, V] = svd(A_current, "econ")
    A_k = U[:,1:m] * S[1:m,1:m] * V[:,1:m]'    # rank-m approximation
    A_current[missing positions] = A_k[missing positions]  # fill only missing
end
```

The known entries are **never modified** — only the missing positions are
updated at each iteration using the low-rank projection.

### Mode selection
The optimal number of modes $m$ is determined by inspecting the singular
value spectrum: modes with $\sigma_j > 10^{-10}$ are retained.
For images, $m = \lfloor 0.1 \cdot m_{op} \rfloor$ (10% of numerical rank).

---

## Key results

### Synthetic function ($p = 40\%$ missing, $z = 0.5$ slice)

| Method | Max error | RMSE |
|--------|-----------|------|
| Zero-fill only | — | — |
| MAKIMA only | reported | reported |
| Zero-fill + SVD ($m=20$, $s=20$) | reported | reported |
| MAKIMA + SVD ($m=14$, $s=20$) | **best** | **best** |

> MAKIMA initialization consistently outperforms zero/white-fill as a
> starting point for the iterative SVD, requiring fewer iterations to converge.

### Image reconstruction ($p = 40\%$ missing)

Both images reconstructed with $m = \lfloor 0.1 \cdot m_{op} \rfloor$ modes
and $s = 50$ iterations starting from MAKIMA/white-fill initialization.

---

## Files

```
project4-gappy-svd/
├── CopiaProyecto4.m    # Main script
├── mansion.jpg           # Test image 1 (required)
├── palm.jpg              # Test image 2 (required)
└── README.md
```

> **Note:** add your own `mansion.jpg` and `palm.jpg` to the folder,
> or replace the filenames in the script with any grayscale-convertible image.

---

## How to run

```matlab
% Place image files in the same folder as the script, then:
run('CopiaProyecto4.m')
```

Generated figures (in order):
1. Original function (unnormalized)
2. Normalized function
3. Function with missing data
4. MAKIMA reconstruction
5. SVD reconstruction (no iteration)
6. Iterative SVD reconstruction (zero-fill init)
7. Iterative SVD reconstruction (MAKIMA init)
8–13. Same sequence for Image 1
14–19. Same sequence for Image 2

---

## Mathematical connection to finance

Matrix completion underpins several core problems in quantitative finance:

- **Covariance matrix estimation:** financial return matrices are often
  incomplete (assets with different trading histories, missing prices).
  Low-rank SVD imputation is used to regularize the covariance estimator
  before portfolio optimization.
- **Yield curve interpolation:** reconstructing interest rate surfaces from
  sparse market quotes.
- **Collaborative filtering / factor models:** recommender-style approaches
  to estimate unobserved asset returns using latent factor structure.

The iterative SVD algorithm here is a direct precursor to more advanced
methods such as **nuclear norm minimization** used in modern compressed sensing.
