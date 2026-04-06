# Inverse Problems & Image Reconstruction — MATLAB Portfolio

**Adriana Díaz Cordero**  
MSc candidate · Industrial Mathematics, Universidade de Santiago de Compostela  
[![ORCID](https://img.shields.io/badge/ORCID-0009--0007--8744--7997-brightgreen)](https://orcid.org/0009-0007-8744-7997)

---

## About this repository

This repository collects four numerical projects developed for the course
*Inverse Problems and Image Reconstruction* (Máster en Matemática Industrial, USC, 2025–2026).

Each project explores a different class of ill-posed problem and implements
regularization strategies to obtain stable, physically meaningful solutions.
The mathematical backbone — compact operators, singular value decomposition,
and Tikhonov-type penalization — connects directly to quantitative finance
applications such as **implied volatility calibration** in the Black-Scholes framework.

---

## Projects

| # | Topic | Key methods | Folder |
|---|-------|-------------|--------|
| 1 | Gravitational inverse problem | Least squares, L-curve, cross-validation | [`project1-gravity/`](./project1-gravity/) |
| 2 | Heat equation — finite differences | Forward/backward FD, inverse crime analysis | [`project2-heat-fd/`](./project2-heat-fd/) |
| 3 | Heat equation — spectral regularization | Tikhonov, Landweber, TSVD, WSVD | [`project3-heat-regularization/`](./project3-heat-regularization/) |
| 4 | Gappy data & image reconstruction | Iterative SVD, matrix completion | [`project4-gappy-svd/`](./project4-gappy-svd/) |

---

## Mathematical background

All four projects deal with **ill-posed problems** in the sense of Hadamard:
the solution does not depend continuously on the data, and direct inversion
amplifies measurement noise catastrophically.

The core decomposition used throughout is the **Singular Value Decomposition**:

$$A = U \Sigma V^T, \qquad F = A^\dagger G = \sum_{j=1}^{r} \frac{1}{\sigma_j} \langle u_j, G \rangle \, v_j$$

Regularization strategies replace or modify this inversion to control
high-frequency noise amplification.

### Connection to quantitative finance

The heat equation

$$u_t = D \, u_{xx}$$

is mathematically equivalent to the **Black-Scholes PDE** after a standard change of variables
($x = \ln S$, $\tau = T - t$, $D = \sigma^2/2$). Consequently:

- **Calibrating implied volatility** from option market prices is an inverse problem
  of exactly the type studied in Projects 2 and 3.
- **Matrix completion** (Project 4) underlies modern approaches to **covariance matrix
  estimation** with missing or noisy financial time series.

---

## Requirements

- MATLAB R2021a or later (tested on R2024a)
- No external toolboxes required (Projects 1–3)
- Image Processing Toolbox: optional for Project 4 image examples

---

## Author

**Adriana Díaz Cordero** | adribarza16@gmail.com  
Sofía Kovalevskaya Prize — best MSc thesis in Mathematics, Cuba, 2025  
[LinkedIn](https://linkedin.com/in/adriana-diaz-cordero) · [ORCID](https://orcid.org/0009-0007-8744-7997)
