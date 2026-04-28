# Anti-correlation-Gaussian

Code accompanying the paper [*Gibbs Sampling using Anti-correlation Gaussian
Data Augmentation, with Applications to L1-ball-type
Models*](https://arxiv.org/abs/2309.09371).

## Files

- **`SourceCode.R`** — `l1ball.linreg(X, y, ...)`: a blocked Gibbs sampler
  for the L1-ball-prior linear regression model that uses the anti-correlation
  Gaussian data-augmentation step described in the paper.
- **`SliceSampler.R`** — scalar and directional slice samplers
  (Neal, [*Slice Sampling*](https://arxiv.org/pdf/physics/0009028.pdf), 2003).
  Used by `l1ball.linreg` to update κ; can be sourced standalone.
- **`LinearReg.R`** — toy variable-selection demo (`n = 300`, `p = 500`)
  that produces estimation, trace, and ACF plots.
- **`TruncatedMVN.R`** — `truncMVN(...)`: extension of the anti-correlation
  trick to sampling a multivariate normal truncated to a box, plus a 2-D demo.
- **`tests/test_samplers.R`** — self-contained smoke + correctness tests
  (see [Testing](#testing)).
- **`100206_tfMRI_MOTOR_LR.nii.gz`** — Git-LFS pointer for the fMRI dataset
  used in the paper. Requires `git lfs install` and `git lfs pull` to fetch
  the actual file.

## Dependencies

Required:

- R ≥ 4.0
- [`truncnorm`](https://cran.r-project.org/package=truncnorm)

Optional (only if you run the corresponding feature/demo):

- `BoomSpikeSlab` &nbsp;— `init_method = "SSVS"` in `l1ball.linreg`
- `tmvtnorm`, `mvtnorm`, `dplyr`, `tidyr`, `tibble`, `ggplot2`, `gridExtra`
  &nbsp;— `TruncatedMVN.R` demo
- `ggplot2` &nbsp;— `LinearReg.R` demo

Install in R:

```r
install.packages(c("truncnorm", "ggplot2"))
# Optional extras:
install.packages(c("BoomSpikeSlab", "tmvtnorm", "mvtnorm",
                   "dplyr", "tidyr", "tibble", "gridExtra"))
```

## Quick start

```r
source("SliceSampler.R")
source("SourceCode.R")

set.seed(1)
n <- 100; p <- 20
X <- matrix(rnorm(n * p), n, p)
beta_true <- c(rep(2, 5), rep(0, p - 5))
y <- as.numeric(X %*% beta_true + rnorm(n))

fit <- l1ball.linreg(X, y, steps = 4000, burnin = 2000, thin = 1,
                     init_method = "random", verbose = FALSE)
round(fit$theta.est, 2)
```

To reproduce the paper's linear-regression figures, run

```bash
Rscript LinearReg.R
```

## Testing

A small test script is provided that

1. checks `truncMVN` reproduces the moments of a known truncated MVN, and
2. checks `l1ball.linreg` recovers the support and signs of the true
   coefficients on a small synthetic problem.

Run from the repo root:

```bash
Rscript tests/test_samplers.R
```

## Notes on the fMRI data

`100206_tfMRI_MOTOR_LR.nii.gz` is tracked with Git LFS. After cloning, run

```bash
git lfs install
git lfs pull
```

to download the binary. None of the R scripts in this repository load it
directly; it is included for reproducibility of the paper's fMRI experiments.

## License

[MIT](LICENSE).
