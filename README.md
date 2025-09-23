# RegCombinBLP

We implement best linear predictions in a context where the outcome of interest and some of the covariates are observed in two different datasets that cannot be matched based on **D'Haultfoeuille, Gaillac, Maurel (2024) <doi:10.48550/arXiv.2412.04816>** (DGM hereafter). The package allows for common regressors observed in both datasets, including auxiliary variables, which researchers do not seek to include in the regression but that appear in both datasets. It also provides asymptotically normal estimators of the bounds.

This README shows how to:

- Install the package from GitHub
- Simulate a small example dataset, used in the simulation section of **D'Haultfoeuille, Gaillac, Maurel (2024) <doi:10.48550/arXiv.2412.04816>**.
- Run `regCombin_BLP()` in three variants of the data-generating process (DGP):
  1) **`DGP_without`** — only the outside regressors (`Xnc`), not observed together with the outcome variable (`Y`).
  2) **`DGP_withXc`** — `Xnc` plus a commonly observed regressor `Xc`
  3) **`DGP_with_all`** — `Xnc`, `Xc`, and an observed auxiliary regressor `Wc`, which does not enters the regression.

---

## Installation

You’ll need a recent R (≥ 4.1 recommended) and the **devtools** package.

```r
# If needed
install.packages("devtools")

# Install core dependency (if your workflow uses it)
devtools::install_github("cgaillac/RegCombin")

# Install this package
devtools::install_github("cgaillac/RegCombinBLP")
```

Other packages used in the examples:

```r
install.packages(c("dplyr", "pracma", "sfsmisc", "Hmisc", "MASS", "snowfall"))
```

> If you plan on parallel runs, make sure your environment supports it, snowfall is installed, and consider setting `nbCores > 1` in `regCombin_BLP()`.

---

## Quick start

Every example below:

1. Simulates data
2. Assembles two data frames:
   - `Ldata`: outcome and “left” (Y-side) regressors you want the BLP for
   - `Rdata`: “right” (X-side) regressors available to combine
3. Calls `regCombin_BLP()` and prints:
   - `out$pt` — point estimates for the  DGM bounds on the BLP (if available for the chosen options)
   - `out$ci` — confidence intervals for the DGM bounds on the BLP

> Notes on key arguments we’ll reuse:
> - `out_var`: name of the outcome in `Ldata`
> - `nc_var`: vector of names of **outside** regressors, not commonly observed with the outcome, here `"Xnc"`
> - `c_var`: vector of **commonly observed** regressors , here `"Xc"` (optional)
> - `w_var`: vector of observed **auxiliary** regressors, commonly observed with the outcome but not entering the regression, here `"Wc"` (optional)
> - We set `DP = FALSE`, `full = FALSE`, `FW = FALSE`, `discrete = TRUE`, `bootstrap = FALSE` for a fast illustrative run. 

---

## Simulation and run helpers

These functions generate the DGP for the simulations, and store data in two separate dataframes `Ldata`, `Rdata`. There are three labels referring to the three different cases described in the simulation part of **D'Haultfoeuille, Gaillac, Maurel (2024) <doi:10.48550/arXiv.2412.04816>** and below.

```r
library(dplyr)
library(MASS)
library(sfsmisc)
library(pracma)
library(RegCombin)
library(RegCombinBLP)

simulate_dgp <- function(
  n = 1200,
  sig_eps = 4, sig_eta = 1,
  a1 = 1, a2 = 10, b1 = 1, b2 = 1, d1 = 1, d2 = 0,
  DGP = c("DGP_without", "DGP_withXc", "DGP_with_all"),
  seed = 3101989
) {
  set.seed(seed)
  DGP <- match.arg(DGP)

  # shocks
  eta_x <- rnorm(n, 0, sig_eta)
  eps_x <- rnorm(n, 0, sig_eps)
  eta_y <- rnorm(n, 0, sig_eta)
  eps_y <- rnorm(n, 0, sig_eps)

  # common components
  Xc_x <- rnorm(n, 0, 1)
  Xc_y <- rnorm(n, 0, 1)

  Wa_x <- runif(n, 0, 1)
  Wa_y <- runif(n, 0, 1)

  # observed noisy "outside" regressor and its Y-side counterpart
  Xo_x <-  Xc_x * a1 + Wa_x * a2 + (1 + Wa_x * d1) * eta_x
  Xo_y <-  Xc_y * a1 + Wa_y * a2 + (1 + Wa_y * d1) * eta_y

  # outcome
  Y <- Xo_y * b1 + Xc_y * b2 + Wa_y * d2 + eps_y

  # rename for role assignment
  out_var <- "Y"
  nc_var  <- "Xnc"      # non-classical regressor name
  c_var   <- NULL       # classical regressor(s)
  w_var   <- NULL       # auxiliary/instrument(s)

  if (DGP == "DGP_withXc") {
    c_var <- "Xc"
  } else if (DGP == "DGP_with_all") {
    c_var <- "Xc"
    w_var <- "Wa"
  }

  # Map raw vectors to role names per DGP
  Xc_y_keep <- if (is.null(c_var)) NULL else Xc_y
  Xc_x_keep <- if (is.null(c_var)) NULL else Xc_x
  Wa_y_keep <- if (is.null(w_var)) NULL else Wa_y
  Wa_x_keep <- if (is.null(w_var)) NULL else Wa_x

  # Build the two data frames with aligned names
  L_cols <- list(Y, Xc_y_keep, Wa_y_keep)
  Ldata <- as.data.frame(do.call(cbind, L_cols[!sapply(L_cols, is.null)]))
  colnames(Ldata) <- c(out_var, c(if (!is.null(c_var)) c_var else NULL),
                                  if (!is.null(w_var)) w_var else NULL)

  R_cols <- list(Xo_x, Xc_x_keep, Wa_x_keep)
  Rdata <- as.data.frame(do.call(cbind, R_cols[!sapply(R_cols, is.null)]))
  colnames(Rdata) <- c(nc_var, c(if (!is.null(c_var)) c_var else NULL),
                                   if (!is.null(w_var)) w_var else NULL)

  list(
    Ldata = Ldata,
    Rdata = Rdata,
    out_var = out_var,
    nc_var = nc_var,
    c_var = c_var,
    w_var = w_var
  )
}

run_blp <- function(Ldata, Rdata, out_var, nc_var, c_var = NULL, w_var = NULL,
                    nbCores = 1, FW = FALSE, discrete = TRUE, pt_est = FALSE,
                    bootstrap = FALSE, dK = 2.5, ASN = TRUE, K = 0,
                    dataset = 0, O_only = TRUE) {

  regCombin_BLP(
    Ldata = Ldata, Rdata = Rdata,
    out_var = out_var, nc_var = nc_var,
    c_var = c_var, w_var = w_var,
    w_x = NULL, w_y = NULL, nbCores = nbCores,
    DP = FALSE, full = FALSE, FW = FW, discrete = discrete,
    pt_est = pt_est, bootstrap = bootstrap, dK = dK,
    ASN = ASN, unchanged = FALSE, factor = FALSE, K_sel = K,
    dataset = dataset, O_only = O_only
  )
}
```

---

## Case 1 — `DGP_without` (only `Xnc`, without (`Xc`,`Wa`))

```r
sim1 <- simulate_dgp(DGP = "DGP_without")
out1 <- run_blp(sim1$Ldata, sim1$Rdata,
                out_var = sim1$out_var,
                nc_var  = sim1$nc_var,
                c_var   = sim1$c_var,   # NULL
                w_var   = sim1$w_var)   # NULL

# Point estimates for the DGM bounds
out1$pt

# Confidence intervals for the DGM bounds
out1$ci
```

**Interpretation.** With only the outside regressor available on the right-hand side, the procedure reports bounds for the BLP.

---

## Case 2 — `DGP_withXc` (`Xnc` + commonly observed `Xc`)

```r
sim2 <- simulate_dgp(DGP = "DGP_withXc")
out2 <- run_blp(sim2$Ldata, sim2$Rdata,
                out_var = sim2$out_var,
                nc_var  = sim2$nc_var,
                c_var   = sim2$c_var,   # "Xc"
                w_var   = sim2$w_var)   # NULL

out2$pt
out2$ci
```

**Interpretation.** Adding a commonly observed regressor `Xc` might help tightening bounds, and in some settings enables point identification, see section "How do common regressors affect identification?" in DGM.

---

## Case 3 — `DGP_with_all` (`Xnc` + `Xc` + auxiliary `Wa`)

```r
sim3 <- simulate_dgp(DGP = "DGP_with_all")
out3 <- run_blp(sim3$Ldata, sim3$Rdata,
                out_var = sim3$out_var,
                nc_var  = sim3$nc_var,
                c_var   = sim3$c_var,   # "Xc"
                w_var   = sim3$w_var)   # "Wa"

out3$pt
out3$ci
```

**Interpretation.** Supplying an observed auxiliary variable `Wa` typically tightens bounds further. On this example  with n=1200, it excludes 0 from the confidence interval on the coefficient of `Xnc'.

---

## Tips & troubleshooting

- **Speed / parallelism.** Increase `nbCores` if your machine has multiple cores and your OS/BLAS supports parallelism.
- **Point estimation.** If you require only the point estimates, consider `pt_est = TRUE`.
- **Bootstrap inference.** Set `bootstrap = TRUE`; this will take longer than using asymptotic normality results.

---

## Citation

If you use **RegCombinBLP** in academic work, please cite D'Haultfoeuille, Gaillac, Maurel (2024) <doi:10.48550/arXiv.2412.04816>, available at https://arxiv.org/abs/2412.04816.

---

## MIT License

---

