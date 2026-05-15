# ============================================================================
# Synthetic data generator for the TVP-SVD replication package
# ----------------------------------------------------------------------------
# Generates a single-equation TVP regression
#
#   y_t = x_t' beta_t + e_t,   e_t ~ N(0, sigma^2)
#
# with M_x = 4 regressors and T = 200 observations.  The coefficient path
# beta_t evolves in K + 1 = 3 piecewise-constant regimes, i.e. it features
# a small number of structural breaks rather than smooth random-walk
# drifts.  This setup is ideally suited for the WN-SVD sparse-mixture
# specification of Hauzenberger, Huber, Koop & Onorante (2021).
#
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>  REPLACE WITH YOUR DATA  <<<<<<<<<<<<<<<<<<<<<
# In a real application, replace simulate_data() with code that returns a
# list with elements y (T x 1 matrix) and X (T x M_x matrix).  For ease of
# downstream diagnostics simulate_data() also returns the true beta path
# and the noise standard deviation; in real-data use cases these can be
# set to NULL.
# ============================================================================

simulate_data <- function(T   = 200,
                          M_x = 4,
                          sigma = 0.5,
                          seed = 1234){
  set.seed(seed)

  # design matrix: M_x standardised regressors (first one acts as intercept)
  X <- cbind(1, matrix(rnorm(T*(M_x - 1)), T, M_x - 1))
  colnames(X) <- c("intercept", paste0("x", 1:(M_x - 1)))

  # true piecewise-constant coefficient paths -----------------------------
  # three regimes with break points at t = 70 and t = 140
  break1 <- 70
  break2 <- 140

  beta_true <- matrix(NA, T, M_x)
  # intercept
  beta_true[1:break1,           1] <-  1.0
  beta_true[(break1 + 1):break2, 1] <- -0.5
  beta_true[(break2 + 1):T,     1] <-  0.7
  # second coefficient
  beta_true[1:break1,           2] <-  0.8
  beta_true[(break1 + 1):break2, 2] <-  0.8   # stable across first break
  beta_true[(break2 + 1):T,     2] <- -1.2
  # third coefficient
  beta_true[1:break1,           3] <- -0.4
  beta_true[(break1 + 1):break2, 3] <-  1.5
  beta_true[(break2 + 1):T,     3] <-  0.2
  # fourth coefficient (essentially zero throughout)
  beta_true[, 4] <- 0.0

  # generate y
  y <- matrix(NA, T, 1)
  for (t in 1:T){
    y[t, 1] <- sum(X[t, ]*beta_true[t, ]) + rnorm(1, 0, sigma)
  }

  list(y = y, X = X, beta_true = beta_true, sigma = sigma)
}
