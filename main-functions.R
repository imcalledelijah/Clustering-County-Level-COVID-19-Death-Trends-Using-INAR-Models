#' @title Binomial Thinning Operator
#' @description Implements the binomial thinning operation for INAR processes
#' @param x Non-negative integer value to be thinned
#' @param alpha Thinning probability (between 0 and 1)
#' @return Thinned count value
#' @examples
#' bin_thinning(10, 0.5)  # Expected value around 5
#' bin_thinning(0, 0.5)   # Always returns 0
#' @export
bin_thinning <- function(x, alpha) {
  # Handle x=0 case explicitly
  if (x == 0) return(0)
  
  # Input validation; checking if the values of x and alpha suffice
  if (!is.numeric(x) || x < 0 || x %% 1 != 0) {
    stop("x must be a non-negative integer")
  }
  if (alpha < 0 || alpha > 1) {
    stop("alpha must be between 0 and 1")
  }
  
  # Perform binomial thinning
  sum(rbinom(x, size = 1, prob = alpha))
}

#' @title Simulate INAR(p*) Process
#' @description Generates a single INAR(p*) time series with specified parameters
#' @param t Length of time series
#' @param alpha Autoregressive parameter (between 0 and 1)
#' @param innov Vector of innovation terms (length >= n)
#' @param p Lag order for INAR process
#' @return Numeric vector containing the simulated time series
#' @examples
#' innov <- rpois(50, lambda = 5); 
#' gives me a vector of 50 values generated from Pois(5)
#' 
#' sim_inar_p_ast(50, 0.5, innov, p = 1); 
#' simulate a time series of length 50 with alpha = 0.5 
#' and a vector of innovations generated from Pois(5) with 1 lag
#' @export
sim_inar_p_ast <- function(t, alpha, innov, p) {
  # Initialize empty vector to store time series values
  x <- numeric(t)
  
  # Set first p values using initial innovations (burn-in period)
  x[1:p] <- innov[1:p]
  
  # Generate subsequent values using INAR(p*) structure
  for (i in (p + 1):t) {
    # Apply binomial thinning to the p-th lagged value
    thinned <- bin_thinning(x[i - p], alpha)
    
    # Current value = thinned p-lag value + current innovation
    x[i] <- thinned + innov[i]
  }
  
  # Return generated time series
  x
}
#' series = sim_inar_p_ast(50, 0.5, innov, p = 3)
#' acf(series)

#' @title Simulate Mixture of INAR Processes
#' @description Generates multiple time series from a mixture of INAR(p*) processes
#' @param n Number of time series to generate
#' @param t Length of each time series
#' @param mix_prop Vector of mixing proportions (should sum to 1)
#' @param params List of parameter lists for each component, containing:
#' \itemize{
#'   \item p - Lag order
#'   \item alpha - Autoregressive parameter
#'   \item lambda - Innovation mean
#'   \item family - Innovation distribution ("pois" or "nbinom")
#'   \item phi - Dispersion parameter (for negative binomial)
#' }
#' @return List with two elements:
#' \itemize{
#'   \item data - Matrix of simulated time series (n x t)
#'   \item z - True cluster assignments
#' }
#' @examples
#' params <- list(
#'   list(p = 1, alpha = 0.5, lambda = 5, family = "pois", phi = NULL),
#'   list(p = 2, alpha = 0.3, lambda = 8, family = "nbinom", phi = 2)
#' )
#' Above we have a list containing two sets of parameters where we do a mixture of 
#' Poisson and Neg Binomial Processes. The first process is a Poisson Process
#' of lag order 1, alpha = 0.5, lambda = 5. The second process is a Negative 
#' Binomial process of lag order 2, alpha = 0.3, lambda = 8, and phi = 2
#' 
#' sim_mixture_inar_p_ast(10, 50, c(0.4, 0.6), params)
#' Here we simulate 10 time series where each is length 50 and the probability
#' of the first one is 40% and the second one is 60% and finally we have an argument
#' of all the parameters
#' 
#' @export

sim_mixture_inar_p_ast <- function(n, t, mix_prop, params) {
  # Determine number of clusters/components from mixing proportions
  G <- length(mix_prop)
  
  # Assign each time series to a cluster (z vector contains true memberships)
  z <- sample(1:G, size = n, prob = mix_prop, replace = TRUE)
  
  # Initialize empty matrix to store all time series (rows = series, columns = time points)
  data <- matrix(0, nrow = n, ncol = t)
  
  # Generate each time series
  for (i in 1:n) {
    # Get cluster assignment for this series
    g <- z[i]
    
    # Extract parameters for this cluster
    p <- params[[g]]$p       # Lag order (how many previous time points affect current)
    alpha <- params[[g]]$alpha  # Autoregressive parameter (survival probability)
    lambda <- params[[g]]$lambda  # Innovation mean (average new events per period)
    family <- params[[g]]$family  # Distribution type: "pois" or "nbinom"
    phi <- params[[g]]$phi    # Dispersion parameter (for negative binomial only)
    
    # Generate innovation sequence (random shocks/noise)
    innov <- if (family == "pois") {
      # Poisson innovations: lambda = mean events per time period
      rpois(t, lambda)
    } else if (family == "nbinom") {
      # Negative binomial innovations: handles overdispersed counts
      # mu = mean events, size = dispersion (smaller = more dispersion)
      rnbinom(t, mu = lambda, size = phi)
    } else {
      stop("Invalid family")
    }
    
    # Generate INAR(p*) time series using binomial thinning
    # Note: sim_inar_p_ast should be the INAR simulation function from previous code
    data[i, ] <- sim_inar_p_ast(t, alpha, innov, p)
  }
  
  # Return simulated data and true cluster assignments
  list(data = data, z = z)
}

#' \itemize{
#'   \item p - Lag order
#'   \item alpha - Autoregressive parameter
#'   \item lambda - Innovation mean
#'   \item family - Innovation distribution ("pois" or "nbinom")
#'   \item phi - Dispersion parameter (for negative binomial)
#' }
#' @return List with two elements:
#' \itemize{
#'   \item data - Matrix of simulated time series (n x t)
#'   \item z - True cluster assignments
#' }

params = list(
   list(p = 1, alpha = 0.5, lambda = 5, family = "pois", phi = NULL),
   list(p = 2, alpha = 0.3, lambda = 8, family = "nbinom", phi = 2))

params = list(
  list(p = 1, alpha = 0.5, lambda = 5, family = "pois", phi = NULL),
  list(p = 2, alpha = 0.3, lambda = 8, family = "pois", phi = NULL))

#' sim_mixture_inar_p_ast(n = 50, t = 200, c(0.4, 0.6), params)


#' @title Compute INAR Log-Likelihood
#' @description Calculates the log-likelihood for a single time series under an INAR(p*) model
#' @param x Time series vector
#' @param p Lag order
#' @param alpha Autoregressive parameter
#' @param lambda Innovation mean parameter
#' @param family Innovation distribution ("pois" or "nbinom")
#' @param phi Dispersion parameter (for negative binomial)
#' @return Log-likelihood value
#' @examples
#' 
#' x <- rpois(50, 5)
#' x here is a vector of 50 values generated from a Pois(5) distribution
#' 
#' compute_loglik(x, p = 1, alpha = 0.5, lambda = 5, family = "pois", phi = NULL)
#' this computes the log likelihood of a vector and seeing if the series fall
#' under a Pois(5) distribution
#' @export
compute_loglik <- function(x, p, alpha, lambda, family, phi) {
  # Get total length of the time series
  T_ <- length(x)
  
  # Initialize log-likelihood accumulator
  loglik <- 0
  
  # Handle INITIAL OBSERVATIONS (first p points):
  # These only depend on the innovation distribution
  if (family == "pois") {
    # For Poisson: sum log-probabilities of first p observations
    loglik <- sum(dpois(x[1:p], lambda, log = TRUE))
  } else if (family == "nbinom") {
    # For Negative Binomial: sum log-probabilities of first p observations
    loglik <- sum(dnbinom(x[1:p], mu = lambda, size = phi, log = TRUE))
  }
  
  # Handle SUBSEQUENT OBSERVATIONS (from t = p + 1 to end):
  # These depend on binomial thinning + innovations
  for (t in (p + 1):T_) {
    # Get previous observation at lag p
    m <- x[t - p]  
    # Get current observation
    k <- x[t]      
    
    # Calculate possible number of survivors (s) from binomial thinning:
    # Can't have more survivors than either current count (k) or previous count (m)
    s_max <- min(k, m)
    # Create sequence of possible survivors (0 to s_max)
    s <- 0:s_max
    
    # Calculate probability for all possible s values:
    if (family == "pois") {
      # Poisson case: sum over binomial thinning probability * innovation probability
      prob <- sum(dbinom(s, m, alpha) *       # Binomial thinning probability
                    dpois(k - s, lambda))       # Innovation probability
    } else if (family == "nbinom") {
      # Negative Binomial case: same structure with different innovation distribution
      prob <- sum(dbinom(s, m, alpha) * 
                    dnbinom(k - s, mu = lambda, size = phi))
    }
    
    # Add log-probability to total log-likelihood:
    # Use 1e-10 to avoid log(0) errors (numerical stability)
    loglik <- loglik + log(prob + 1e-10)
  }
  
  # Return final computed log-likelihood
  loglik
}

set.seed(539)
x <- rpois(50, 5)
# x here is a vector of 50 values generated from a Pois(5) distribution

compute_loglik(x, p = 1, alpha = 0.5, lambda = 5, family = "pois", phi = NULL)
# this computes the log likelihood of our vector above

loglike = numeric(20)
for(i in 1:20){
  loglike[i] = compute_loglik(x, p = 1, alpha = 0.5, lambda = i, family = "pois", phi = NULL)
}
loglike

# this log likelihood takes care of patients coming in AND out of the hospital

#' @title Initialize Parameters using k-means
#' @description Provides initial parameter estimates for EM algorithm using k-means clustering
#' @param data Matrix of time series (n x t)
#' @param G Number of clusters
#' @param p_vec Vector of lag orders for each cluster
#' @param family_list List of innovation families for each cluster
#' @return List containing:
#' \itemize{
#'   \item params - List of parameter estimates for each cluster
#'   \item pi - Initial mixing proportions
#' }
#' @examples
#' data <- matrix(rpois(100*50, 5), nrow = 100)
#' initialize_params(data, G = 2, p_vec = c(1,1), family_list = list("pois", "pois"))
#' @export
initialize_params <- function(data, G, p_vec, family_list) {
  # Get dimensions: n = number of time series, t = length of each series
  n <- nrow(data)
  t <- ncol(data)
  
  # Create feature matrix (means and autocorrelations at different lags)
  # Features help distinguish different time series patterns
  features <- matrix(0, nrow = n, ncol = G)
  
  # Calculate autocorrelation features for each cluster's specified lag
  for (g in 1:G) {
    p <- p_vec[g]  # Get lag order for this cluster
    features[, g] <- apply(data, 1, function(x) {
      if (p >= t) return(0)  # Handle cases where lag exceeds series length
      cor(x[(p + 1):t], x[1:(t - p)])  # Calculate autocorrelation at lag p
    })
  }
  
  # Add mean of each series as first feature column
  features <- cbind(apply(data, 1, mean), features)
  
  # Cluster time series using k-means based on these features
  km <- kmeans(features, centers = G)
  
  # Prepare to store parameters for each cluster
  params <- list()
  
  # Estimate parameters for each cluster
  for (g in 1:G) {
    idx <- km$cluster == g  # Identify series in this cluster
    x_g <- data[idx, ]      # Extract cluster members
    p_g <- p_vec[g]        # Get lag order for this cluster
    
    # Estimate autocorrelation parameter (alpha)
    acf_g <- apply(x_g, 1, function(x) {
      if (p_g >= t) return(0)
      cor(x[(p_g + 1):t], x[1:(t - p_g)])  # Calculate autocorrelations
    })
    alpha_init <- ifelse(p_g >= t, 0.5, mean(acf_g, na.rm = TRUE))  # Average autocorrelation
    alpha_init <- pmax(pmin(alpha_init, 0.99), 0.01)  # Keep between 0.01-0.99
    
    # Estimate innovation mean (lambda)
    mean_x <- mean(x_g)  # Average count in cluster
    lambda_init <- (1 - alpha_init) * mean_x  # Derived from INAR properties
    
    # Estimate dispersion parameter (phi) for negative binomial
    if (family_list[[g]] == "nbinom") {
      var_x <- var(x_g)  # Calculate variance
      phi_init <- mean_x^2 / (var_x - mean_x)  # Dispersion formula
      phi_init <- pmax(phi_init, 1e-3)  # Ensure positive value
    } else {
      phi_init <- NULL  # No dispersion for Poisson
    }
    
    # Store parameter estimates for this cluster
    params[[g]] <- list(
      p = p_g,         # Lag order
      alpha = alpha_init,  # Autoregressive parameter
      lambda = lambda_init,  # Innovation mean
      phi = phi_init,   # Dispersion (if negative binomial)
      family = family_list[[g]]  # Distribution family
    )
  }
  
  # Calculate mixing proportions (cluster sizes)
  pi <- km$size / n
  
  # Return both parameters and mixing proportions
  list(params = params, pi = pi)
}

data <- matrix(rpois(100*50, 5), nrow = 100)
# gets us a matrix of 100 time series that have length 50

initialize_params(data, G = 2, p_vec = c(1,1), family_list = list("pois", "pois"))
# gives us what our initial parameters are of the mixture model we want to use
# wants 2 clusters where each has lag 1 and they are both PINAR Models


#' ask about doing this function for negative binomial
#' data <- matrix(rnbinom(100*50, 5), nrow = 100)
#' gets us a matrix of 100 time series that have length 50

initialize_params(data, G = 2, p_vec = c(1,1), family_list = list("pois", "pois"))

#' @title EM Algorithm for INAR Mixture Models with Aitken Acceleration
#' @description Implements the Expectation-Maximization algorithm for finite mixtures of 
#' INAR(p*) processes using Aitken's acceleration for convergence detection. Follows the 
#' methodology from Roick et al. (2021) with numerical stability enhancements.
#' 
#' @param data A matrix of discrete-valued time series (n x t) where rows represent 
#' individuals and columns represent time points
#' @param G Integer specifying the number of mixture components/clusters
#' @param p_vec Vector of length G specifying the lag order (p*) for each component
#' @param family_list List of length G specifying innovation distributions:
#' \itemize{
#'   \item "pois" - Poisson distribution
#'   \item "nbinom" - Negative binomial distribution
#' }
#' @param max_iter Maximum number of EM iterations (default: 100)
#' @param tol Convergence tolerance for Aitken acceleration criterion (default: 1e-4)
#' 
#' @return A list containing:
#' \itemize{
#'   \item params - List of estimated parameters for each component (alpha, lambda, phi, p, family)
#'   \item pi - Estimated mixing proportions
#'   \item resp - Final responsibility matrix (n x G)
#'   \item log_lik - Final log-likelihood value
#'   \item converged - Logical indicating convergence status
#' }
#' 
#' @examples
#' \dontrun{
#' # Example using simulated Poisson INAR(1) data
#' set.seed(123)
#' params <- list(
#'   list(p = 1, alpha = 0.5, lambda = 5, family = "pois"),
#'   list(p = 1, alpha = 0.3, lambda = 8, family = "pois")
#' )
#' sim_data <- sim_mixture_inar_p_ast(100, 50, c(0.4, 0.6), params)
#' em_result <- em_inar(sim_data$data, G = 2, p_vec = c(1,1), 
#'                family_list = list("pois", "pois"))
#' }
#' 
#' @references 
#' Roick T, Karlis D, McNicholas PD (2021). Clustering discrete-valued time series. 
#' Advances in Data Analysis and Classification 15:209-229.\cr
#' McNicholas PD, Murphy TB, McDaid AF, Frost D (2010). Serial and parallel implementations 
#' of model-based clustering via parsimonious Gaussian mixture models. 
#' Computational Statistics & Data Analysis 54(3):711-723.
#' 
#' @importFrom matrixStats rowLogSumExps
#' @export
em_inar <- function(data, G, p_vec, family_list, max_iter = 2000, tol = 1e-3) {
  n <- nrow(data)
  t <- ncol(data)
  
  # Enhanced initialization with fallback values
  init <- tryCatch(
    initialize_params(data, G, p_vec, family_list),
    error = function(e) {
      warning("Initialization failed, using fallback: ", e$message)
      params <- lapply(1:G, function(g) {
        list(
          p = p_vec[g],
          alpha = runif(1, 0.1, 0.9),  # Random valid alpha
          lambda = pmax(mean(as.numeric(data)), 1e-3),  # Convert matrix to vector
          phi = if(family_list[[g]] == "nbinom") 1 else NULL,
          family = family_list[[g]]
        )
      })
      list(params = params, pi = rep(1/G, G))
    }
  )
  
  params <- init$params
  pi <- pmax(init$pi, 1e-3)
  pi <- pi / sum(pi)
  log_lik <- matrix(-Inf, nrow = n, ncol = G)
  
  # Calculate initial log-likelihoods with stability checks
  for (g in 1:G) {
    p_g <- params[[g]]$p
    alpha_g <- pmax(pmin(params[[g]]$alpha, 0.99), 0.01)
    lambda_g <- pmax(params[[g]]$lambda, 1e-3)
    phi_g <- if (family_list[[g]] == "nbinom") pmax(params[[g]]$phi, 1e-3) else NULL
    
    for (i in 1:n) {
      log_lik[i, g] <- tryCatch(
        compute_loglik(data[i, ], p_g, alpha_g, lambda_g, family_list[[g]], phi_g),
        error = function(e) -1e10
      )
    }
  }
  
  # EM Algorithm Core
  converged <- FALSE
  prev_loglik <- -Inf
  prev_prev_loglik <- -Inf
  iter <- 1
  
  while (!converged && iter <= max_iter) {
    # E-step: Calculate responsibilities
    log_resp <- sweep(log_lik, 2, log(pi), "+")
    log_resp_norm <- matrixStats::rowLogSumExps(log_resp)
    log_resp <- log_resp - log_resp_norm
    resp <- exp(log_resp)
    resp <- pmax(resp, .Machine$double.eps)
    resp <- resp / rowSums(resp)
    
    # M-step: Update parameters
    for (g in 1:G) {
      family_g <- family_list[[g]]
      p_g <- params[[g]]$p
      
      # Set optimization bounds
      lower <- if (family_g == "pois") c(1e-3, 1e-3) else c(1e-3, 1e-3, 1e-3)
      upper <- if (family_g == "pois") c(0.999, Inf) else c(0.999, Inf, Inf)
      init_par <- if (family_g == "pois") {
        c(params[[g]]$alpha, params[[g]]$lambda)
      } else {
        c(params[[g]]$alpha, params[[g]]$lambda, params[[g]]$phi)
      }
      
      # Optimization
      opt <- optim(
        init_par,
        fn = function(par) {
          alpha <- par[1]
          lambda <- par[2]
          phi <- if (family_g == "nbinom") par[3] else NULL
          -sum(resp[, g] * sapply(1:n, function(i) {
            compute_loglik(data[i, ], p_g, alpha, lambda, family_g, phi)
          }))
        },
        method = "L-BFGS-B",
        lower = lower,
        upper = upper,
        control = list(maxit = max_iter)
      )
      
      # Update parameters with validation
      params[[g]]$alpha <- pmax(pmin(opt$par[1], 0.99), 0.01)
      params[[g]]$lambda <- pmax(opt$par[2], 1e-3)
      if (family_g == "nbinom") params[[g]]$phi <- pmax(opt$par[3], 1e-3)
    }
    
    # Update mixing proportions
    pi <- colMeans(resp, na.rm = TRUE) %>% pmax(1e-3) %>% {. / sum(.)}
    
    # Convergence check with Aitken acceleration
    current_loglik <- sum(matrixStats::rowLogSumExps(sweep(log_lik, 2, log(pi), "+")))
    
    if (iter >= 2) {
      delta_prev <- prev_loglik - prev_prev_loglik
      delta_current <- current_loglik - prev_loglik
      
      if (delta_prev != 0) {
        a <- delta_current / delta_prev
        if (a < 1) {
          l_inf <- prev_loglik + delta_current / (1 - a)
          if (l_inf - current_loglik < tol) converged <- TRUE
        }
      }
    }
    
    if (abs(current_loglik - prev_loglik) < tol) converged <- TRUE
    
    prev_prev_loglik <- prev_loglik
    prev_loglik <- current_loglik
    iter <- iter + 1
  }
  
  list(params = params, pi = pi, resp = round(resp, 3), log_lik = current_loglik, converged = converged)
}


#' @title Calculate BIC for INAR Mixture Model
#' @description Implements model selection using Bayesian Information Criterion
#' @param em_result Output from em_inar()
#' @param n_obs Number of observations (total time points across all series)
#' @return BIC value
#' @export
calculate_bic <- function(em_result, n_obs) {
  params <- em_result$params  # Extract params from result object
  G <- length(em_result$pi)
  
  # Count parameters: (G-1) mixing proportions + 2*G (alpha, lambda) + G*phi (if nbinom)
  n_params <- (G - 1) + 2 * G + 
    sum(sapply(params, function(p) !is.null(p$phi)))
  
  bic <- -2 * em_result$log_lik + n_params * log(n_obs)
  bic
}


#' @title Automated Model Selection for INAR Mixture Models (Serial Version)
#' @description Performs exhaustive model selection across different numbers of components (G)
#'              and lag combinations using BIC. Tests all unique combinations of p-values for
#'              multi-component models.
#' 
#' @param data Matrix of time series (n x t) with rows as individuals and columns as time points
#' @param max_g Maximum number of clusters/components to consider
#' @param p_candidates Vector of candidate lag orders (p*) to combine
#' @param family Innovation distribution family ("pois" or "nbinom")
#' @param max_iter Maximum EM iterations per model (default: 100)
#' 
#' @return List containing:
#' \itemize{
#'   \item best_model - Model with lowest BIC
#'   \item all_models - List of all tested models with BIC values
#' }
#' 
#' @examples
#' \dontrun{
#' # Test up to 3 components with lags 5 and 10
#' data <- matrix(rpois(50*200, 5), nrow = 50)
#' result <- auto_model_selection_serial(data, max_g = 3, 
#'                                      p_candidates = c(5, 10),
#'                                      family = "pois")
#' }
#' 
#' @references 
#' Roick et al. (2021). Clustering discrete-valued time series. ADAC 15:209-229.
#' @export
auto_model_selection_serial <- function(
    data, max_g, p_candidates, family, max_iter = 100) 
{
  t <- ncol(data)
  feasible_p <- p_candidates[p_candidates < t]
  all_models <- list()
  
  for(g in 1:max_g) {
    # Generate p combinations
    if(g == 1) {
      p_combinations <- as.list(feasible_p)
    } else {
      # Generate all combinations (including permutations)
      p_combinations <- expand.grid(replicate(g, feasible_p, simplify = FALSE))
      # Convert to list of vectors
      p_combinations <- split(as.matrix(p_combinations), seq(nrow(p_combinations)))
    }
    
    for(p_vec in p_combinations) {
      tryCatch({
        family_list <- rep(list(family), g)
        p_vec <- unlist(p_vec)  # Ensure numeric vector
        
        message(sprintf("\nTesting G=%d with p=%s", g, paste(p_vec, collapse = ",")))
        
        fit <- em_inar(
          data = data,
          G = g,
          p_vec = p_vec,
          family_list = family_list,
          max_iter = max_iter
        )
        
        if(fit$converged) {
          bic <- calculate_bic(fit, nrow(data) * ncol(data))
          all_models[[length(all_models) + 1]] <- list(
            G = g,
            p = p_vec,
            bic = bic,
            fit = fit
          )
          message(sprintf("Converged! BIC=%.2f", bic))
        }
      }, error = function(e) {
        message(sprintf("Error: %s", e$message))
      })
    }
  }
  
  if(length(all_models) == 0) stop("All model configurations failed")
  
  bic_values <- sapply(all_models, function(x) x$bic)
  best_idx <- which.min(bic_values)
  
  list(
    best_model = all_models[[best_idx]],
    all_models = all_models
  )
}

auto_model_selection_parallel <- function(data, max_g, p_candidates, family, max_iter = 100) {
  library(doParallel)
  library(foreach)
  
  # Optional: Log output to file
  sink("parallel_model_log.txt", split = TRUE)
  
  t <- ncol(data)
  feasible_p <- as.numeric(p_candidates[p_candidates < t])
  all_models <- list()
  
  # Set up parallel backend
  n_cores <- parallel::detectCores() - 1
  cl <- makeCluster(n_cores)
  registerDoParallel(cl)
  
  # Export data and anything in global environment
  clusterExport(cl, varlist = c("data"))
  
  for (g in 1:max_g) {
    # Generate p combinations
    if (g == 1) {
      p_combinations <- as.list(feasible_p)
    } else {
      p_combinations <- expand.grid(replicate(g, feasible_p, simplify = FALSE), stringsAsFactors = FALSE)
      p_combinations <- split(as.matrix(p_combinations), seq(nrow(p_combinations)))
    }
    
    models_g <- foreach(p_vec = p_combinations,
                        .combine = list,
                        .multicombine = TRUE,
                        .packages = c("dplyr")) %dopar% {
                          
                          # ✅ Load your custom function definitions inside each worker
                          source("/Users/elijahamirianfar/My Drive (elijahl2l3l4@gmail.com)/04. CSU FULLERTON 2023-2025/RESEARCH/CODE/main-functions.R")
                          # Optional: source additional helper functions if needed
                          # source("/Users/elijahamirianfar/.../plot-functions-elijah.R")
                          
                          result <- NULL
                          tryCatch({
                            p_vec <- as.numeric(unlist(p_vec))
                            cat(sprintf("TRY: G = %d, p = %s\n", g, paste(p_vec, collapse = ",")))
                            
                            if (any(is.na(p_vec))) stop("p_vec contains NA")
                            
                            family_list <- rep(list(family), g)
                            
                            fit <- tryCatch({
                              em_inar(
                                data = data,
                                G = g,
                                p_vec = p_vec,
                                family_list = family_list,
                                max_iter = max_iter
                              )
                            }, error = function(e_fit) {
                              cat(sprintf("em_inar() failed for G=%d, p=%s: %s\n", g, paste(p_vec, collapse = ","), e_fit$message))
                              return(NULL)
                            })
                            
                            if (!is.null(fit)) {
                              cat("FIT STRUCTURE:\n")
                              print(str(fit))
                              
                              if (!is.null(fit$converged) && fit$converged) {
                                bic <- tryCatch({
                                  calculate_bic(fit, nrow(data) * ncol(data))
                                }, error = function(e_bic) {
                                  cat(sprintf("calculate_bic() failed for G=%d, p=%s: %s\n", g, paste(p_vec, collapse = ","), e_bic$message))
                                  return(NA_real_)
                                })
                                
                                if (!is.na(bic)) {
                                  result <- list(
                                    G = g,
                                    p = p_vec,
                                    bic = bic,
                                    fit = fit
                                  )
                                } else {
                                  cat(sprintf("BIC is NA for G = %d, p = %s\n", g, paste(p_vec, collapse = ",")))
                                }
                              } else {
                                cat(sprintf("Fit did not converge for G = %d, p = %s\n", g, paste(p_vec, collapse = ",")))
                              }
                            }
                          }, error = function(e) {
                            cat(sprintf("Outer error at G=%d, p=%s: %s\n", g, paste(p_vec, collapse = ","), e$message))
                          })
                          
                          result
                        }
    
    models_g <- Filter(Negate(is.null), models_g)
    all_models <- c(all_models, models_g)
  }
  
  # Stop parallel backend
  stopCluster(cl)
  sink()
  
  if (length(all_models) == 0) stop("All model configurations failed")
  
  bic_values <- sapply(all_models, function(x) {
    if (!is.null(x$bic)) x$bic else NA_real_
  })
  
  valid_idx <- which(!is.na(bic_values))
  if (length(valid_idx) == 0) stop("No valid BICs found")
  
  best_idx <- valid_idx[which.min(bic_values[valid_idx])]
  
  list(
    best_model = all_models[[best_idx]],
    all_models = all_models
  )
}


#' @title Assess Clustering Performance
#' @description Calculates the Adjusted Rand Index (ARI) and cross-tabulation between true and predicted cluster memberships.
#' @param z True cluster assignments (vector of integers)
#' @param resp Responsibility matrix from EM algorithm (n x G matrix)
#' @return List containing:
#' \itemize{
#'   \item ARI - Adjusted Rand Index value
#'   \item ContingencyTable - Cross-tabulation of true vs predicted clusters
#' }
#' @examples
#' \dontrun{
#' # Example with perfect agreement
#' z <- c(1,1,2,2)
#' resp <- matrix(c(0.9,0.1, 0.8,0.2, 0.1,0.9, 0.2,0.8), nrow=4, byrow=TRUE)
#' assess_performance(z, resp)  # Should return ARI = 1
#' 
#' # Example with random assignment
#' resp_random <- matrix(0.5, nrow=4, ncol=2)
#' assess_performance(z, resp_random)  # ARI near 0
#' }
#' @export
assess_performance <- function(z, resp) {
  # Get MAP classification from responsibility matrix
  map <- apply(resp, 1, which.max)
  
  # Create contingency table
  ct <- table(True = z, Predicted = map)
  
  # Calculate components for ARI
  n <- length(z)
  ct_matrix <- as.matrix(ct)
  
  sum_nij <- sum(choose(ct_matrix, 2))  # Sum of combinations for each cell
  a <- rowSums(ct_matrix)               # Row sums (true cluster sizes)
  sum_a <- sum(choose(a, 2))            # Sum of combinations for rows
  b <- colSums(ct_matrix)               # Column sums (predicted cluster sizes)
  sum_b <- sum(choose(b, 2))            # Sum of combinations for columns
  total_pairs <- choose(n, 2)           # Total possible pairs
  
  # Calculate expected agreement
  expected <- (sum_a * sum_b) / total_pairs
  
  # Calculate numerator and denominator for ARI
  numerator <- sum_nij - expected
  denominator <- 0.5 * (sum_a + sum_b) - expected
  
  # Handle edge case where denominator is 0
  if (denominator == 0) {
    ari <- 0
  } else {
    ari <- numerator / denominator
  }
  
  # Return results
  list(ARI = ari, ContingencyTable = ct)
}





