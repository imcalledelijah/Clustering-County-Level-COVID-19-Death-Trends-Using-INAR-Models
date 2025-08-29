#' @title INAR(p*) Process Simulation and Analysis
#' @description Generates and analyzes INAR(p*) process with flexible innovation distribution
#' @param t Total time points (excluding burn-in)
#' @param burn_in Number of initial observations to discard
#' @param lambda Mean of innovation distribution
#' @param alpha Autoregressive parameter (0 < alpha < 1)
#' @param p Lag order for INAR process
#' @param family Innovation distribution ("pois" or "nbinom")
#' @param phi Dispersion parameter for Negative Binomial (required if family = "nbinom")
#' @return List containing:
#' \itemize{
#'   \item plots - Combined visualization of TS and ACF
#'   \item ts_data - Time series values after burn-in
#'   \item acf_data - ACF values
#'   \item innovations - Full innovation series
#'   \item parameters - Theoretical vs observed parameters
#' }
#' @export
simulate_inar_analysis <- function(
    t = 1000, burn_in = 100, lambda = 5, 
    alpha = 0.5, p = 1, family = "pois",
    phi = NULL) {
  
  # Input validation
  if (alpha <= 0 | alpha >= 1) stop("alpha must be between 0 and 1")
  if (burn_in >= t) stop("burn_in must be smaller than total time points")
  if (family == "nbinom" && (is.null(phi) || phi <= 0)) {
    stop("phi must be provided and positive for Negative Binomial innovations")
  }
  
  t <- t + burn_in + p  # Add burn-in and lag order to total time points
  
  # Generate innovations
  innov <- switch(family,
                  "pois" = rpois(t, lambda = lambda),
                  "nbinom" = rnbinom(t, mu = lambda, size = phi),
                  stop("Invalid innovation distribution"))
  
  # Simulate INAR process
  x <- sim_inar_p_ast(t, alpha, innov, p)
  
  # Remove burn-in period
  x <- x[(burn_in + 1):t]
  kept_innov <- innov[(burn_in + 1 - p):t]  # Innovations used in final series
  
  # Create dataframes
  ts_df <- data.frame(time = seq_along(x), value = x)
  acf_df <- data.frame(lag = 0:30, acf = acf(x, lag.max = 30, plot = FALSE)$acf[1:31])
  
  # Calculate statistics
  mu_x <- lambda / (1 - alpha)
  mu_x_obs <- mean(x)
  alpha_obs <- acf(x, plot = FALSE)$acf[p + 1]
  lambda_obs <- mu_x_obs * (1 - alpha_obs)
  
  # Calculate dispersion parameter (if applicable)
  if (family == "nbinom") {
    var_innov <- var(kept_innov)
    phi_obs <- ifelse(var_innov > lambda, 
                      lambda^2 / (var_innov - lambda), 
                      NA)
  }
  
  # Parameter matrix - FIXED VERSION
  if (family == "nbinom") {
    params <- rbind(
      Theoretical = c(mu_x, alpha, lambda, phi),
      Observed = c(mu_x_obs, alpha_obs, lambda_obs, phi_obs)
    )
    colnames(params) <- c("mu", "alpha", "lambda", "phi")
  } else {
    params <- rbind(
      Theoretical = c(mu_x, alpha, lambda),
      Observed = c(mu_x_obs, alpha_obs, lambda_obs)
    )
    colnames(params) <- c("mu", "alpha", "lambda")
  }
  
  family_name <- ifelse(family == "pois", "Poisson", "Negative Binomial")
  # Create plots
  ts_plot <- ggplot(ts_df, aes(time, value)) +
    geom_line(color = "steelblue", alpha = 0.7) +
    geom_hline(yintercept = c(mu_x, mu_x_obs), 
               linetype = c("dashed", "dotted"), color = c("red", "darkgreen")) +
    labs(
      title = paste0("INAR(", p, "*) Process with ", family_name, " Innovations"),
      y = "Value", x = "Time") +
    theme_minimal()
  
  acf_plot <- ggplot(acf_df, aes(lag, acf)) +
    geom_col(width = 0.1, fill = "steelblue") +
    geom_hline(
      yintercept = qnorm(c(0.025, 0.975))/sqrt(length(x)), 
      linetype = "dashed", color = "red") +
    scale_x_continuous(breaks = c(0, 5, 10, 15, 20, 25, 30)) +
    labs(
      title = "Autocorrelation Function", y = "ACF", x = "Lag"
    ) +
    theme_minimal()
  
  # Return results
  list(
    plots = ts_plot / acf_plot,
    ts_data = ts_df,
    acf_data = acf_df,
    innovations = kept_innov,
    parameters = params
  )
}

theme_set(theme_bw(base_family = "Helvetica") + 
            theme(
              plot.title = element_text(size = 24, face = "bold"),        
              plot.subtitle = element_text(size = 20),                    
              axis.title = element_text(size = 20),                       
              axis.text = element_text(size = 16),                        
              legend.title = element_text(size = 18),                     
              legend.text = element_text(size = 16),                      
              legend.position = "bottom",                                 
              strip.text = element_text(size = 18, face = "bold"),
              plot.margin = margin(t = 10, r = 40, b = 10, l = 10)
            )
)

#' @title Plot INAR Time Series
#' @description Creates a time series plot of all simulated series
#' @param sim_result List containing:
#' \itemize{
#'   \item data - Matrix of time series (n x t)
#' }
#' @return ggplot object showing all time series in a single color
#' @examples
#' params <- list(
#'   list(p = 1, alpha = 0.5, lambda = 5, family = "pois", phi = NULL)
#' )
#' sim_data <- sim_mixture_inar_p_ast(50, 30, c(1), params)
#' plot_inar_ts(sim_data)
#' @export
#' @importFrom ggplot2 ggplot aes geom_line labs theme_bw
#' @importFrom tidyr pivot_longer
plot_inar_ts <- function(sim_result, x_lim = c(0,0)) {
  
  # Get total number of the time series
  n_series <- nrow(sim_result)
  
  # Convert matrix to long format
  df <- as.data.frame(sim_result) %>%
    mutate(Series = factor(row.names(.))) %>%
    pivot_longer(cols = -Series, 
                 names_to = "Time", 
                 values_to = "Value") %>%
    mutate(Time = as.Date(Time, format = "%Y-%m-%d")) # this changed

    
  # Create plot with single color
  ggplot(df, aes(x = Time, y = Value, group = Series)) +
    geom_line(color = "steelblue", alpha = 0.4) +
    scale_x_date(expand = c(0,0)) +
    coord_cartesian(xlim = as.Date(x_lim),
                    ylim = c(0,max(sim_result)+1)) + 
    labs(title = "Average Daily Deaths for every 10,000 Inhabitants",
         subtitle = paste(n_series, "Counties"),
         x = "Time Point",
         y = "Observed Value") +
    theme_bw()
}

plot_inar_ts_membership <- function(sim_result, clusters, x_lim = c(0, 0)) {
  
  # Ensure input is in correct format
  n_series <- nrow(sim_result)
  
  # Reshape and add cluster labels
  df <- as.data.frame(sim_result) %>%
    mutate(
      Series = factor(row.names(.)),
      Cluster = factor(clusters)  # works for any number of clusters
    ) %>%
    tidyr::pivot_longer(cols = -c(Series, Cluster),
                        names_to = "Time",
                        values_to = "Value") %>%
    mutate(Time = as.Date(Time, format = "%Y-%m-%d"))
  
  # Plot
  ggplot(df, aes(x = Time, y = Value, group = Series, color = Cluster)) +
    geom_line(alpha = 0.4) +
    scale_x_date(expand = c(0, 0)) +
    coord_cartesian(xlim = as.Date(x_lim),
                    ylim = c(0, max(sim_result, na.rm = TRUE) + 1)) +
    labs(
      title = "Time Series by Cluster",
      subtitle = paste("All", n_series, "Series Colored by Cluster"),
      x = "Date",
      y = "Observed Value"
    ) +
    scale_color_viridis_d(name = "Cluster") +  # Adapts automatically
    theme_bw() +
    theme(legend.position = "bottom")
}

#' @title Plot INAR Time Series by Cluster
#' @description Creates a time series plot colored by cluster membership
#' @param sim_result List containing:
#' \itemize{
#'   \item data - Matrix of time series (n x t)
#'   \item z - Cluster assignments vector
#' }
#' @return ggplot object showing time series colored by cluster
#' @examples
#' params <- list(
#'   list(p = 1, alpha = 0.5, lambda = 5, family = "pois", phi = NULL),
#'   list(p = 2, alpha = 0.3, lambda = 8, family = "nbinom", phi = 2)
#' )
#' sim_data <- sim_mixture_inar_p_ast(50, 30, c(0.4, 0.6), params)
#' plot_inar_clusters(sim_data)
#' @export
plot_inar_clusters_elijah <- function(sim_result, clusters, x_lim = c(0,0)) {
  
  # Get total number of the time series
  n_series <- nrow(sim_result)
  
  # Convert matrix to data frame and add cluster information
  df <- as.data.frame(sim_result) %>%
    mutate(Series = factor(row.names(.)),
           Cluster = factor(clusters)) %>%
    pivot_longer(cols = -c(Series, Cluster), 
                 names_to = "Time", 
                 values_to = "Value") %>%
    mutate(Time = as.Date(Time, format = "%Y-%m-%d")) #this changed
  
  # Summarize: average daily value for each cluster
  df_avg <- df %>%
    group_by(Time, Cluster) %>%
    summarize(Avg_Value = mean(Value, na.rm = TRUE), .groups = "drop")
    
  # Create plot
  p <- ggplot(df_avg, aes(x = Time, y = Avg_Value)) +
    geom_line(size = 1.2) +
    facet_wrap(~ Cluster, scales = "free")+
    scale_x_date(expand = c(0, 0)) +
    coord_cartesian(xlim = as.Date(x_lim),
                    ylim = c(0,1.1)) + 
    scale_color_brewer(palette = "Set2") +
    labs(title = "Average Daily Deaths by Cluster for every 10,000 Inhabitants",
         subtitle = paste(n_series, "time series (clustered)"),
         x = "Date",
         y = "Average Deaths") +
    theme_bw() +
    theme(legend.position = "bottom")
  
  
  return(p)
}

plot_inar_clusters <- function(sim_result, clusters, x_lim = c(0,0)) {
  
  # Define consistent cluster colors
  cluster_colors <- c(
    "1" = "#3b0f70",  # purple
    "2" = "#1a9850",  # green
    "4" = "#fee08b"   # yellow
  )
  
  # Get total number of the time series
  n_series <- nrow(sim_result)
  
  # Convert matrix to data frame and add cluster information
  df <- as.data.frame(sim_result) %>%
    mutate(Series = factor(row.names(.)),
           Cluster = factor(clusters)) %>%
    pivot_longer(cols = -c(Series, Cluster), 
                 names_to = "Time", 
                 values_to = "Value") %>%
    mutate(Time = as.Date(Time, format = "%Y-%m-%d"))
  
  # Summarize: average daily value for each cluster
  df_avg <- df %>%
    group_by(Time, Cluster) %>%
    summarize(Avg_Value = mean(Value, na.rm = TRUE), .groups = "drop")
  
  # Create plot
  p <- ggplot(df_avg, aes(x = Time, y = Avg_Value, color = Cluster)) +
    geom_line(size = 1.2) +
    scale_x_date(expand = c(0, 0)) +
    coord_cartesian(xlim = as.Date(x_lim),
                    ylim = c(0, 1.5)) + 
    scale_color_manual(values = cluster_colors, name = "Cluster") +
    labs(title = "Average Daily Deaths by Cluster for every 10,000 Inhabitants",
         subtitle = paste(n_series, "counties (clustered)"),
         x = "Date",
         y = "Average Deaths") +
    theme_bw() +
    theme(legend.position = "bottom")
  
  return(p)
}

plot_inar_ts_and_cluster_averages <- function(sim_result, clusters){
  # Define consistent cluster colors
  cluster_colors <- c(
    "Cluster 1" = "#3b0f70",  # Purple
    "Cluster 2" = "#1a9850",  # Green
    "Cluster 3" = "#f6c55f"   # Yellow-orange
  )
  
  # Recode cluster factors
  clusters_recoded <- recode(as.character(clusters),
                             "1" = "Cluster 1",
                             "2" = "Cluster 2",
                             "4" = "Cluster 3")
  
  # Major COVID-19 events in the USA
  covid_events <- data.frame(
    date = as.Date(c("2020-03-13",  # National Emergency Declared
                     "2020-12-11",  # FDA approves Pfizer vaccine
                     "2021-01-20",  # Biden Inauguration
                     "2021-04-19",  # All adults eligible for vaccines
                     "2021-12-01")), # Omicron variant detected in USA
    event = c("National Emergency", 
              "Pfizer Vaccine Approval",
              "Biden Inauguration",
              "All Adults Vaccine Eligible",
              "Omicron Variant Detected")
  )
  
  # Reshape and add cluster labels
  df <- as.data.frame(sim_result) %>%
    mutate(
      Series = factor(row.names(.)),
      Cluster = factor(clusters_recoded, levels = c("Cluster 1", "Cluster 2", "Cluster 3"))
    ) %>%
    pivot_longer(cols = -c(Series, Cluster),
                 names_to = "Time",
                 values_to = "Value") %>%
    mutate(Time = as.Date(Time, format = "%Y-%m-%d"))
  
  # Compute cluster averages
  df_avg <- df %>%
    group_by(Time, Cluster) %>%
    summarize(Avg_Value = mean(Value, na.rm = TRUE), .groups = "drop")
  
  ## Set legend order explicitly
  legend_levels <- c("Cluster 1", "Cluster 2", "Cluster 3")
  
  # Arrange data to plot Cluster 3 first (behind)
  # Reorder factor to send Cluster 3 to the back
  df_avg$Cluster <- factor(df_avg$Cluster, levels = c("Cluster 3", "Cluster 1", "Cluster 2"))
  
  # Plot
  ggplot() +
    geom_line(data = df_avg, aes(x = Time, y = Avg_Value, color = Cluster), size = 0.8) +
    geom_vline(data = covid_events, aes(xintercept = date), linetype = "dashed", color = "red") +
    geom_text(data = covid_events, aes(x = date, label = event), y = Inf, angle = 90, 
              vjust = 1.25, hjust = 1.1, size = 5, color = "red") +
    scale_x_date(
      expand = c(0, 0),
      date_labels = "%b %Y", date_breaks = "6 months"
    ) +
    labs(
      title = "Time Series of Cluster Profiles",
      subtitle = "Deaths per 10,000 People",
      x = "Date",
      y = "Daily COVID-19 Deaths"
    ) +
    scale_color_manual(values = cluster_colors, name = "Clusters:", breaks = legend_levels)
}


#' @title Plot ACF Boxplots for INAR Mixture Simulations
#' @description Creates boxplots of autocorrelation values across multiple lags for simulated INAR data
#' @param sim_data Output from sim_mixture_inar_p_ast containing:
#' \itemize{
#'   \item data - Matrix of simulated time series (n x t)
#'   \item z - True cluster assignments
#' }
#' @param max_lag Maximum lag to display (default = 20)
#' @return ggplot object showing boxplots of ACF values by lag
#' @examples
#' params <- list(
#'   list(p = 1, alpha = 0.5, lambda = 5, family = "pois", phi = NULL),
#'   list(p = 2, alpha = 0.3, lambda = 8, family = "nbinom", phi = 2)
#' )
#' sim_data <- sim_mixture_inar_p_ast(100, 50, c(0.4, 0.6), params)
#' plot_inar_acf(sim_data)
#' @export
#' @importFrom ggplot2 ggplot aes geom_boxplot scale_x_continuous labs theme_bw
#' @importFrom stats acf
plot_inar_acf <- function(sim_data, max_lag = 20) {
  # Extract time series data
  data_matrix <- sim_data
  n_series <- nrow(data_matrix)
  
  # Initialize data frame to store ACF values
  acf_results <- vector("list", length = n_series)
  
  # Compute ACF values for each series
  for (i in seq_len(n_series)) {
    ts <- data_matrix[i, ]
    
    # Compute ACF values up to max_lag
    acf_vals <- acf(ts, lag.max = max_lag, plot = FALSE)$acf[-1] # exclude lag 0
    
    # Skip series if ACF is undefined (constant series)
    if (all(is.na(acf_vals))) next
    
    # Store computed ACF values
    acf_results[[i]] <- data.frame(
      Lag = seq_len(max_lag),
      ACF = as.numeric(acf_vals),
      Series = i
    )
  }
  
  # Combine all series into a single data frame
  acf_data <- do.call(rbind, acf_results)
  
  # Plot ACF boxplots
  ggplot(acf_data, aes(x = factor(Lag), y = ACF)) +
    geom_boxplot(fill = "skyblue", color = "darkblue") +
    scale_x_discrete(breaks = seq_len(max_lag)) +
    labs(
      title = "Distribution of Autocorrelation by Lag",
      subtitle = paste(n_series, "time series"),
      x = "Lag",
      y = "Autocorrelation"
    ) 
}


#' @title Plot Mean-Variance Relationship for Count Time Series
#' @description Creates a scatter plot of mean vs variance with x=y reference line
#' @param sim_data Output from sim_mixture_inar_p_ast containing:
#' \itemize{
#'   \item data - Matrix of simulated time series (n x t)
#' }
#' @return ggplot object showing mean-variance relationship
#' @examples
#' params <- list(
#'   list(p = 1, alpha = 0.5, lambda = 5, family = "pois", phi = NULL),
#'   list(p = 2, alpha = 0.3, lambda = 8, family = "nbinom", phi = 2)
#' )
#' sim_data <- sim_mixture_inar_p_ast(100, 50, c(0.4, 0.6), params)
#' plot_mean_variance(sim_data)
#' @export
#' @importFrom ggplot2 ggplot aes geom_point geom_abline labs theme_bw
plot_mean_variance <- function(sim_data, membership) {
  # Define cluster colors
  cluster_colors <- c(
    "Cluster 1" = "#3b0f70",  # Purple
    "Cluster 2" = "#1a9850",  # Green
    "Cluster 3" = "#f6c55f"   # Yellow-orange
  )
  
  # Recode cluster factors
  membership_recode <- recode(as.character(membership),
                              "1" = "Cluster 1",
                              "2" = "Cluster 2",
                              "4" = "Cluster 3")
  
  # Compute mean and variance for each time series
  plot_data <- data.frame(
    Mean = rowMeans(sim_data, na.rm = TRUE),
    Variance = apply(sim_data, 1, var, na.rm = TRUE),
    Cluster = factor(membership_recode, levels = c("Cluster 1", "Cluster 2", "Cluster 3"))
  ) %>%
    # Arrange data to plot Cluster 1 on top
    mutate(ordering = if_else(Cluster == "Cluster 1", 1, 0)) %>%
    arrange(ordering) %>%
    select(-ordering)
  
  # Generate plot
  ggplot(plot_data, aes(x = Mean, y = Variance, color = Cluster)) +
    geom_point(alpha = 0.9, size = 2.5) +
    geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
    labs(
      title = "Mean-Variance Relationship",
      subtitle = "Red line represents equidispersion (Mean = Variance)",
      x = "Mean of Time Series",
      y = "Variance of Time Series"
    ) +
    scale_x_continuous(limits = c(0, max(plot_data$Mean))) +
    scale_y_continuous(limits = c(0, max(plot_data$Variance))) +
    scale_color_manual(values = cluster_colors, name = "Clusters: ")
}


plot_mean_variance_simple <- function(sim_data) {
  # Compute mean and variance for each time series
  plot_data <- data.frame(
    Mean = rowMeans(sim_data, na.rm = TRUE),
    Variance = apply(sim_data, 1, var, na.rm = TRUE)
  )
  
  # Generate plot
  ggplot(plot_data, aes(x = Mean, y = Variance)) +
    geom_point(alpha = 0.9, size = 2.5, color = "#3b0f70") +  # Solid purple for all
    geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
    labs(
      title = "Mean-Variance Relationship",
      subtitle = "Red line represents equidispersion (Mean = Variance)",
      x = "Mean of Time Series",
      y = "Variance of Time Series"
    ) +
    scale_x_continuous(limits = c(0, max(plot_data$Mean))) +
    scale_y_continuous(limits = c(0, max(plot_data$Variance))) +
    theme_minimal(base_family = "Helvetica")
}
