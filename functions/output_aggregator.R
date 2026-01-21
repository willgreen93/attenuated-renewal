setwd("~/attenuated-renewal")

# Fast dependencies
library(dplyr)
library(tibble)
library(matrixStats)  # for colQuantiles (much faster than quantile in groups)

#files <- sorted_files_sim[1:2000]
files_sim <- list.files("fits/sim/fit_sim3", full.names = TRUE)

info_sim <- file.info(files_sim)
files_sim <- files_sim[!info_sim$isdir]

sorted_files_sim <- files_sim[order(info_sim$mtime[!info_sim$isdir], decreasing = FALSE)]

files <- sorted_files_sim


# Preallocate result holders
plot_list     <- vector("list", length(files))
lastweek_list <- vector("list", length(files))

# Thin every 10th sample (same logic as your filter(sample %% 10 == 0))
thin_every <- 10

for (ii in 1675:1750) {
  print(ii)
  i <- files[ii]
  # message(ii, "/", length(files), "  ", basename(i))  # optional light logging
  
  p <- proc.time()
  fit_output <- readRDS(i)
  (proc.time()-p)["elapsed"]
  
  N  <- fit_output$input_list$N
  h  <- fit_output$input_list$h
  Nd <- N + h
  
  # True infections vector (length N + h) and totals
  true_vec         <- fit_output$simulation$incidence$incidence[seq_len(Nd)]
  total_infection  <- sum(true_vec)
  
  # Sum across strata (keep [samples x days]) using apply margins (keeps your original intent)
  # infections: [samples x days] (keep first N cols)
  inf_m  <- apply(fit_output$infections, c(1, 2), sum)[, seq_len(N), drop = FALSE]
  # forecast:   [samples x h]
  fore_m <- apply(fit_output$forecast,   c(1, 2), sum)
  
  # Combine past + forecast → [samples x (N+h)]
  pred_mat <- cbind(inf_m, fore_m)
  
  # Thin samples by row index (same as filter(sample %% 10 == 0))
  nsamp    <- nrow(pred_mat)
  thin_idx <- which(seq_len(nsamp) %% thin_every == 0)
  pred_thin <- pred_mat[thin_idx, , drop = FALSE]
  
  # ---------- Fast 'inf_plot' (per-day quantiles across samples) ----------
  # Compute column quantiles (lower/median/upper) across rows (samples)
  lower  <- matrixStats::colQuantiles(pred_thin, probs = 0.025, na.rm = TRUE)
  median <- matrixStats::colQuantiles(pred_thin, probs = 0.50,  na.rm = TRUE)
  upper  <- matrixStats::colQuantiles(pred_thin, probs = 0.975, na.rm = TRUE)
  
  # Add params once, then build the 'simul' string from those cols
  pl <- fit_output$simulation$param_list  # named list of scalar params (assortativity, spatial, etc.)
  
  inf_plot <- tibble(
    day          = seq_len(Nd),
    true_value   = true_vec,
    lower        = as.numeric(lower),
    median       = as.numeric(median),
    upper        = as.numeric(upper),
    scenario     = fit_output$tag,
    cutoff       = N,
    model        = fit_output$model_name,
    epi_phase    = fit_output$input_list$epi_phase[1],
    total_infection = total_infection
  ) %>%
    mutate(!!!pl) %>%  # splice parameters into columns (cheap)
    mutate(cutoff2   = factor(as.character(cutoff)),
           epi_phase = as.numeric(epi_phase),
           scenario  = as.character(scenario))
  
  # ---------- Fast 'inf_last_week' (sum over days N+22 : N+28) ----------
  week_idx <- intersect((N + 22):(N + 28), seq_len(Nd))  # guard if h < 28
  # Per-sample sums for last week (rows = samples)
  pred_week_sum <- if (length(week_idx)) {
    rowSums(pred_thin[, week_idx, drop = FALSE])
  } else {
    rep(NA_real_, nrow(pred_thin))
  }
  true_week_sum <- if (length(week_idx)) sum(true_vec[week_idx]) else NA_real_
  
  inf_last_week <- tibble(
    scenario        = fit_output$tag,
    sample          = thin_idx,  # preserve original sample ids
    model           = fit_output$model_name,
    cutoff          = N,
    epi_phase       = fit_output$input_list$epi_phase[1],
    total_infection = total_infection,
    predicted      = as.numeric(pred_week_sum),
    true_value      = true_week_sum
  ) %>%
    mutate(epi_phase = as.numeric(epi_phase),
           scenario = as.character(scenario)) %>%
    mutate(!!!pl)
  
  # Accumulate (no rbind inside loop)
  plot_list[[ii]]     <- inf_plot
  lastweek_list[[ii]] <- inf_last_week
}

# Bind once at the end (fast)
saveRDS(plot_list, "fits/outputs/plot_list.rds")
saveRDS(lastweek_list, "fits/outputs/lastweek_list.rds")

skeleton_plot_sim      <- dplyr::bind_rows(plot_list)     %>% rename(place=scenario) %>% mutate(sr="SIM") %>% select(day, true_value, cutoff, model, epi_phase, total_infection, lower, median, upper, cutoff2, place, sr)
skeleton_last_week_sim <- dplyr::bind_rows(lastweek_list) %>% rename(place=scenario) %>% mutate(sr="SIM") %>% select(sample, model, cutoff, epi_phase, total_infection, predicted, true_value, place, sr)

saveRDS(skeleton_plot_sim, "fits/outputs/skeleton_plot_sim_all2.rds")
saveRDS(skeleton_last_week_sim, "fits/outputs/skeleton_last_week_sim_all2.rds")


files_cit <- list.files("fits/cities/", full.names = TRUE) 

info_cit <- file.info(files_cit)
files_cit <- files_cit[!info_cit$isdir]

plot_list_cit     <- vector("list", length(files_cit))
lastweek_list_cit <- vector("list", length(files_cit))

for (ii in seq_along(files_cit)) {
  print(ii)
  i <- files_cit[ii]
  # message(ii, "/", length(files), "  ", basename(i))  # optional light logging
  
  fit_output <- readRDS(i)
  N  <- fit_output$input_list$N
  h  <- fit_output$input_list$h
  Nd <- N + h
  
  # True infections vector (length N + h) and totals
  true_vec         <- fit_output$incidence_full[seq_len(Nd)]
  total_infection  <- sum(true_vec)
  
  # Sum across strata (keep [samples x days]) using apply margins (keeps your original intent)
  # infections: [samples x days] (keep first N cols)
  inf_m  <- apply(fit_output$infections, c(1, 2), sum)[, seq_len(N), drop = FALSE]
  # forecast:   [samples x h]
  fore_m <- apply(fit_output$forecast,   c(1, 2), sum)
  
  # Combine past + forecast → [samples x (N+h)]
  pred_mat <- cbind(inf_m, fore_m)
  
  # Thin samples by row index (same as filter(sample %% 10 == 0))
  nsamp    <- nrow(pred_mat)
  thin_idx <- which(seq_len(nsamp) %% thin_every == 0)
  pred_thin <- pred_mat[thin_idx, , drop = FALSE]
  
  # ---------- Fast 'inf_plot' (per-day quantiles across samples) ----------
  # Compute column quantiles (lower/median/upper) across rows (samples)
  lower  <- matrixStats::colQuantiles(pred_thin, probs = 0.025, na.rm = TRUE)
  median <- matrixStats::colQuantiles(pred_thin, probs = 0.50,  na.rm = TRUE)
  upper  <- matrixStats::colQuantiles(pred_thin, probs = 0.975, na.rm = TRUE)
  
  # Add params once, then build the 'simul' string from those cols
  #pl <- fit_output$simulation$param_list  # named list of scalar params (assortativity, spatial, etc.)
  
  inf_plot <- tibble(
    day          = seq_len(Nd),
    true_value   = true_vec,
    lower        = as.numeric(lower),
    median       = as.numeric(median),
    upper        = as.numeric(upper),
    cutoff2      = factor(as.character(cutoff)),
    place        = fit_output$tag,
    cutoff       = N,
    model        = fit_output$model_name,
    epi_phase    = as.numeric(fit_output$input_list$epi_phase[1]),
    total_infection = total_infection
  ) 
  
  # ---------- Fast 'inf_last_week' (sum over days N+22 : N+28) ----------
  week_idx <- intersect((N + 22):(N + 28), seq_len(Nd))  # guard if h < 28
  # Per-sample sums for last week (rows = samples)
  pred_week_sum <- if (length(week_idx)) {
    rowSums(pred_thin[, week_idx, drop = FALSE])
  } else {
    rep(NA_real_, nrow(pred_thin))
  }
  true_week_sum <- if (length(week_idx)) sum(true_vec[week_idx]) else NA_real_
  
  inf_last_week <- tibble(
    place           = fit_output$tag,
    sample          = thin_idx,  # preserve original sample ids
    model           = fit_output$model_name,
    cutoff          = N,
    epi_phase       = fit_output$input_list$epi_phase[1],
    total_infection = total_infection,
    predicted       = as.numeric(pred_week_sum),
    true_value      = true_week_sum
  ) %>%
    mutate(epi_phase = as.numeric(epi_phase)) #%>%
  #mutate(!!!pl)
  
  # Accumulate (no rbind inside loop)
  plot_list_cit[[ii]]     <- inf_plot
  lastweek_list_cit[[ii]] <- inf_last_week
}

skeleton_plot_cit      <- dplyr::bind_rows(plot_list_cit)     %>% mutate(sr="CIT") %>% select(day, true_value, cutoff, model, epi_phase, total_infection, lower, median, upper, cutoff2, place, sr)
skeleton_last_week_cit <- dplyr::bind_rows(lastweek_list_cit) %>% mutate(sr="CIT") %>% select(sample, model, cutoff, epi_phase, total_infection, predicted, true_value, place, sr)

saveRDS(skeleton_plot_cit, "fits/outputs/skeleton_plot_cit.rds")
saveRDS(skeleton_last_week_cit, "fits/outputs/skeleton_last_week_cit.rds")

files_homo <- list.files("fits/sim/fit_homo", full.names = TRUE) 

plot_list_homo     <- vector("list", length(files_homo))
lastweek_list_homo <- vector("list", length(files_homo))

for (ii in seq_along(files_homo)) {
  i <- files_homo[ii]
  print(i)
  # message(ii, "/", length(files), "  ", basename(i))  # optional light logging
  
  fit_output <- readRDS(i)
  N  <- fit_output$input_list$N
  h  <- fit_output$input_list$h
  Nd <- N + h
  
  # True infections vector (length N + h) and totals
  true_vec         <- fit_output$simulation$incidence$incidence[seq_len(Nd)]
  total_infection  <- sum(true_vec)
  
  # Sum across strata (keep [samples x days]) using apply margins (keeps your original intent)
  # infections: [samples x days] (keep first N cols)
  inf_m  <- apply(fit_output$infections, c(1, 2), sum)[, seq_len(N), drop = FALSE]
  # forecast:   [samples x h]
  fore_m <- apply(fit_output$forecast,   c(1, 2), sum)
  
  # Combine past + forecast → [samples x (N+h)]
  pred_mat <- cbind(inf_m, fore_m)
  
  # Thin samples by row index (same as filter(sample %% 10 == 0))
  nsamp    <- nrow(pred_mat)
  thin_idx <- which(seq_len(nsamp) %% thin_every == 0)
  pred_thin <- pred_mat[thin_idx, , drop = FALSE]
  
  # ---------- Fast 'inf_plot' (per-day quantiles across samples) ----------
  # Compute column quantiles (lower/median/upper) across rows (samples)
  lower  <- matrixStats::colQuantiles(pred_thin, probs = 0.025, na.rm = TRUE)
  median <- matrixStats::colQuantiles(pred_thin, probs = 0.50,  na.rm = TRUE)
  upper  <- matrixStats::colQuantiles(pred_thin, probs = 0.975, na.rm = TRUE)
  
  # Add params once, then build the 'simul' string from those cols
  pl <- fit_output$simulation$param_list  # named list of scalar params (assortativity, spatial, etc.)
  
  inf_plot <- tibble(
    day          = seq_len(Nd),
    true_value   = true_vec,
    lower        = as.numeric(lower),
    median       = as.numeric(median),
    upper        = as.numeric(upper),
    scenario     = fit_output$tag,
    cutoff       = N,
    model        = fit_output$model_name,
    epi_phase    = fit_output$input_list$epi_phase[1],
    total_infection = total_infection
  ) %>%
    mutate(!!!pl) %>%  # splice parameters into columns (cheap)
    mutate(
      simul = paste0(
        cutoff, "_", model, "_",
        assortativity, "_", spatial, "_",
        isolation_time_mean, "_", isolation_time_sd, "_",
        round(ass_v_param, 2), "_", round(transmission_prob, 2), "_",
        round(directionality, 2), "_", incubation_period
      ),
      cutoff2  = factor(as.character(cutoff)),
      epi_phase = as.factor(epi_phase)
    )
  # ---------- Fast 'inf_last_week' (sum over days N+22 : N+28) ----------
  week_idx <- intersect((N + 22):(N + 28), seq_len(Nd))  # guard if h < 28
  # Per-sample sums for last week (rows = samples)
  pred_week_sum <- if (length(week_idx)) {
    rowSums(pred_thin[, week_idx, drop = FALSE])
  } else {
    rep(NA_real_, nrow(pred_thin))
  }
  true_week_sum <- if (length(week_idx)) sum(true_vec[week_idx]) else NA_real_
  
  inf_last_week <- tibble(
    scenario        = fit_output$tag,
    sample          = thin_idx,  # preserve original sample ids
    model           = fit_output$model_name,
    cutoff          = N,
    epi_phase       = fit_output$input_list$epi_phase[1],
    total_infection = total_infection,
    predicted      = as.numeric(pred_week_sum),
    true_value      = true_week_sum
  ) %>%
    mutate(epi_phase = as.factor(epi_phase)) %>%
    mutate(!!!pl)
  
  # Accumulate (no rbind inside loop)
  plot_list_homo[[ii]]     <- inf_plot
  lastweek_list_homo[[ii]] <- inf_last_week
  
  
}

saveRDS(dplyr::bind_rows(plot_list_homo)     %>% mutate(place=as.character(scenario+130), sr="HOMO", epi_phase=as.numeric(as.character(epi_phase))) %>% select(day, true_value, cutoff, model, epi_phase, total_infection, lower, median, upper, cutoff2, place, sr), "fits/outputs/skeleton_plot_homo.rds")      
saveRDS(dplyr::bind_rows(lastweek_list_homo) %>% mutate(place=as.character(scenario+130), sr="HOMO", epi_phase=as.numeric(as.character(epi_phase))) %>% select(sample, model, cutoff, epi_phase, total_infection, predicted, true_value, place, sr), "fits/outputs/skeleton_last_week_homo.rds") 


####

files_het <- list.files("fits/sim/fit_het", full.names = TRUE) 

plot_list_het     <- vector("list", length(files_het))
lastweek_list_het <- vector("list", length(files_het))

for (ii in seq_along(files_het)) {
  i <- files_het[ii]
  print(i)
  # message(ii, "/", length(files), "  ", basename(i))  # optional light logging
  
  fit_output <- readRDS(i)
  N  <- fit_output$input_list$N
  h  <- fit_output$input_list$h
  Nd <- N + h
  
  # True infections vector (length N + h) and totals
  true_vec         <- fit_output$simulation$incidence$incidence[seq_len(Nd)]
  total_infection  <- sum(true_vec)
  
  # Sum across strata (keep [samples x days]) using apply margins (keeps your original intent)
  # infections: [samples x days] (keep first N cols)
  inf_m  <- apply(fit_output$infections, c(1, 2), sum)[, seq_len(N), drop = FALSE]
  # forecast:   [samples x h]
  fore_m <- apply(fit_output$forecast,   c(1, 2), sum)
  
  # Combine past + forecast → [samples x (N+h)]
  pred_mat <- cbind(inf_m, fore_m)
  
  # Thin samples by row index (same as filter(sample %% 10 == 0))
  nsamp    <- nrow(pred_mat)
  thin_idx <- which(seq_len(nsamp) %% thin_every == 0)
  pred_thin <- pred_mat[thin_idx, , drop = FALSE]
  
  # ---------- Fast 'inf_plot' (per-day quantiles across samples) ----------
  # Compute column quantiles (lower/median/upper) across rows (samples)
  lower  <- matrixStats::colQuantiles(pred_thin, probs = 0.025, na.rm = TRUE)
  median <- matrixStats::colQuantiles(pred_thin, probs = 0.50,  na.rm = TRUE)
  upper  <- matrixStats::colQuantiles(pred_thin, probs = 0.975, na.rm = TRUE)
  
  # Add params once, then build the 'simul' string from those cols
  pl <- fit_output$simulation$param_list  # named list of scalar params (assortativity, spatial, etc.)
  
  inf_plot <- tibble(
    day          = seq_len(Nd),
    true_value   = true_vec,
    lower        = as.numeric(lower),
    median       = as.numeric(median),
    upper        = as.numeric(upper),
    scenario     = fit_output$tag,
    cutoff       = N,
    model        = fit_output$model_name,
    epi_phase    = fit_output$input_list$epi_phase[1],
    total_infection = total_infection
  ) %>%
    mutate(!!!pl) %>%  # splice parameters into columns (cheap)
    mutate(
      simul = paste0(
        cutoff, "_", model, "_",
        assortativity, "_", spatial, "_",
        isolation_time_mean, "_", isolation_time_sd, "_",
        round(ass_v_param, 2), "_", round(transmission_prob, 2), "_",
        round(directionality, 2), "_", incubation_period
      ),
      cutoff2  = factor(as.character(cutoff)),
      epi_phase = as.factor(epi_phase)
    )
  # ---------- Fast 'inf_last_week' (sum over days N+22 : N+28) ----------
  week_idx <- intersect((N + 22):(N + 28), seq_len(Nd))  # guard if h < 28
  # Per-sample sums for last week (rows = samples)
  pred_week_sum <- if (length(week_idx)) {
    rowSums(pred_thin[, week_idx, drop = FALSE])
  } else {
    rep(NA_real_, nrow(pred_thin))
  }
  true_week_sum <- if (length(week_idx)) sum(true_vec[week_idx]) else NA_real_
  
  inf_last_week <- tibble(
    scenario        = fit_output$tag,
    sample          = thin_idx,  # preserve original sample ids
    model           = fit_output$model_name,
    cutoff          = N,
    epi_phase       = fit_output$input_list$epi_phase[1],
    total_infection = total_infection,
    predicted      = as.numeric(pred_week_sum),
    true_value      = true_week_sum
  ) %>%
    mutate(epi_phase = as.factor(epi_phase)) %>%
    mutate(!!!pl)
  
  # Accumulate (no rbind inside loop)
  plot_list_het[[ii]]     <- inf_plot
  lastweek_list_het[[ii]] <- inf_last_week
}

saveRDS(dplyr::bind_rows(plot_list_het)     %>% mutate(place=as.character(scenario+100), sr="HET", epi_phase=as.numeric(as.character(epi_phase))) %>% select(day, true_value, cutoff, model, epi_phase, total_infection, lower, median, upper, cutoff2, place, sr), "fits/outputs/skeleton_plot_het.rds")      
saveRDS(dplyr::bind_rows(lastweek_list_het) %>% mutate(place=as.character(scenario+100), sr="HET", epi_phase=as.numeric(as.character(epi_phase))) %>% select(sample, model, cutoff, epi_phase, total_infection, predicted, true_value, place, sr), "fits/outputs/skeleton_last_week_het.rds") 

