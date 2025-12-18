#=============================================================================--
# 0. General Settings & Libraries  ----
#=============================================================================--

set.seed(20251202)

req <- c(
  "data.table", "dplyr", "tibble", "tidyr", "purrr",
  "igraph", "tictoc",
  "plm", "sandwich", "lmtest", "EconGeo"
)

to_install <- setdiff(req, rownames(installed.packages()))
if (length(to_install)) install.packages(to_install)

suppressPackageStartupMessages(invisible(lapply(req, library, character.only = TRUE)))

setwd(this.path::here())
getwd()

#=============================================================================--
# 1. Data Loading & Preprocessing  ----
#=============================================================================--

# 1.1 Regions (EU selection) 
regions <- data.table::fread(
  "../data/REGPAT_REGIONS.txt",
  sep = "|", header = TRUE, encoding = "UTF-8"
) %>%
  as_tibble() %>%
  rename_with(tolower)

eu_ctry <- c(
  "AT","BE","BG","CH","CY","CZ","DE","DK","EE","ES",
  "FI","FR","GB","GR","HR","HU","IE","IS","IT","LI","LT",
  "LU","LV","MT","NL","NO","PL","PT","RO","SE","SI","SK"
)

regions_eu <- regions %>% filter(ctry_code %in% eu_ctry)

eu_nuts3_codes <- regions_eu %>%
  distinct(reg_code) %>%
  pull(reg_code)

nuts3_to_nuts2 <- regions_eu %>%
  distinct(reg_code, up_reg_code)

rm(regions, regions_eu); gc()

# 1.2 IPC codes (global, later restricted to EU patents)
ipc_epo <- data.table::fread(
  "../data/202401_EPO_IPC.txt",
  sep = "|", header = TRUE, encoding = "UTF-8"
)

tech_df <- ipc_epo %>%
  as_tibble() %>%
  transmute(
    pat_id   = as.character(appln_id),
    app_year = as.integer(app_year),
    ipc4     = substr(IPC, 1, 4)
  ) %>%
  filter(!is.na(app_year), !is.na(ipc4), ipc4 != "") %>%
  distinct(pat_id, app_year, ipc4)

rm(ipc_epo); gc()

# 1.3 Inventors (EU only)
inv_epo <- data.table::fread(
  "../data/202401_EPO_Inv_reg.txt",
  sep = "|", header = TRUE, encoding = "UTF-8"
)

inv_epo_eu <- inv_epo %>%
  as_tibble() %>%
  filter(ctry_code %in% eu_ctry, reg_code %in% eu_nuts3_codes) %>%
  transmute(
    person_id = as.character(person_id),
    appln_id  = as.character(appln_id),
    reg_share,
    inv_share,
    reg_code,
    ctry_code
  )

rm(inv_epo); gc()

# 1.4 Merge year info into inventors
pat_years <- tech_df %>%
  select(pat_id, app_year) %>%
  distinct() %>%
  rename(appln_id = pat_id)

inv_epo_eu <- inv_epo_eu %>%
  inner_join(pat_years, by = "appln_id")

# 1.5 Restrict IPC universe to EU patents only for Z-scores
eu_pat_ids_all <- inv_epo_eu %>% distinct(appln_id) %>% pull(appln_id)

tech_df_eu <- tech_df %>%
  filter(pat_id %in% eu_pat_ids_all)

rm(tech_df); gc()

#=============================================================================--
# 2. Time Windows (5-year non-overlapping)  ----
#=============================================================================--

min_year_data <- min(tech_df_eu$app_year, na.rm = TRUE)
max_year_data <- max(tech_df_eu$app_year, na.rm = TRUE)

min_year <- max(1979L, min_year_data)
max_year <- min(2023L, max_year_data)

starts <- seq(min_year, max_year - 4L, by = 5L)
ends   <- pmin(starts + 4L, max_year)

windows <- tibble(
  period_id = seq_along(starts),
  start_y   = starts,
  end_y     = ends,
  window    = paste0(starts, "-", ends)
)

print(windows)

#=============================================================================--
# 3. IPC Z-scores + Patent-level atypicality per window (rolling: prev+current)  ----
#=============================================================================-- 

compute_ipc_window <- function(period_id, windows_df, tech_df_eu, inv_epo_eu) {
  
  w <- windows_df %>% filter(period_id == !!period_id)
  start_y <- w$start_y
  end_y   <- w$end_y
  w_label <- w$window
  
  cat(sprintf("IPC Z-scores for window %s (%d-%d)\n", w_label, start_y, end_y))
  
  # Estimation window: preceding 5 years + current (no future)
  real_min <- min(tech_df_eu$app_year, na.rm = TRUE)
  est_start <- max(start_y - 5L, real_min)
  
  est_codes <- tech_df_eu %>%
    filter(app_year >= est_start, app_year <= end_y) %>%
    distinct(pat_id, ipc4)
  
  if (nrow(est_codes) == 0L) return(list(pair_stats = tibble(), pat_flags = tibble()))
  
  pat_sizes <- est_codes %>% count(pat_id, name = "n_codes")
  est_codes <- est_codes %>% inner_join(pat_sizes %>% filter(n_codes >= 2L), by = "pat_id")
  
  N <- n_distinct(est_codes$pat_id)
  if (N < 2L) return(list(pair_stats = tibble(), pat_flags = tibble()))
  
  # Pair generation in estimation window
  est_pairs <- est_codes %>%
    inner_join(est_codes, by = "pat_id", suffix = c(".i", ".j")) %>%
    filter(ipc4.i < ipc4.j)
  
  pair_stats <- est_pairs %>%
    count(ipc4.i, ipc4.j, name = "O_obs") %>%
    rename(ipc_i = ipc4.i, ipc_j = ipc4.j)
  
  code_counts <- est_codes %>% count(ipc4, name = "n_i")
  
  pair_stats <- pair_stats %>%
    left_join(code_counts %>% rename(ipc_i = ipc4), by = "ipc_i") %>%
    left_join(code_counts %>% rename(ipc_j = ipc4, n_j = n_i), by = "ipc_j") %>%
    mutate(
      N_pat  = as.numeric(N),
      E_ij   = (as.numeric(n_i) * as.numeric(n_j)) / N_pat,
      sigma2 = E_ij * (1 - n_i / N_pat) * ((N_pat - n_j) / (N_pat - 1)),
      sigma  = sqrt(pmax(sigma2, 0)),
      Z_score = if_else(sigma > 0, (O_obs - E_ij) / sigma, NA_real_)
    ) %>%
    filter(!is.na(Z_score)) %>%
    mutate(
      period_id   = period_id,
      window      = w_label,
      is_atypical = (Z_score < 0)
    ) %>%
    select(ipc_i, ipc_j, O_obs, n_i, n_j, E_ij, sigma, Z_score,
           period_id, window, is_atypical)
  
  # Patent-level atypicality (target patents only)
  eu_pat_ids <- inv_epo_eu %>%
    filter(app_year >= start_y, app_year <= end_y) %>%
    distinct(appln_id) %>%
    pull(appln_id)
  
  if (length(eu_pat_ids) == 0L) return(list(pair_stats = pair_stats, pat_flags = tibble()))
  
  target_codes <- tech_df_eu %>%
    filter(pat_id %in% eu_pat_ids) %>%
    distinct(pat_id, ipc4)
  
  pat_sizes_t <- target_codes %>% count(pat_id, name = "n_codes")
  target_codes <- target_codes %>% inner_join(pat_sizes_t %>% filter(n_codes >= 2L), by = "pat_id")
  
  if (nrow(target_codes) == 0L) return(list(pair_stats = pair_stats, pat_flags = tibble()))
  
  target_pairs <- target_codes %>%
    inner_join(target_codes, by = "pat_id", suffix = c(".i", ".j")) %>%
    filter(ipc4.i < ipc4.j) %>%
    rename(ipc_i = ipc4.i, ipc_j = ipc4.j)
  
  pat_flags <- target_pairs %>%
    left_join(pair_stats %>% select(ipc_i, ipc_j, is_atypical), by = c("ipc_i", "ipc_j")) %>%
    mutate(is_atypical = replace_na(is_atypical, FALSE)) %>%
    group_by(pat_id) %>%
    summarise(
      n_pairs            = n(),
      n_atyp_pairs       = sum(is_atypical),
      is_atypical_patent = any(is_atypical),
      .groups = "drop"
    ) %>%
    mutate(period_id = period_id, window = w_label)
  
  list(pair_stats = pair_stats, pat_flags = pat_flags)
}

#=============================================================================--
# 4. Run IPC Z-scores + patent atypicality (parallel if unix)  ----
#=============================================================================--

n_cores <- max(1L, parallel::detectCores() - 3L)
cat(sprintf("Using %d cores.\n", n_cores))

tic("IPC Z-scores + patent atypicality")

ipc_results_list <-
  if (.Platform$OS.type == "unix") {
    parallel::mclapply(
      windows$period_id,
      compute_ipc_window,
      windows_df = windows,
      tech_df_eu = tech_df_eu,
      inv_epo_eu = inv_epo_eu,
      mc.cores   = n_cores
    )
  } else {
    lapply(
      windows$period_id,
      compute_ipc_window,
      windows_df = windows,
      tech_df_eu = tech_df_eu,
      inv_epo_eu = inv_epo_eu
    )
  }

ipc_pair_stats <- purrr::map(ipc_results_list, "pair_stats") %>% bind_rows()
pat_atyp_flags <- purrr::map(ipc_results_list, "pat_flags") %>% bind_rows()

toc()

data.table::fwrite(ipc_pair_stats, "ipc_pair_zscores_eu.csv")

#=============================================================================--
# 5. Build Inventor Panel per window (backboned places) with betweenness ----
#=============================================================================--

build_window_panel <- function(curr_period_id, windows_df, inv_epo_eu, pat_flags_df) {
  
  w <- windows_df %>% filter(period_id == !!curr_period_id)
  w_start <- w$start_y
  w_end   <- w$end_y
  w_label <- w$window
  
  cat(sprintf("Inventor panel for window %s (%d-%d)\n", w_label, w_start, w_end))
  
  window_invs <- inv_epo_eu %>%
    filter(app_year >= w_start, app_year <= w_end) %>%
    transmute(
      person_id,
      pat_id   = appln_id,
      reg_code,
      ctry_code
    ) %>%
    distinct()
  
  if (nrow(window_invs) == 0L) return(tibble())
  
  all_inventors <- window_invs %>% distinct(person_id) %>% pull(person_id)
  
  # Co-inventor edges (unique, unweighted)
  co_pairs <- window_invs %>%
    inner_join(window_invs, by = "pat_id", suffix = c(".i", ".j")) %>%
    filter(person_id.i < person_id.j) %>%
    distinct(person_id.i, person_id.j)
  
  # Place signatures (structural equivalence by identical neighbor sets)
  if (nrow(co_pairs) > 0L) {
    neighbors <- bind_rows(
      co_pairs %>% transmute(person_id = person_id.i, neighbor = person_id.j),
      co_pairs %>% transmute(person_id = person_id.j, neighbor = person_id.i)
    ) %>%
      distinct() %>%
      arrange(person_id, neighbor)
    
    place_map <- neighbors %>%
      group_by(person_id) %>%
      summarise(place_sig = paste(sort(unique(neighbor)), collapse = "|"), .groups = "drop")
  } else {
    place_map <- tibble(person_id = character(), place_sig = character())
  }
  
  isolates <- setdiff(all_inventors, place_map$person_id)
  if (length(isolates) > 0L) {
    place_map <- bind_rows(
      place_map,
      tibble(person_id = isolates, place_sig = NA_character_)
    )
  }
  
  place_map <- place_map %>%
    mutate(place_id = if_else(is.na(place_sig), NA_integer_, as.integer(factor(place_sig))))
  
  # Place graph + network stats mapped back to inventors
  if (nrow(co_pairs) > 0L) {
    
    place_map_noniso <- place_map %>% filter(!is.na(place_sig))
    
    co_pairs_places <- co_pairs %>%
      left_join(place_map_noniso %>% select(person_id, place_id),
                by = c("person_id.i" = "person_id")) %>%
      rename(place_i = place_id) %>%
      left_join(place_map_noniso %>% select(person_id, place_id),
                by = c("person_id.j" = "person_id")) %>%
      rename(place_j = place_id) %>%
      filter(!is.na(place_i), !is.na(place_j), place_i != place_j) %>%
      distinct(place_i, place_j)
    
    if (nrow(co_pairs_places) > 0L) {
      g_place <- graph_from_data_frame(
        co_pairs_places %>% transmute(from = as.character(place_i), to = as.character(place_j)),
        directed = FALSE
      )
    } else {
      g_place <- make_empty_graph()
    }
    
    all_places <- unique(place_map_noniso$place_id)
    missing_places <- setdiff(as.character(all_places), V(g_place)$name)
    if (length(missing_places) > 0L) {
      g_place <- add_vertices(g_place, nv = length(missing_places), name = missing_places)
    }
    
    if (vcount(g_place) > 0L) {
      deg_place  <- igraph::degree(g_place)
      cons_place <- igraph::constraint(g_place)
      btw_place  <- igraph::betweenness(g_place, directed = FALSE, normalized = TRUE)
      
      place_stats <- tibble(
        place_id    = as.integer(names(deg_place)),
        degree      = as.numeric(deg_place),
        constraint  = as.numeric(cons_place),
        betweenness = as.numeric(btw_place)
      )
      
      network_stats_noniso <- place_map_noniso %>%
        left_join(place_stats, by = "place_id") %>%
        transmute(person_id, degree, constraint, betweenness)
      
    } else {
      network_stats_noniso <- place_map_noniso %>%
        transmute(person_id, degree = 0, constraint = NA_real_, betweenness = 0)
    }
    
    isolates_vec <- place_map %>% filter(is.na(place_sig)) %>% pull(person_id)
    network_stats_iso <- tibble(person_id = isolates_vec, degree = 0, constraint = NA_real_, betweenness = 0)
    
    network_stats <- bind_rows(network_stats_noniso, network_stats_iso)
    
  } else {
    network_stats <- tibble(person_id = all_inventors, degree = 0, constraint = NA_real_, betweenness = 0)
  }
  
  # Patent atypicality -> inventor outcomes
  pat_flags_win <- pat_flags_df %>%
    filter(period_id == curr_period_id) %>%
    select(pat_id, is_atypical_patent)
  
  inv_period <- window_invs %>%
    distinct(person_id, pat_id) %>%
    left_join(pat_flags_win, by = "pat_id") %>%
    mutate(is_atypical_patent = replace_na(is_atypical_patent, FALSE)) %>%
    group_by(person_id) %>%
    summarise(
      n_pat      = n_distinct(pat_id),
      n_atyp     = sum(is_atypical_patent),
      share_atyp = if_else(n_pat > 0, mean(is_atypical_patent), NA_real_),
      .groups    = "drop"
    ) %>%
    mutate(period_id = curr_period_id)
  
  # Dominant region & country per inventor (window-local)
  region_info <- window_invs %>%
    count(person_id, reg_code, ctry_code, name = "n_obs") %>%
    arrange(person_id, desc(n_obs)) %>%
    group_by(person_id) %>%
    slice_head(n = 1) %>%
    ungroup() %>%
    select(person_id, reg_code, ctry_code)
  
  # Final panel slice
  tibble(person_id = all_inventors) %>%
    mutate(period_id = curr_period_id, window = w_label) %>%
    left_join(place_map %>% select(person_id, place_id, place_sig), by = "person_id") %>%
    left_join(network_stats, by = "person_id") %>%
    left_join(inv_period, by = c("person_id", "period_id")) %>%
    left_join(region_info, by = "person_id") %>%
    mutate(
      n_pat      = replace_na(n_pat, 0L),
      n_atyp     = replace_na(n_atyp, 0L),
      share_atyp = replace_na(share_atyp, 0)
    )
}

#=============================================================================--
# 6. Run panel construction (parallel if unix)  ----
#=============================================================================--

tic("Inventor panel construction")

panel_list <-
  if (.Platform$OS.type == "unix") {
    parallel::mclapply(
      windows$period_id,
      build_window_panel,
      windows_df   = windows,
      inv_epo_eu   = inv_epo_eu,
      pat_flags_df = pat_atyp_flags,
      mc.cores     = n_cores
    )
  } else {
    lapply(
      windows$period_id,
      build_window_panel,
      windows_df   = windows,
      inv_epo_eu   = inv_epo_eu,
      pat_flags_df = pat_atyp_flags
    )
  }

final_panel_df <- bind_rows(panel_list)

toc()

data.table::fwrite(final_panel_df, "final_inventor_panel_eu.csv")
cat("Panel construction complete.\n")
print(head(final_panel_df))

#=============================================================================--
# 7. Lightweight summaries (ONE VERSION EACH, no duplicates) ----
#=============================================================================--

# 7.1 IPC Z-score summary
cat("Summary of IPC pair Z-scores:\n")
print(summary(ipc_pair_stats$Z_score))

# 7.2 Network summary per window (inventor graph, not place graph)
summarize_window_network <- function(curr_period_id, windows_df, inv_epo_eu, pat_flags_df) {
  
  w <- windows_df %>% filter(period_id == !!curr_period_id)
  w_start <- w$start_y
  w_end   <- w$end_y
  w_label <- w$window
  
  cat(sprintf("Network summary for window %s (%d-%d)\n", w_label, w_start, w_end))
  
  window_invs <- inv_epo_eu %>%
    filter(app_year >= w_start, app_year <= w_end) %>%
    transmute(person_id, pat_id = appln_id, ctry_code) %>%
    distinct()
  
  if (nrow(window_invs) == 0L) {
    return(tibble(period_id = curr_period_id, window = w_label,
                  n_nodes = 0L, n_edges = 0L, n_atypical_edges = 0L,
                  n_countries = 0L, n_components = 0L))
  }
  
  all_inventors <- window_invs %>% distinct(person_id) %>% pull(person_id)
  
  co_pairs_all <- window_invs %>%
    inner_join(window_invs, by = "pat_id", suffix = c(".i", ".j"), relationship = "many-to-many") %>%
    filter(person_id.i < person_id.j) %>%
    distinct(person_id.i, person_id.j)
  
  if (nrow(co_pairs_all) > 0L) {
    g_inv <- graph_from_data_frame(
      co_pairs_all %>% transmute(from = person_id.i, to = person_id.j),
      directed = FALSE
    )
    missing_vertices <- setdiff(all_inventors, V(g_inv)$name)
    if (length(missing_vertices) > 0L) {
      g_inv <- add_vertices(g_inv, nv = length(missing_vertices), name = missing_vertices)
    }
    n_nodes      <- length(all_inventors)
    n_edges      <- ecount(g_inv)
    n_components <- components(g_inv)$no
  } else {
    n_nodes      <- length(all_inventors)
    n_edges      <- 0L
    n_components <- n_nodes
  }
  
  pat_flags_win_atyp <- pat_flags_df %>%
    filter(period_id == curr_period_id, is_atypical_patent) %>%
    select(pat_id)
  
  if (nrow(pat_flags_win_atyp) == 0L) {
    n_atypical_edges <- 0L
  } else {
    window_invs_atyp <- window_invs %>%
      semi_join(pat_flags_win_atyp, by = c("pat_id" = "pat_id"))
    
    if (nrow(window_invs_atyp) == 0L) {
      n_atypical_edges <- 0L
    } else {
      co_pairs_atyp <- window_invs_atyp %>%
        inner_join(window_invs_atyp, by = "pat_id", suffix = c(".i", ".j"),
                   relationship = "many-to-many") %>%
        filter(person_id.i < person_id.j) %>%
        distinct(person_id.i, person_id.j)
      
      n_atypical_edges <- nrow(co_pairs_atyp)
    }
  }
  
  n_countries <- window_invs %>% distinct(ctry_code) %>% nrow()
  
  tibble(
    period_id        = curr_period_id,
    window           = w_label,
    n_nodes          = n_nodes,
    n_edges          = n_edges,
    n_atypical_edges = n_atypical_edges,
    n_countries      = n_countries,
    n_components     = n_components
  )
}

network_summary <- lapply(
  windows$period_id,
  summarize_window_network,
  windows_df   = windows,
  inv_epo_eu   = inv_epo_eu,
  pat_flags_df = pat_atyp_flags
) %>%
  bind_rows() %>%
  arrange(period_id)

print(network_summary)

# 7.3 Share of inventors with no external co-inventor ties in place graph (degree==0)
no_external_co_inventor_ties <- final_panel_df %>%
  group_by(period_id, window) %>%
  summarise(
    n_inventors        = n_distinct(person_id),
    n_single_inventors = sum(degree == 0, na.rm = TRUE),
    share_single       = n_single_inventors / n_inventors,
    .groups = "drop"
  ) %>%
  arrange(period_id)

print(no_external_co_inventor_ties)

# 7.4 Backboning diagnostics (collapse + place size distribution)
place_sizes <- final_panel_df %>%
  filter(!is.na(place_id)) %>%
  count(period_id, window, place_id, name = "n_inventors_place")

backbone_summary <- place_sizes %>%
  group_by(period_id, window) %>%
  summarise(
    n_places              = n_distinct(place_id),
    n_inventors_in_places = sum(n_inventors_place),
    n_places_multi        = sum(n_inventors_place > 1L),
    n_inventors_multi     = sum(if_else(n_inventors_place > 1L, n_inventors_place, 0L)),
    mean_place_size       = mean(n_inventors_place),
    median_place_size     = median(n_inventors_place),
    max_place_size        = max(n_inventors_place),
    .groups = "drop"
  ) %>%
  left_join(
    final_panel_df %>%
      group_by(period_id, window) %>%
      summarise(n_inventors = n_distinct(person_id), .groups = "drop"),
    by = c("period_id", "window")
  ) %>%
  mutate(
    n_isolates         = n_inventors - n_inventors_in_places,
    collapse_ratio     = n_places / n_inventors,
    share_in_multi     = n_inventors_multi / n_inventors,
    share_isolates     = n_isolates / n_inventors,
    share_places_multi = n_places_multi / n_places
  ) %>%
  arrange(period_id)

print(backbone_summary)

# 7.5 Betweenness summary (mapped back to inventors)
backbone_brokerage_summary <- final_panel_df %>%
  group_by(period_id, window) %>%
  summarise(
    n_inventors        = n_distinct(person_id),
    mean_betweenness   = mean(betweenness, na.rm = TRUE),
    median_betweenness = median(betweenness, na.rm = TRUE),
    share_zero_btw     = mean(betweenness == 0, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(period_id)

print(backbone_brokerage_summary)

#=============================================================================--
# 8. Region-window panel + FE regression (NUTS2) ----
#=============================================================================--

# 8.1 Dependent variable: share of atypical patents per NUTS2 & window
pat_region_window <- inv_epo_eu %>%
  left_join(nuts3_to_nuts2, by = "reg_code") %>%
  filter(!is.na(up_reg_code)) %>%
  transmute(pat_id = appln_id, region = up_reg_code) %>%
  distinct() %>%
  inner_join(
    pat_atyp_flags %>% select(pat_id, period_id, window, is_atypical_patent),
    by = "pat_id"
  ) %>%
  distinct()

reg_atyp_panel <- pat_region_window %>%
  group_by(region, period_id, window) %>%
  summarise(
    n_pat_reg  = n_distinct(pat_id),
    n_atyp_reg = sum(is_atypical_patent, na.rm = TRUE),
    share_atyp = if_else(n_pat_reg > 0, n_atyp_reg / n_pat_reg, NA_real_),
    .groups = "drop"
  )

# 8.2 Main X's: (i) share of brokers (window-specific cutoff) and
#               (ii) mean betweenness per NUTS2-region & window

inv_panel_reg <- final_panel_df %>%
  left_join(nuts3_to_nuts2, by = "reg_code") %>%
  filter(!is.na(up_reg_code)) %>%
  rename(region = up_reg_code)

# window-specific broker cutoff: bottom 25% of constraint within each window
inv_panel_reg <- inv_panel_reg %>%
  group_by(period_id, window) %>%
  mutate(
    broker_cutoff_win = quantile(constraint, probs = 0.25, na.rm = TRUE),
    broker_win = if_else(!is.na(constraint) & constraint <= broker_cutoff_win, 1L, 0L)
  ) %>%
  ungroup()

brokers_panel <- inv_panel_reg %>%
  group_by(region, period_id, window) %>%
  summarise(
    n_inv               = n(),
    n_brokers_win       = sum(broker_win),
    share_brokers_win   = if_else(n_inv > 0, n_brokers_win / n_inv, NA_real_),
    mean_deg            = mean(degree, na.rm = TRUE),
    mean_constr         = mean(constraint, na.rm = TRUE),
    mean_betweenness    = mean(betweenness, na.rm = TRUE),      # 2nd proxy
    median_betweenness  = median(betweenness, na.rm = TRUE),
    share_zero_btw      = mean(betweenness == 0, na.rm = TRUE),
    .groups = "drop"
  )


# 8.3 Network controls (DENSITY, ISOLATE, COMMUNITY, INTERREGIONAL)
compute_place_network_controls <- function(curr_period_id, windows_df, inv_epo_eu, nuts3_to_nuts2) {
  
  w <- windows_df %>% filter(period_id == !!curr_period_id)
  w_start <- w$start_y
  w_end   <- w$end_y
  w_label <- w$window
  
  cat(sprintf("Place-graph controls for window %s (%d-%d)\n", w_label, w_start, w_end))
  
  # Window inventors with NUTS2 region
  window_invs <- inv_epo_eu %>%
    filter(app_year >= w_start, app_year <= w_end) %>%
    left_join(nuts3_to_nuts2, by = "reg_code") %>%
    filter(!is.na(up_reg_code)) %>%
    transmute(
      person_id,
      pat_id  = appln_id,
      region  = up_reg_code
    ) %>%
    distinct()
  
  if (nrow(window_invs) == 0L) {
    return(tibble(
      region = character(), period_id = integer(), window = character(),
      DENSITY = numeric(), ISOLATE = numeric(), COMMUNITY = numeric(), INTERREGIONAL = numeric()
    ))
  }
  
  # Dominant region per inventor (mode within window)
  person_region <- window_invs %>%
    count(person_id, region, name = "n_obs") %>%
    arrange(person_id, desc(n_obs)) %>%
    group_by(person_id) %>%
    slice_head(n = 1) %>%
    ungroup() %>%
    select(person_id, region)
  
  all_inventors <- person_region$person_id
  
  # Co-inventor edges
  co_pairs <- window_invs %>%
    inner_join(window_invs, by = "pat_id", suffix = c(".i", ".j"), relationship = "many-to-many") %>%
    filter(person_id.i < person_id.j) %>%
    distinct(person_id.i, person_id.j)
  
  # If no edges: all inventors isolated -> place graph empty
  if (nrow(co_pairs) == 0L) {
    return(
      person_region %>%
        distinct(region) %>%
        transmute(
          region,
          period_id     = curr_period_id,
          window        = w_label,
          DENSITY       = 0,
          ISOLATE       = 1,          # all places are isolates (effectively)
          COMMUNITY     = 0,
          INTERREGIONAL = NA_real_
        )
    )
  }
  
  # Build "place signatures" (structural equivalence)
  neighbors <- bind_rows(
    co_pairs %>% transmute(person_id = person_id.i, neighbor = person_id.j),
    co_pairs %>% transmute(person_id = person_id.j, neighbor = person_id.i)
  ) %>%
    distinct() %>%
    arrange(person_id, neighbor)
  
  place_map <- neighbors %>%
    group_by(person_id) %>%
    summarise(place_sig = paste(sort(unique(neighbor)), collapse = "|"), .groups = "drop")
  
  isolates <- setdiff(all_inventors, place_map$person_id)
  if (length(isolates) > 0L) {
    place_map <- bind_rows(place_map, tibble(person_id = isolates, place_sig = NA_character_))
  }
  
  place_map <- place_map %>%
    mutate(place_id = if_else(is.na(place_sig), NA_integer_, as.integer(factor(place_sig))))
  
  # Non-isolate inventors -> place graph vertices
  place_map_noniso <- place_map %>% filter(!is.na(place_id))
  
  # Map inventor edges to place edges
  co_pairs_places <- co_pairs %>%
    left_join(place_map_noniso %>% select(person_id, place_id),
              by = c("person_id.i" = "person_id")) %>%
    rename(place_i = place_id) %>%
    left_join(place_map_noniso %>% select(person_id, place_id),
              by = c("person_id.j" = "person_id")) %>%
    rename(place_j = place_id) %>%
    filter(!is.na(place_i), !is.na(place_j), place_i != place_j) %>%
    distinct(place_i, place_j)
  
  if (nrow(co_pairs_places) == 0L) {
    # Places exist but no between-place edges
    return(
      person_region %>%
        distinct(region) %>%
        transmute(
          region,
          period_id     = curr_period_id,
          window        = w_label,
          DENSITY       = 0,
          ISOLATE       = 1,
          COMMUNITY     = 0,
          INTERREGIONAL = NA_real_
        )
    )
  }
  
  g_place <- graph_from_data_frame(
    co_pairs_places %>% transmute(from = as.character(place_i), to = as.character(place_j)),
    directed = FALSE
  )
  
  # Ensure all places are vertices
  all_places <- unique(place_map_noniso$place_id)
  missing_places <- setdiff(as.character(all_places), V(g_place)$name)
  if (length(missing_places) > 0L) {
    g_place <- add_vertices(g_place, nv = length(missing_places), name = missing_places)
  }
  
  # Assign each place to a region (modal region among its inventors)
  place_region <- place_map_noniso %>%
    left_join(person_region, by = "person_id") %>%
    filter(!is.na(region)) %>%
    count(place_id, region, name = "n") %>%
    arrange(place_id, desc(n)) %>%
    group_by(place_id) %>%
    slice_head(n = 1) %>%
    ungroup() %>%
    transmute(place_id, region)
  
  reg_list <- sort(unique(place_region$region))
  
  # Edge list with region info for INTERREGIONAL
  edges_df <- as_data_frame(g_place, what = "edges") %>%
    transmute(
      place_from = as.integer(from),
      place_to   = as.integer(to)
    ) %>%
    left_join(place_region %>% rename(region_from = region), by = c("place_from" = "place_id")) %>%
    left_join(place_region %>% rename(region_to   = region), by = c("place_to"   = "place_id")) %>%
    filter(!is.na(region_from), !is.na(region_to))
  
  bind_rows(lapply(reg_list, function(r) {
    
    places_r <- place_region %>% filter(region == r) %>% pull(place_id)
    places_r_chr <- as.character(places_r)
    
    subg <- induced_subgraph(g_place, vids = V(g_place)[name %in% places_r_chr])
    
    v_r <- vcount(subg)
    e_r <- ecount(subg)
    
    dens_r <- if (v_r > 1) 2 * e_r / (v_r * (v_r - 1)) else 0
    
    deg_r <- degree(subg)
    iso_r <- if (length(deg_r) > 0) mean(deg_r == 0) else 1
    
    comm_norm <- if (e_r > 0 && v_r > 0) {
      comm <- cluster_louvain(subg)
      length(unique(membership(comm))) / v_r
    } else {
      0
    }
    
    edges_r <- edges_df %>% filter(region_from == r | region_to == r)
    inter_share <- if (nrow(edges_r) == 0L) NA_real_ else mean(edges_r$region_from != edges_r$region_to)
    
    tibble(
      region        = r,
      period_id     = curr_period_id,
      window        = w_label,
      DENSITY       = dens_r,
      ISOLATE       = iso_r,
      COMMUNITY     = comm_norm,
      INTERREGIONAL = inter_share
    )
  }))
}

network_controls_panel <- lapply(
  windows$period_id,
  compute_place_network_controls,
  windows_df     = windows,
  inv_epo_eu     = inv_epo_eu,
  nuts3_to_nuts2 = nuts3_to_nuts2
) %>%
  bind_rows()


# 8.4 Relatedness + Complexity (safe)
compute_relatedness_complexity <- function(curr_period_id, windows_df, inv_epo_eu, tech_df_eu, nuts3_to_nuts2) {
  
  w <- windows_df %>% filter(period_id == !!curr_period_id)
  w_start <- w$start_y
  w_end   <- w$end_y
  w_label <- w$window
  
  cat(sprintf("RELATEDNESS/COMPLEXITY for window %s (%d-%d)\n", w_label, w_start, w_end))
  
  pats_win <- tech_df_eu %>%
    filter(app_year >= w_start, app_year <= w_end) %>%
    distinct(pat_id, ipc4) %>%
    filter(!is.na(ipc4), ipc4 != "", !grepl("^\\s*$", ipc4))
  
  if (nrow(pats_win) == 0L) return(tibble())
  
  pat_reg <- inv_epo_eu %>%
    filter(app_year >= w_start, app_year <= w_end) %>%
    left_join(nuts3_to_nuts2, by = "reg_code") %>%
    filter(!is.na(up_reg_code)) %>%
    transmute(pat_id = appln_id, region = up_reg_code) %>%
    distinct()
  
  reg_tech <- pat_reg %>%
    inner_join(pats_win, by = "pat_id", relationship = "many-to-many") %>%
    count(region, ipc4, name = "n") %>%
    filter(!is.na(region), !is.na(ipc4), ipc4 != "", !grepl("^\\s*$", ipc4))
  
  if (nrow(reg_tech) == 0L) return(tibble())
  
  mat_df <- reg_tech %>%
    pivot_wider(names_from = ipc4, values_from = n, values_fill = 0)
  
  region_vec <- mat_df$region
  mat_counts <- as.matrix(mat_df[, -1, drop = FALSE])
  rownames(mat_counts) <- region_vec
  
  mat_bin <- (mat_counts > 0) * 1
  if (ncol(mat_bin) < 2) {
    return(tibble(region = region_vec, RELATEDNESS = NA_real_, COMPLEXITY = NA_real_,
                  period_id = curr_period_id, window = w_label))
  }
  
  denom <- sqrt(colSums(mat_bin^2))
  denom[denom == 0] <- NA_real_
  mat_norm <- sweep(mat_bin, 2, denom, "/")
  
  rel_tt <- t(mat_norm) %*% mat_norm
  rel_tt <- as.matrix(rel_tt)
  diag(rel_tt) <- NA_real_
  
  tech_names <- colnames(mat_bin)
  
  relatedness_vec <- sapply(seq_len(nrow(mat_bin)), function(i) {
    techs <- which(mat_bin[i, ] > 0)
    if (length(techs) <= 1) return(NA_real_)
    sub_mat <- rel_tt[techs, techs, drop = FALSE]
    mean(sub_mat[upper.tri(sub_mat)], na.rm = TRUE)
  })
  
  relatedness_df <- tibble(
    region      = region_vec,
    RELATEDNESS = as.numeric(relatedness_vec),
    period_id   = curr_period_id,
    window      = w_label
  )
  
  comp_tech <- tryCatch(EconGeo::mort(mat_bin, rca = FALSE),
                        error = function(e) rep(NA_real_, ncol(mat_bin)))
  
  comp_tech_df <- tibble(ipc4 = tech_names, COMPLEXITY_tech = as.numeric(comp_tech))
  
  reg_comp <- reg_tech %>%
    left_join(comp_tech_df, by = "ipc4") %>%
    group_by(region) %>%
    summarise(COMPLEXITY = median(COMPLEXITY_tech, na.rm = TRUE), .groups = "drop") %>%
    mutate(period_id = curr_period_id, window = w_label)
  
  relatedness_df %>%
    full_join(reg_comp, by = c("region", "period_id", "window"))
}

relatedness_complexity_panel <- lapply(
  windows$period_id,
  compute_relatedness_complexity,
  windows_df     = windows,
  inv_epo_eu     = inv_epo_eu,
  tech_df_eu     = tech_df_eu,
  nuts3_to_nuts2 = nuts3_to_nuts2
) %>% bind_rows()

# 8.5 Merge region-window panel
panel_reg <- reg_atyp_panel %>%
  left_join(brokers_panel, by = c("region", "period_id", "window")) %>%
  left_join(network_controls_panel, by = c("region", "period_id", "window")) %>%
  left_join(relatedness_complexity_panel, by = c("region", "period_id", "window"))


cat("\nPanel summary (raw):\n")
print(summary(panel_reg))

#=============================================================================--
# 8.6 Region-window panel + lags (clean) ----
#=============================================================================--

# --- Outcome: share of atypical patents in region-window (already built as reg_atyp_panel) ---
# reg_atyp_panel: region, period_id, window, n_pat_reg, n_atyp_reg, share_atyp

# --- Build region-window brokerage aggregates (both dummy-share and continuous) ---
# inv_panel_reg is final_panel_df + NUTS2 mapping (region) with inventor-level degree/constraint/betweenness

# (A) window-specific broker dummy at inventor level -> region share
inv_panel_reg2 <- inv_panel_reg %>%
  group_by(period_id, window) %>%
  mutate(
    broker_cutoff_win = quantile(constraint, probs = 0.25, na.rm = TRUE),
    broker_win = if_else(!is.na(constraint) & constraint <= broker_cutoff_win, 1L, 0L)
  ) %>%
  ungroup()

brokers_panel <- inv_panel_reg2 %>%
  group_by(region, period_id, window) %>%
  summarise(
    n_inv             = n(),
    share_brokers_win = mean(broker_win, na.rm = TRUE),
    mean_constr       = mean(constraint, na.rm = TRUE),
    mean_brokerage    = mean(1 - constraint, na.rm = TRUE),  # Burt-style brokerage (higher = more brokerage)
    mean_betweenness  = mean(betweenness, na.rm = TRUE),
    .groups = "drop"
  )

# --- Merge everything into one region-window panel ---
panel_reg <- reg_atyp_panel %>%
  left_join(brokers_panel, by = c("region", "period_id", "window")) %>%
  left_join(network_controls_panel, by = c("region", "period_id", "window")) %>%
  left_join(relatedness_complexity_panel, by = c("region", "period_id", "window"))

# --- Optional: cumulative DV (stock), for robustness ---
panel_reg <- panel_reg %>%
  arrange(region, period_id) %>%
  group_by(region) %>%
  mutate(
    cum_pat   = cumsum(replace_na(n_pat_reg, 0)),
    cum_atyp  = cumsum(replace_na(n_atyp_reg, 0)),
    share_atyp_cum = if_else(cum_pat > 0, cum_atyp / cum_pat, NA_real_)
  ) %>%
  ungroup()

# --- Create 1-window lags (your causal timing choice) ---
panel_reg_lagged <- panel_reg %>%
  arrange(region, period_id) %>%
  group_by(region) %>%
  mutate(
    L_share_brokers_win = lag(share_brokers_win, 1),
    L_mean_betweenness  = lag(mean_betweenness, 1),
    L_mean_constr       = lag(mean_constr, 1),
    L_mean_brokerage    = lag(mean_brokerage, 1),
    
    L_ISOLATE       = lag(ISOLATE, 1),
    L_COMMUNITY     = lag(COMMUNITY, 1),
    L_DENSITY       = lag(DENSITY, 1),
    L_RELATEDNESS   = lag(RELATEDNESS, 1),
    L_COMPLEXITY    = lag(COMPLEXITY, 1),
    L_INTERREGIONAL = lag(INTERREGIONAL, 1)
  ) %>%
  ungroup() %>%
  filter(!is.na(share_atyp))   # keep Y defined; model-specific filtering happens below


cat("\nPanel summary (after lagging, before model-specific filtering):\n")
print(summary(panel_reg_lagged))
cat("\nRows after lagging (Y only):", nrow(panel_reg_lagged), "\n")

#=============================================================================--
# 8.7 Helper: TWFE + cluster-robust SE by region (safe) ----
#=============================================================================--

run_twfe <- function(df, fml, model_name) {
  pdata <- plm::pdata.frame(df, index = c("region", "period_id"))
  fit   <- plm::plm(fml, data = pdata, model = "within", effect = "twoways")
  vc    <- sandwich::vcovHC(fit, type = "HC1", cluster = "group")
  
  cat("\n============================================================\n")
  cat(model_name, "\n")
  cat("N:", nrow(df),
      "| Regions:", dplyr::n_distinct(df$region),
      "| Periods:", dplyr::n_distinct(df$period_id), "\n")
  cat("============================================================\n")
  print(summary(fit))
  cat("\n--- Robust SE (HC1, clustered by region) ---\n")
  print(lmtest::coeftest(fit, vcov. = vc))
  
  invisible(list(model = fit, vcov = vc))
}

controls_str <- paste(
  "L_ISOLATE", "L_COMMUNITY", "L_DENSITY",
  "L_RELATEDNESS", "L_COMPLEXITY", "L_INTERREGIONAL",
  sep = " + "
)

#=============================================================================--
# 8.8 Main models (flow DV) + stable interaction ----
#=============================================================================--

# Model 1: broker share (dummy based)
df_m1 <- panel_reg_lagged %>%
  filter(
    !is.na(L_share_brokers_win),
    !is.na(L_ISOLATE), !is.na(L_COMMUNITY), !is.na(L_DENSITY),
    !is.na(L_RELATEDNESS), !is.na(L_COMPLEXITY), !is.na(L_INTERREGIONAL)
  )
fml_1 <- as.formula(paste("share_atyp ~ L_share_brokers_win +", controls_str))
res_m1 <- run_twfe(df_m1, fml_1, "TWFE Model 1 (flow DV): L_share_brokers_win")

# Model 2: mean betweenness (raw)
df_m2 <- panel_reg_lagged %>%
  filter(
    !is.na(L_mean_betweenness),
    !is.na(L_ISOLATE), !is.na(L_COMMUNITY), !is.na(L_DENSITY),
    !is.na(L_RELATEDNESS), !is.na(L_COMPLEXITY), !is.na(L_INTERREGIONAL)
  )
fml_2 <- as.formula(paste("share_atyp ~ L_mean_betweenness +", controls_str))
res_m2 <- run_twfe(df_m2, fml_2, "TWFE Model 2 (flow DV): L_mean_betweenness")

# Model 3: continuous Burt-style brokerage (1 - constraint)
df_m3 <- panel_reg_lagged %>%
  filter(
    !is.na(L_mean_brokerage),
    !is.na(L_ISOLATE), !is.na(L_COMMUNITY), !is.na(L_DENSITY),
    !is.na(L_RELATEDNESS), !is.na(L_COMPLEXITY), !is.na(L_INTERREGIONAL)
  )
fml_3 <- as.formula(paste("share_atyp ~ L_mean_brokerage +", controls_str))
res_m3 <- run_twfe(df_m3, fml_3, "TWFE Model 3 (flow DV): L_mean_brokerage = 1 - constraint")

# Model 4: interaction (numerically stable)
#   betweenness is tiny & zero-inflated -> log1p then z-score both terms before interaction
df_int <- panel_reg_lagged %>%
  mutate(
    L_btw_log = log1p(L_mean_betweenness),
    L_btw_z   = as.numeric(scale(L_btw_log)),
    L_rel_z   = as.numeric(scale(L_RELATEDNESS))
  ) %>%
  filter(
    !is.na(L_btw_z), !is.na(L_rel_z),
    !is.na(L_ISOLATE), !is.na(L_COMMUNITY), !is.na(L_DENSITY),
    !is.na(L_COMPLEXITY), !is.na(L_INTERREGIONAL)
  )

controls_int_str <- paste(
  "L_ISOLATE", "L_COMMUNITY", "L_DENSITY",
  "L_COMPLEXITY", "L_INTERREGIONAL",
  sep = " + "
)
fml_int <- as.formula(paste("share_atyp ~ L_btw_z * L_rel_z +", controls_int_str))
res_int <- run_twfe(df_int, fml_int, "TWFE Model 4 (flow DV): z(log1p(btw)) * z(relatedness)")

#=============================================================================--
# 8.9 Optional robustness: cumulative DV (stock) ----
#=============================================================================--

df_cum <- panel_reg_lagged %>%
  filter(
    !is.na(share_atyp_cum),
    !is.na(L_mean_brokerage),
    !is.na(L_ISOLATE), !is.na(L_COMMUNITY), !is.na(L_DENSITY),
    !is.na(L_RELATEDNESS), !is.na(L_COMPLEXITY), !is.na(L_INTERREGIONAL)
  )
fml_cum <- as.formula(paste("share_atyp_cum ~ L_mean_brokerage +", controls_str))
res_cum <- run_twfe(df_cum, fml_cum, "TWFE Robustness (cumulative DV): share_atyp_cum ~ L_mean_brokerage")

#=============================================================================--
# 9. Save workspace (optional) ----
#=============================================================================--

# save.image("atypical_backboned_9windows_fullpipeline.RData")


#=============================================================================--
# 10. panel reg2 ----
#=============================================================================--
# ============================================================================-
# Region-window FE regressions (clean, consistent, numerically stable)
# Assumes you already have:
#   panel_reg  (region-window panel with share_atyp, share_brokers_win,
#               mean_betweenness, ISOLATE, COMMUNITY, DENSITY,
#               RELATEDNESS, COMPLEXITY, INTERREGIONAL)
# and packages: dplyr, plm, sandwich, lmtest
# ==============================================================================-

panel_reg_lagged <- panel_reg %>%
  arrange(region, period_id) %>%
  group_by(region) %>%
  mutate(
    L_share_brokers_win = dplyr::lag(share_brokers_win, 1),
    L_mean_betweenness  = dplyr::lag(mean_betweenness, 1),
    L_ISOLATE           = dplyr::lag(ISOLATE, 1),
    L_COMMUNITY         = dplyr::lag(COMMUNITY, 1),
    L_DENSITY           = dplyr::lag(DENSITY, 1),
    L_RELATEDNESS       = dplyr::lag(RELATEDNESS, 1),
    L_COMPLEXITY        = dplyr::lag(COMPLEXITY, 1),
    L_INTERREGIONAL     = dplyr::lag(INTERREGIONAL, 1)
  ) %>%
  ungroup() %>%
  filter(!is.na(share_atyp))

cat("\nPanel summary (after lagging, before model-specific filtering):\n")
print(summary(panel_reg_lagged))
cat("\nRows after lagging (Y only):", nrow(panel_reg_lagged), "\n")


# 2) Helper: run TWFE + cluster-robust SE by region
run_twfe <- function(df, fml, model_name) {
  pdata <- plm::pdata.frame(df, index = c("region", "period_id"))
  fit   <- plm::plm(fml, data = pdata, model = "within", effect = "twoways")
  vc    <- sandwich::vcovHC(fit, type = "HC1", cluster = "group")
  
  cat("\n============================================================\n")
  cat(model_name, "\n")
  cat("N:", nrow(df),
      "| Regions:", dplyr::n_distinct(df$region),
      "| Periods:", dplyr::n_distinct(df$period_id), "\n")
  cat("============================================================\n")
  print(summary(fit))
  cat("\n--- Robust SE (HC1, clustered by region) ---\n")
  print(lmtest::coeftest(fit, vcov. = vc))
  
  invisible(list(model = fit, vcov = vc))
}

# Controls as a STRING (safe to paste into formula)
controls_str <- paste(
  "L_ISOLATE", "L_COMMUNITY", "L_DENSITY",
  "L_RELATEDNESS", "L_COMPLEXITY", "L_INTERREGIONAL",
  sep = " + "
)

# ---------------------------------------------------------------------- -
# Model 1: L_share_brokers_win
# ---------------------------------------------------------------------- -
df_m1 <- panel_reg_lagged %>%
  filter(
    !is.na(L_share_brokers_win),
    !is.na(L_ISOLATE), !is.na(L_COMMUNITY), !is.na(L_DENSITY),
    !is.na(L_RELATEDNESS), !is.na(L_COMPLEXITY), !is.na(L_INTERREGIONAL)
  )

fml_1 <- as.formula(paste("share_atyp ~ L_share_brokers_win +", controls_str))
res_m1 <- run_twfe(df_m1, fml_1, "TWFE Model 1: L_share_brokers_win")

# ---------------------------------------------------------------------- -
# Model 2: L_mean_betweenness
# ---------------------------------------------------------------------- -
df_m2 <- panel_reg_lagged %>%
  filter(
    !is.na(L_mean_betweenness),
    !is.na(L_ISOLATE), !is.na(L_COMMUNITY), !is.na(L_DENSITY),
    !is.na(L_RELATEDNESS), !is.na(L_COMPLEXITY), !is.na(L_INTERREGIONAL)
  )

fml_2 <- as.formula(paste("share_atyp ~ L_mean_betweenness +", controls_str))
res_m2 <- run_twfe(df_m2, fml_2, "TWFE Model 2: L_mean_betweenness")

# ---------------------------------------------------------------------- - 
# Model 3: both proxies
# ---------------------------------------------------------------------- -
df_m3 <- panel_reg_lagged %>%
  filter(
    !is.na(L_share_brokers_win), !is.na(L_mean_betweenness),
    !is.na(L_ISOLATE), !is.na(L_COMMUNITY), !is.na(L_DENSITY),
    !is.na(L_RELATEDNESS), !is.na(L_COMPLEXITY), !is.na(L_INTERREGIONAL)
  )

fml_3 <- as.formula(paste(
  "share_atyp ~ L_share_brokers_win + L_mean_betweenness +", controls_str
))
res_m3 <- run_twfe(df_m3, fml_3, "TWFE Model 3: both proxies")

# ---------------------------------------------------------------------- - 
# Model 4: interaction (numerically stable): z(log1p(betweenness)) * z(relatedness)
#   Reason: your betweenness is extremely tiny and often zero -> interaction ill-conditioned.
# ---------------------------------------------------------------------- -
df_int <- panel_reg_lagged %>%
  mutate(
    L_btw_log = log1p(L_mean_betweenness),
    L_btw_z   = as.numeric(scale(L_btw_log)),
    L_rel_z   = as.numeric(scale(L_RELATEDNESS))
  ) %>%
  filter(
    !is.na(L_btw_z), !is.na(L_rel_z),
    !is.na(L_ISOLATE), !is.na(L_COMMUNITY), !is.na(L_DENSITY),
    !is.na(L_COMPLEXITY), !is.na(L_INTERREGIONAL)
  )

controls_int_str <- paste(
  "L_ISOLATE", "L_COMMUNITY", "L_DENSITY",
  "L_COMPLEXITY", "L_INTERREGIONAL",
  sep = " + "
)

fml_int <- as.formula(paste(
  "share_atyp ~ L_btw_z * L_rel_z +", controls_int_str
))
res_int <- run_twfe(df_int, fml_int, "TWFE Model 4: interaction L_btw_z * L_rel_z")



#=============================================================================--
# +1 Cumulatie atypical patent ----
#=============================================================================--

#=============================================================================--
# 8.6–8.9 Regressions: FLOW + CUMULATIVE atypical shares (TWFE) 
#=============================================================================--

# Assumes you already built:
#   panel_reg with columns:
#     region, period_id, window,
#     n_pat_reg, n_atyp_reg, share_atyp,
#     share_brokers_win, mean_betweenness, mean_constr, mean_brokerage,
#     DENSITY, ISOLATE, COMMUNITY, INTERREGIONAL, RELATEDNESS, COMPLEXITY
#
# and packages loaded: dplyr, plm, sandwich, lmtest

# ------------------------------------------------------------------------------ -
# 1) Build cumulative DV ("stock"): cumulative atypical share up to current window
# ------------------------------------------------------------------------------ -
panel_reg <- panel_reg %>%
  arrange(region, period_id) %>%
  group_by(region) %>%
  mutate(
    cum_pat       = cumsum(replace_na(n_pat_reg, 0)),
    cum_atyp      = cumsum(replace_na(n_atyp_reg, 0)),
    share_atyp_cum = if_else(cum_pat > 0, cum_atyp / cum_pat, NA_real_)
  ) %>%
  ungroup()

# ------------------------------------------------------------------------------ -
# 2) Create 1-window lags for all RHS variables (same for flow + cum models)
# ------------------------------------------------------------------------------ -
panel_reg_lagged <- panel_reg %>%
  arrange(region, period_id) %>%
  group_by(region) %>%
  mutate(
    L_share_brokers_win = dplyr::lag(share_brokers_win, 1),
    L_mean_betweenness  = dplyr::lag(mean_betweenness, 1),
    L_mean_constr       = dplyr::lag(mean_constr, 1),
    L_mean_brokerage    = dplyr::lag(mean_brokerage, 1),
    
    L_ISOLATE           = dplyr::lag(ISOLATE, 1),
    L_COMMUNITY         = dplyr::lag(COMMUNITY, 1),
    L_DENSITY           = dplyr::lag(DENSITY, 1),
    L_RELATEDNESS       = dplyr::lag(RELATEDNESS, 1),
    L_COMPLEXITY        = dplyr::lag(COMPLEXITY, 1),
    L_INTERREGIONAL     = dplyr::lag(INTERREGIONAL, 1)
  ) %>%
  ungroup()

cat("\nPanel summary (after lagging; contains both flow and cumulative DV):\n")
print(summary(panel_reg_lagged))
cat("\nRows total:", nrow(panel_reg_lagged), "\n")

# ------------------------------------------------------------------------------ -
# 3) Helper: TWFE + cluster-robust SE by region
# ------------------------------------------------------------------------------ -
run_twfe <- function(df, fml, model_name) {
  pdata <- plm::pdata.frame(df, index = c("region", "period_id"))
  fit   <- plm::plm(fml, data = pdata, model = "within", effect = "twoways")
  vc    <- sandwich::vcovHC(fit, type = "HC1", cluster = "group")
  
  cat("\n============================================================\n")
  cat(model_name, "\n")
  cat("N:", nrow(df),
      "| Regions:", dplyr::n_distinct(df$region),
      "| Periods:", dplyr::n_distinct(df$period_id), "\n")
  cat("============================================================\n")
  print(summary(fit))
  cat("\n--- Robust SE (HC1, clustered by region) ---\n")
  print(lmtest::coeftest(fit, vcov. = vc))
  
  invisible(list(model = fit, vcov = vc))
}

# Controls (string-safe for formula paste)
controls_str <- paste(
  "L_ISOLATE", "L_COMMUNITY", "L_DENSITY",
  "L_RELATEDNESS", "L_COMPLEXITY", "L_INTERREGIONAL",
  sep = " + "
)

# ------------------------------------------------------------------------------ -
# 4) Model runner: runs SAME model for flow DV and cumulative DV
# ------------------------------------------------------------------------------ -
run_flow_and_cum <- function(df_base, rhs_str, tag) {
  
  # ---------------- FLOW ---------------- -
  df_flow <- df_base %>%
    filter(!is.na(share_atyp))
  
  fml_flow <- as.formula(paste("share_atyp ~", rhs_str))
  run_twfe(df_flow, fml_flow, paste0(tag, " (FLOW DV: share_atyp)"))
  
  # -------------- CUMULATIVE ------------ -
  df_cum <- df_base %>%
    filter(!is.na(share_atyp_cum))
  
  fml_cum <- as.formula(paste("share_atyp_cum ~", rhs_str))
  run_twfe(df_cum, fml_cum, paste0(tag, " (CUMULATIVE DV: share_atyp_cum)"))
}

# ------------------------------------------------------------------------------ -
# 5) Base filtering: keep rows where controls exist (then model-specific filters)
# ------------------------------------------------------------------------------ -
df_base <- panel_reg_lagged %>%
  filter(
    !is.na(L_ISOLATE), !is.na(L_COMMUNITY), !is.na(L_DENSITY),
    !is.na(L_RELATEDNESS), !is.na(L_COMPLEXITY), !is.na(L_INTERREGIONAL)
  )

#=============================================================================--
# 8.8 Main models: FLOW + CUMULATIVE versions
#=============================================================================--

# Model 1: broker share (dummy-based upstream)
df_m1 <- df_base %>% filter(!is.na(L_share_brokers_win))
rhs_1 <- paste("L_share_brokers_win +", controls_str)
run_flow_and_cum(df_m1, rhs_1, "TWFE Model 1: L_share_brokers_win")

# Model 2: mean betweenness
df_m2 <- df_base %>% filter(!is.na(L_mean_betweenness))
rhs_2 <- paste("L_mean_betweenness +", controls_str)
run_flow_and_cum(df_m2, rhs_2, "TWFE Model 2: L_mean_betweenness")

# Model 3: continuous Burt brokerage (1 - constraint)
df_m3 <- df_base %>% filter(!is.na(L_mean_brokerage))
rhs_3 <- paste("L_mean_brokerage +", controls_str)
run_flow_and_cum(df_m3, rhs_3, "TWFE Model 3: L_mean_brokerage = 1 - constraint")

# Model 3b (optional): mean constraint instead (interpretation flips)
df_m3b <- df_base %>% filter(!is.na(L_mean_constr))
rhs_3b <- paste("L_mean_constr +", controls_str)
run_flow_and_cum(df_m3b, rhs_3b, "TWFE Model 3b: L_mean_constr (Burt constraint)")

# Model 4: stable interaction (z(log1p(betweenness)) * z(relatedness))
df_int <- df_base %>%
  mutate(
    L_btw_log = log1p(L_mean_betweenness),
    L_btw_z   = as.numeric(scale(L_btw_log)),
    L_rel_z   = as.numeric(scale(L_RELATEDNESS))
  ) %>%
  filter(!is.na(L_btw_z), !is.na(L_rel_z))

controls_int_str <- paste(
  "L_ISOLATE", "L_COMMUNITY", "L_DENSITY",
  "L_COMPLEXITY", "L_INTERREGIONAL",
  sep = " + "
)

rhs_4 <- paste("L_btw_z * L_rel_z +", controls_int_str)
run_flow_and_cum(df_int, rhs_4, "TWFE Model 4: interaction z(log1p(btw)) * z(relatedness)")

#=============================================================================--
# End
#=============================================================================--


#=============================================================================--
# ROBUSTNESS: vary brokerage threshold + atypical cutoff, rerun FLOW + CUM TWFE ----
#=============================================================================--

# --- If ipc_pair_stats is not in memory, load it ---
# ipc_pair_stats <- data.table::fread("ipc_pair_zscores_eu.csv") %>% as_tibble()

# --- USER-SET GRIDS ---
broker_ps <- c(0.10, 0.25, 0.40)          # brokerage threshold: bottom p of constraint
z_cuts    <- c(0, -0.5, -1.0, -1.5)       # atypical cutoff: Z_score < z_cut

# --- Helper: build patent atypical flags for a given z_cut ---
make_pat_flags_for_zcut <- function(z_cut, ipc_pair_stats, tech_df_eu, inv_epo_eu, windows) {
  
  # mark atypical pairs under this cutoff
  pair_atyp <- ipc_pair_stats %>%
    mutate(is_atypical = (Z_score < z_cut)) %>%
    select(ipc_i, ipc_j, period_id, window, is_atypical)
  
  # build patent-level atypicality within each period using same target logic
  pat_flags_list <- lapply(windows$period_id, function(pid) {
    
    w <- windows %>% filter(period_id == pid)
    start_y <- w$start_y; end_y <- w$end_y
    
    eu_pat_ids <- inv_epo_eu %>%
      filter(app_year >= start_y, app_year <= end_y) %>%
      distinct(appln_id) %>%
      pull(appln_id)
    
    if (length(eu_pat_ids) == 0L) return(tibble())
    
    target_codes <- tech_df_eu %>%
      filter(pat_id %in% eu_pat_ids) %>%
      distinct(pat_id, ipc4)
    
    pat_sizes <- target_codes %>% count(pat_id, name = "n_codes")
    target_codes <- target_codes %>% inner_join(pat_sizes %>% filter(n_codes >= 2L), by = "pat_id")
    if (nrow(target_codes) == 0L) return(tibble())
    
    target_pairs <- target_codes %>%
      inner_join(target_codes, by = "pat_id", suffix = c(".i", ".j")) %>%
      filter(ipc4.i < ipc4.j) %>%
      transmute(pat_id, ipc_i = ipc4.i, ipc_j = ipc4.j)
    
    pat_flags <- target_pairs %>%
      left_join(pair_atyp %>% filter(period_id == pid) %>% select(ipc_i, ipc_j, is_atypical),
                by = c("ipc_i", "ipc_j")) %>%
      mutate(is_atypical = replace_na(is_atypical, FALSE)) %>%
      group_by(pat_id) %>%
      summarise(
        is_atypical_patent = any(is_atypical),
        .groups = "drop"
      ) %>%
      mutate(period_id = pid, window = w$window)
    
    pat_flags
  })
  
  bind_rows(pat_flags_list)
}

# --- Helper: run FLOW + CUM models given a panel_reg built for this scenario ---
run_flow_cum_suite <- function(panel_reg, scenario_tag) {
  
  # cumulative DV
  panel_reg <- panel_reg %>%
    arrange(region, period_id) %>%
    group_by(region) %>%
    mutate(
      cum_pat        = cumsum(replace_na(n_pat_reg, 0)),
      cum_atyp       = cumsum(replace_na(n_atyp_reg, 0)),
      share_atyp_cum = if_else(cum_pat > 0, cum_atyp / cum_pat, NA_real_)
    ) %>%
    ungroup()
  
  # lags
  panel_reg_lagged <- panel_reg %>%
    arrange(region, period_id) %>%
    group_by(region) %>%
    mutate(
      L_share_brokers_win = dplyr::lag(share_brokers_win, 1),
      L_mean_betweenness  = dplyr::lag(mean_betweenness, 1),
      L_mean_brokerage    = dplyr::lag(mean_brokerage, 1),
      
      L_ISOLATE       = dplyr::lag(ISOLATE, 1),
      L_COMMUNITY     = dplyr::lag(COMMUNITY, 1),
      L_DENSITY       = dplyr::lag(DENSITY, 1),
      L_RELATEDNESS   = dplyr::lag(RELATEDNESS, 1),
      L_COMPLEXITY    = dplyr::lag(COMPLEXITY, 1),
      L_INTERREGIONAL = dplyr::lag(INTERREGIONAL, 1)
    ) %>%
    ungroup()
  
  run_twfe <- function(df, fml, model_name) {
    pdata <- plm::pdata.frame(df, index = c("region", "period_id"))
    fit   <- plm::plm(fml, data = pdata, model = "within", effect = "twoways")
    vc    <- sandwich::vcovHC(fit, type = "HC1", cluster = "group")
    out   <- lmtest::coeftest(fit, vcov. = vc)
    list(fit = fit, ct = out)
  }
  
  controls_str <- paste(
    "L_ISOLATE", "L_COMMUNITY", "L_DENSITY",
    "L_RELATEDNESS", "L_COMPLEXITY", "L_INTERREGIONAL",
    sep = " + "
  )
  
  df_base <- panel_reg_lagged %>%
    filter(
      !is.na(L_ISOLATE), !is.na(L_COMMUNITY), !is.na(L_DENSITY),
      !is.na(L_RELATEDNESS), !is.na(L_COMPLEXITY), !is.na(L_INTERREGIONAL)
    )
  
  # interaction stable
  df_int <- df_base %>%
    mutate(
      L_btw_log = log1p(L_mean_betweenness),
      L_btw_z   = as.numeric(scale(L_btw_log)),
      L_rel_z   = as.numeric(scale(L_RELATEDNESS))
    ) %>%
    filter(!is.na(L_btw_z), !is.na(L_rel_z))
  
  controls_int_str <- paste(
    "L_ISOLATE", "L_COMMUNITY", "L_DENSITY",
    "L_COMPLEXITY", "L_INTERREGIONAL",
    sep = " + "
  )
  
  models <- list(
    M1 = paste("L_share_brokers_win +", controls_str),
    M2 = paste("L_mean_betweenness +", controls_str),
    M3 = paste("L_mean_brokerage +", controls_str),
    M4 = paste("L_btw_z * L_rel_z +", controls_int_str)
  )
  
  # results holder
  res_tbl <- list()
  
  # FLOW + CUM for each model
  for (m in names(models)) {
    
    if (m == "M4") {
      df_use <- df_int
    } else {
      df_use <- df_base
    }
    
    rhs <- models[[m]]
    
    # FLOW
    df_flow <- df_use %>% filter(!is.na(share_atyp))
    fml_f   <- as.formula(paste("share_atyp ~", rhs))
    r_f     <- run_twfe(df_flow, fml_f, paste0(scenario_tag, " | ", m, " | FLOW"))
    
    # CUM
    df_cum  <- df_use %>% filter(!is.na(share_atyp_cum))
    fml_c   <- as.formula(paste("share_atyp_cum ~", rhs))
    r_c     <- run_twfe(df_cum, fml_c, paste0(scenario_tag, " | ", m, " | CUM"))
    
    # extract key coefficient(s)
    get_coef <- function(ct, term) {
      if (!term %in% rownames(ct)) return(tibble(term = term, est = NA_real_, se = NA_real_, p = NA_real_))
      tibble(term = term, est = ct[term, 1], se = ct[term, 2], p = ct[term, 4])
    }
    
    key_terms <- switch(
      m,
      M1 = c("L_share_brokers_win"),
      M2 = c("L_mean_betweenness"),
      M3 = c("L_mean_brokerage"),
      M4 = c("L_btw_z", "L_rel_z", "L_btw_z:L_rel_z")
    )
    
    res_tbl[[paste0(m, "_FLOW")]] <- bind_rows(lapply(key_terms, \(t) get_coef(r_f$ct, t))) %>%
      mutate(model = m, dv = "flow", scenario = scenario_tag)
    
    res_tbl[[paste0(m, "_CUM")]] <- bind_rows(lapply(key_terms, \(t) get_coef(r_c$ct, t))) %>%
      mutate(model = m, dv = "cum", scenario = scenario_tag)
  }
  
  bind_rows(res_tbl)
}

#=============================================================================--
# MAIN ROBUSTNESS LOOP
#=============================================================================--

all_results <- list()

for (p in broker_ps) {
  for (zc in z_cuts) {
    
    cat("\n\n================ ROBUSTNESS RUN =================\n")
    cat("Broker p =", p, " | Z cutoff =", zc, "\n")
    cat("================================================\n")
    
    # 1) recompute patent atypical flags for this Z cutoff
    pat_flags_z <- make_pat_flags_for_zcut(
      z_cut = zc,
      ipc_pair_stats = ipc_pair_stats,
      tech_df_eu = tech_df_eu,
      inv_epo_eu = inv_epo_eu,
      windows = windows
    )
    
    # 2) rebuild region-window DV under this cutoff
    pat_region_window_z <- inv_epo_eu %>%
      left_join(nuts3_to_nuts2, by = "reg_code") %>%
      filter(!is.na(up_reg_code)) %>%
      transmute(pat_id = appln_id, region = up_reg_code) %>%
      distinct() %>%
      inner_join(
        pat_flags_z %>% select(pat_id, period_id, window, is_atypical_patent),
        by = "pat_id"
      )
    
    reg_atyp_panel_z <- pat_region_window_z %>%
      group_by(region, period_id, window) %>%
      summarise(
        n_pat_reg  = n_distinct(pat_id),
        n_atyp_reg = sum(is_atypical_patent, na.rm = TRUE),
        share_atyp = if_else(n_pat_reg > 0, n_atyp_reg / n_pat_reg, NA_real_),
        .groups = "drop"
      )
    
    # 3) rebuild brokerage variables under this broker threshold p
    inv_panel_reg2 <- inv_panel_reg %>%
      group_by(period_id, window) %>%
      mutate(
        broker_cutoff_win = quantile(constraint, probs = p, na.rm = TRUE),
        broker_win = if_else(!is.na(constraint) & constraint <= broker_cutoff_win, 1L, 0L)
      ) %>%
      ungroup()
    
    brokers_panel_zp <- inv_panel_reg2 %>%
      group_by(region, period_id, window) %>%
      summarise(
        n_inv             = n(),
        share_brokers_win = mean(broker_win, na.rm = TRUE),
        mean_brokerage    = mean(1 - constraint, na.rm = TRUE),
        mean_betweenness  = mean(betweenness, na.rm = TRUE),
        .groups = "drop"
      )
    
    # 4) scenario panel
    panel_s <- reg_atyp_panel_z %>%
      left_join(brokers_panel_zp, by = c("region", "period_id", "window")) %>%
      left_join(network_controls_panel, by = c("region", "period_id", "window")) %>%
      left_join(relatedness_complexity_panel, by = c("region", "period_id", "window"))
    
    # 5) run suite (flow + cum)
    tag <- paste0("p=", p, " | z<", zc)
    all_results[[tag]] <- run_flow_cum_suite(panel_s, tag)
  }
}

robustness_results <- bind_rows(all_results)

cat("\nRobustness results (key terms only) preview:\n")
print(head(robustness_results, 30))

data.table::fwrite(robustness_results, "robustness_broker_threshold_x_atyp_cutoff.csv")
cat("\nSaved: robustness_broker_threshold_x_atyp_cutoff.csv\n")


#=============================================================================--
# Visualizing the diff models
#=============================================================================--

# ============================================================
# Paper-style theme + caption (apply to your existing plots)
# ============================================================

theme_req <- c("ggplot2", "scales")
to_install <- setdiff(theme_req, rownames(installed.packages()))
if (length(to_install)) install.packages(to_install)
suppressPackageStartupMessages(invisible(lapply(theme_req, library, character.only = TRUE)))

caption_text <- "Source: own calculation on OECD REGPAT data"

theme_paper_serif <- theme_minimal(base_family = "serif") +
  theme(
    plot.background  = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    panel.grid.major = element_line(colour = "grey90", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    
    axis.title = element_text(size = 12),
    axis.text  = element_text(size = 10),
    
    plot.title = element_text(
      size = 16, face = "bold", hjust = 0,
      margin = margin(t = 6, r = 0, b = 2, l = 6)
    ),
    plot.subtitle = element_text(
      size = 11, colour = "grey30", hjust = 0,
      margin = margin(t = 0, r = 0, b = 6, l = 6)
    ),
    plot.caption = element_text(
      size = 8, colour = "grey40", hjust = 1,
      margin = margin(t = 6, r = 6, b = 4, l = 0)
    ),
    
    legend.position = "right",
    legend.title    = element_text(size = 9, face = "bold"),
    legend.text     = element_text(size = 8)
  )

# ============================================================
# 1) Apply to your EXISTING line plot object (keep your code)
#    Example assumes your object is called: gg_line
# ============================================================

gg_line <- gg_line +
  labs(caption = caption_text) +
  theme_paper_serif +
  # nicer default styling for line charts:
  theme(
    panel.grid.major.x = element_blank(),
    axis.line = element_line(colour = "grey70", linewidth = 0.25)
  ) +
  # if your line plot maps colour to scenarios, this improves readability:
  scale_colour_viridis_d(option = "viridis", direction = -1, name = "Scenario") +
  scale_fill_viridis_d(option = "viridis", direction = -1, name = "Scenario") +
  guides(
    colour = guide_legend(override.aes = list(linewidth = 1.0, alpha = 1)),
    fill   = guide_legend(override.aes = list(alpha = 0.9))
  )

print(gg_line)

# ============================================================
# 2) Apply to your EXISTING boxplot object (keep your code)
#    Example assumes your object is called: gg_box
# ============================================================

gg_box <- gg_box +
  labs(caption = caption_text) +
  theme_paper_serif +
  theme(
    panel.grid.major.x = element_blank(),
    axis.line = element_line(colour = "grey70", linewidth = 0.25)
  ) +
  # if your boxplot uses fill for scenarios/thresholds:
  scale_fill_viridis_d(option = "viridis", direction = -1, name = "Scenario") +
  scale_colour_viridis_d(option = "viridis", direction = -1, name = "Scenario")

print(gg_box)

# ============================================================
# Optional: save with consistent sizing
# ============================================================
ggsave("robustness_line_paper.png", gg_line, width = 12, height = 6, dpi = 300)
ggsave("robustness_box_paper.png",  gg_box,  width = 12, height = 6, dpi = 300)
