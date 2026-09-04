library(cppRouting)
library(tidyr)
library(dplyr)
library(readr)
library(sf)
library(reshape2)
library(fixest)
library(DescTools)
library(data.table)
RcppParallel::setThreadOptions(numThreads = 12)
setwd('D:/urban_hierarchy_congestion')

run_mode_choice_ue_iter_4modes <- function(
    od,
    drive_edges,
    drive_sgr,
    o, d,
    asc_car,
    asc_walk,
    asc_bike,
    beta,
    phi_bus = 0.25,
    flow_col  = "flow",
    drive_col = "drive_time",
    pt_col    = "pt_time",
    walk_col  = "walk_time",
    bike_col  = "bike_time",
    from_col  = "o_node_id",
    to_col    = "d_node_id",
    scale_factor = 1 / 4,
    max_outer_it = 15,
    ue_max_gap = 0.001,
    ue_max_it  = 30,
    ue_algorithm = "bfw",
    aon_method = "d",
    clip_u = 50,
    tol_att = 0.001,
    tol_share = 1e-4,
    use_msa = TRUE,
    demand_eps = 1e-8,
    walk_time_cap = NA_real_,
    bike_time_cap = NA_real_,
    verbose = TRUE
) {
  
  eps <- 1e-9
  
  EF_CO2 <- function(v_kmh) {
    
    # km/h → mph
    v_mph <- v_kmh / 1.609
    
    ef_mile <- ifelse(
      v_mph <= 30,
      -0.1152 * v_mph^3 + 7.9865 * v_mph^2 - 187.62 * v_mph + 1943.3,
      ifelse(
        v_mph <= 60,
        -0.0011 * v_mph^3 + 0.1333 * v_mph^2 - 6.6627 * v_mph + 509.17,
        350
      )
    )
    
    # g/mile → g/km
    ef_km <- ef_mile / 1.609
    
    return(ef_km)
  }
  
  softmax_prob <- function(U_mat) {
    U_max <- apply(U_mat, 1, max)
    U_shift <- U_mat - U_max
    expU <- exp(U_shift)
    denom <- rowSums(expU)
    expU / denom
  }
  
  clip_utility <- function(x, clip_u) {
    pmin(pmax(x, -clip_u), clip_u)
  }
  
  # =========================================================
  # Basic checks
  # =========================================================
  n <- nrow(od)
  if (n == 0) stop("od is empty.")
  
  cols_needed <- c(flow_col, drive_col, pt_col, walk_col, bike_col, from_col, to_col)
  missing <- setdiff(cols_needed, colnames(od))
  if (length(missing) > 0) {
    stop(paste("od missing columns:", paste(missing, collapse = ", ")))
  }
  
  if (!is.finite(phi_bus) || phi_bus < 0 || phi_bus > 1) {
    stop("phi_bus must be within [0,1].")
  }
  if (!is.finite(beta)) stop("beta must be finite.")
  if (!is.finite(asc_car) || !is.finite(asc_walk) || !is.finite(asc_bike)) {
    stop("ASC parameters must be finite.")
  }
  
  od[[flow_col]][!is.finite(od[[flow_col]])] <- 0
  od[[flow_col]][od[[flow_col]] < 0] <- 0
  
  if (sum(od[[flow_col]]) <= 0) stop("Total flow is zero.")
  
  # =========================================================
  # Fixed / initial times
  # =========================================================
  drive_time_ff  <- as.numeric(od[[drive_col]])
  pt_time_ideal  <- as.numeric(od[[pt_col]])
  walk_time_base <- as.numeric(od[[walk_col]])
  bike_time_base <- as.numeric(od[[bike_col]])
  
  drive_time_curr <- drive_time_ff
  
  # =========================================================
  # Output containers
  # =========================================================
  history <- data.frame(
    iter = integer(),
    phi_bus = numeric(),
    share_car = numeric(),
    share_pt = numeric(),
    share_walk = numeric(),
    share_bike = numeric(),
    total_demand_car_peak = numeric(),
    n_od_used = integer(),
    FFT = numeric(),
    ATT = numeric(),
    TTI = numeric(),
    TD = numeric(),
    ue_gap = numeric(),
    att_rel = numeric(),
    share_change_max = numeric(),
    stringsAsFactors = FALSE
  )
  
  ue_data <- NULL
  ATT_prev <- NA_real_
  
  P_car  <- rep(NA_real_, n)
  P_pt   <- rep(NA_real_, n)
  P_walk <- rep(NA_real_, n)
  P_bike <- rep(NA_real_, n)
  
  total_demand_car_peak <- NA_real_
  ATT <- NA_real_
  FFT <- NA_real_
  TTI <- NA_real_
  TD  <- NA_real_
  CO2 <- NA_real_
  gap_val <- NA_real_
  outer_converged <- FALSE
  outer_iterations <- 0L
  final_drive_time <- drive_time_curr
  final_att_rel <- NA_real_
  final_share_change <- NA_real_
  share_prev <- rep(NA_real_, 4)
  
  share_car_k  <- NA_real_
  share_pt_k   <- NA_real_
  share_walk_k <- NA_real_
  share_bike_k <- NA_real_
  
  # =========================================================
  # Outer iteration
  # =========================================================
  for (k in seq_len(max_outer_it)) {
    outer_iterations <- k
    
    # -------------------------------------------------------
    # 0) Effective times
    # -------------------------------------------------------
    pt_time_eff <- pt_time_ideal + phi_bus * pmax(0, drive_time_curr - drive_time_ff)
    
    t_car  <- drive_time_curr
    t_pt   <- pt_time_eff
    t_walk <- walk_time_base
    t_bike <- bike_time_base
    
    # -------------------------------------------------------
    # 1) Availability
    # -------------------------------------------------------
    avail_car  <- is.finite(t_car)  & (t_car  >= 0)
    avail_pt   <- is.finite(t_pt)   & (t_pt   >= 0)
    avail_walk <- is.finite(t_walk) & (t_walk >= 0)
    avail_bike <- is.finite(t_bike) & (t_bike >= 0)
    
    if (is.finite(walk_time_cap)) {
      avail_walk <- avail_walk & (t_walk <= walk_time_cap)
    }
    if (is.finite(bike_time_cap)) {
      avail_bike <- avail_bike & (t_bike <= bike_time_cap)
    }
    
    any_avail <- avail_car | avail_pt | avail_walk | avail_bike
    if (!all(any_avail)) {
      stop("Some OD pairs have no available modes.")
    }
    
    # -------------------------------------------------------
    # 2) Utilities
    # PT is the base mode with ASC = 0
    # -------------------------------------------------------
    U <- matrix(-1e12, nrow = n, ncol = 4)
    colnames(U) <- c("car", "pt", "walk", "bike")
    
    U_car_raw  <- asc_car  + beta * t_car
    U_pt_raw   <- 0.0      + beta * t_pt
    U_walk_raw <- asc_walk + beta * t_walk
    U_bike_raw <- asc_bike + beta * t_bike
    
    U[avail_car,  "car"]  <- clip_utility(U_car_raw[avail_car],   clip_u)
    U[avail_pt,   "pt"]   <- clip_utility(U_pt_raw[avail_pt],     clip_u)
    U[avail_walk, "walk"] <- clip_utility(U_walk_raw[avail_walk], clip_u)
    U[avail_bike, "bike"] <- clip_utility(U_bike_raw[avail_bike], clip_u)
    
    # -------------------------------------------------------
    # 3) Multinomial logit
    # -------------------------------------------------------
    P <- softmax_prob(U)
    
    P_car  <- P[, "car"]
    P_pt   <- P[, "pt"]
    P_walk <- P[, "walk"]
    P_bike <- P[, "bike"]
    
    P_car[!is.finite(P_car)]   <- 0
    P_pt[!is.finite(P_pt)]     <- 0
    P_walk[!is.finite(P_walk)] <- 0
    P_bike[!is.finite(P_bike)] <- 0
    
    P_sum <- P_car + P_pt + P_walk + P_bike
    bad_prob <- !is.finite(P_sum) | (P_sum <= 0)
    if (any(bad_prob)) stop("Invalid mode-choice probabilities encountered.")
    
    P_car  <- P_car  / P_sum
    P_pt   <- P_pt   / P_sum
    P_walk <- P_walk / P_sum
    P_bike <- P_bike / P_sum
    
    # -------------------------------------------------------
    # 4) Car demand for UE
    # -------------------------------------------------------
    demand_car_peak <- od[[flow_col]] * P_car * scale_factor
    demand_car_peak[!is.finite(demand_car_peak)] <- 0
    demand_car_peak[demand_car_peak < 0] <- 0
    
    flow_sum <- sum(od[[flow_col]])
    
    share_car_k  <- sum(od[[flow_col]] * P_car)  / max(flow_sum, eps)
    share_pt_k   <- sum(od[[flow_col]] * P_pt)   / max(flow_sum, eps)
    share_walk_k <- sum(od[[flow_col]] * P_walk) / max(flow_sum, eps)
    share_bike_k <- sum(od[[flow_col]] * P_bike) / max(flow_sum, eps)
    
    keep_idx <- is.finite(demand_car_peak) & (demand_car_peak > demand_eps)
    n_od_used <- sum(keep_idx)
    total_demand_car_peak <- sum(demand_car_peak[keep_idx])
    
    if (n_od_used == 0 || total_demand_car_peak <= 0) {
      stop("No OD pairs with positive driving demand after mode choice.")
    }
    
    od_use <- od[keep_idx, , drop = FALSE]
    demand_use <- demand_car_peak[keep_idx]
    
    # -------------------------------------------------------
    # 5) Free-flow time for car users
    # -------------------------------------------------------
    FFT <- sum(demand_use * drive_time_ff[keep_idx]) / max(total_demand_car_peak, eps)
    
    # -------------------------------------------------------
    # 6) UE assignment
    # -------------------------------------------------------
    ue <- assign_traffic(
      Graph = drive_sgr,
      from = od_use[[from_col]],
      to   = od_use[[to_col]],
      demand = demand_use,
      max_gap = ue_max_gap,
      max_it  = ue_max_it,
      algorithm = ue_algorithm,
      aon_method = aon_method,
      verbose = verbose
    )
    
    ue_data <- ue$data
    
    ATT <- sum(ue$data$flow * ue$data$cost) / max(total_demand_car_peak, eps)
    TTI <- ATT / max(FFT, eps)
    TD  <- sum(ue$data$flow * (ue$data$cost - ue$data$ftt))
    
    ue_data <- ue$data %>%
      inner_join(
        drive_edges %>% select(u, v, length, highway),
        by = c("from" = "u", "to" = "v")
      )
    
    speed_kmh <- ue_data$length / pmax(ue_data$cost, eps) * 60
    
    CO2 <- sum(EF_CO2(speed_kmh) * ue_data$flow * ue_data$length, na.rm = TRUE)
    
    if (!"ftt" %in% names(ue_data)) {
      stop("ue_data does not contain ftt. Cannot calculate same-link free-flow emissions.")
    }
    
    speed_ff_same_link_kmh <- ue_data$length / pmax(ue_data$ftt, eps) * 60
    
    CO2_ff_same_link <- sum(EF_CO2(speed_ff_same_link_kmh) * ue_data$flow * ue_data$length, na.rm = TRUE)
    CO2_cong_excess <- CO2 - CO2_ff_same_link
    
    gap_val <- if (!is.null(ue$gap)) as.numeric(ue$gap) else NA_real_
    
    # -------------------------------------------------------
    # 7) Update congested car time
    # -------------------------------------------------------
    congested_drive_sgr <- makegraph(
      df = ue_data[, c("from", "to", "cost")],
      directed = TRUE
    )
    
    tij <- get_distance_matrix(
      Graph = congested_drive_sgr,
      from = o$node_id,
      to   = d$node_id
    )
    drive_time_new <- as.vector(t(tij))
    final_drive_time <- drive_time_new
    
    # -------------------------------------------------------
    # 8) Convergence
    # -------------------------------------------------------
    if (k == 1) {
      att_rel <- Inf
    } else {
      att_rel <- abs(ATT - ATT_prev) / pmax(abs(ATT_prev), eps)
    }
    
    share_vec <- c(share_car_k, share_pt_k, share_walk_k, share_bike_k)
    if (k == 1 || any(!is.finite(share_prev))) {
      share_change_max <- Inf
    } else {
      share_change_max <- max(abs(share_vec - share_prev))
    }
    
    # -------------------------------------------------------
    # 9) MSA update
    # -------------------------------------------------------
    if (use_msa) {
      lambda <- 1 / k
      drive_time_next <- (1 - lambda) * drive_time_curr + lambda * drive_time_new
    } else {
      drive_time_next <- drive_time_new
    }
    
    history <- rbind(
      history,
      data.frame(
        iter = k,
        phi_bus = phi_bus,
        share_car = share_car_k,
        share_pt = share_pt_k,
        share_walk = share_walk_k,
        share_bike = share_bike_k,
        total_demand_car_peak = total_demand_car_peak,
        n_od_used = n_od_used,
        FFT = FFT,
        ATT = ATT,
        TTI = TTI,
        TD = TD,
        ue_gap = gap_val,
        att_rel = att_rel,
        share_change_max = share_change_max,
        stringsAsFactors = FALSE
      )
    )
    
    if (verbose) {
      message(sprintf(
        "Iter %d: car=%.4f, pt=%.4f, walk=%.4f, bike=%.4f, ATT=%.3f, TTI=%.4f, UE gap=%.3e, att_rel=%.3e, share_change=%.3e",
        k, share_car_k, share_pt_k, share_walk_k, share_bike_k, ATT, TTI, gap_val, att_rel, share_change_max
      ))
    }
    
    final_att_rel <- att_rel
    final_share_change <- share_change_max
    
    if (att_rel < tol_att && share_change_max < tol_share) {
      outer_converged <- TRUE
      if (verbose) message("Converged (outer ATT and modal-share criteria).")
      break
    }
    
    ATT_prev <- ATT
    share_prev <- share_vec
    drive_time_curr <- drive_time_next
  }
  
  # =========================================================
  # Final OD table
  # =========================================================
  od$P_car  <- P_car
  od$P_pt   <- P_pt
  od$P_walk <- P_walk
  od$P_bike <- P_bike
  
  # Use the shortest-path OD times from the final UE assignment.
  # These correspond to the final modal demand actually assigned in the last outer iteration.
  od$con_drive_time <- final_drive_time
  od$pt_time_eff <- pt_time_ideal + phi_bus * pmax(0, final_drive_time - drive_time_ff)
  
  # =========================================================
  # Average travel times by mode
  # =========================================================
  denom_car  <- sum(od[[flow_col]] * od$P_car)
  denom_pt   <- sum(od[[flow_col]] * od$P_pt)
  denom_walk <- sum(od[[flow_col]] * od$P_walk)
  denom_bike <- sum(od[[flow_col]] * od$P_bike)
  
  ATT_car <- if (denom_car > 0) {
    sum(od[[flow_col]] * od$P_car * od$con_drive_time) / denom_car
  } else {
    NA_real_
  }
  
  PTT <- if (denom_pt > 0) {
    sum(od[[flow_col]] * od$P_pt * od$pt_time_eff) / denom_pt
  } else {
    NA_real_
  }
  
  WTT <- if (denom_walk > 0) {
    sum(od[[flow_col]] * od$P_walk * od[[walk_col]]) / denom_walk
  } else {
    NA_real_
  }
  
  BTT <- if (denom_bike > 0) {
    sum(od[[flow_col]] * od$P_bike * od[[bike_col]]) / denom_bike
  } else {
    NA_real_
  }
  
  # =========================================================
  # System-wide average travel time
  # =========================================================
  TTT <- (
    sum(od[[flow_col]] * od$P_car  * od$con_drive_time) +
      sum(od[[flow_col]] * od$P_pt   * od$pt_time_eff) +
      sum(od[[flow_col]] * od$P_walk * od[[walk_col]]) +
      sum(od[[flow_col]] * od$P_bike * od[[bike_col]])
  ) / max(sum(od[[flow_col]]), eps)
  
  DIST <- if ("rdist" %in% names(od)) {
    sum(od[[flow_col]] * od[["rdist"]]) / max(sum(od[[flow_col]]), eps)
  } else {
    NA_real_
  }
  
  # Total free-Flow time for all modes
  TTT_ff_like <- (
    sum(od[[flow_col]] * od$P_car  * drive_time_ff) +
      sum(od[[flow_col]] * od$P_pt   * pt_time_ideal) +
      sum(od[[flow_col]] * od$P_walk * od[[walk_col]]) +
      sum(od[[flow_col]] * od$P_bike * od[[bike_col]])
  ) / max(sum(od[[flow_col]]), eps)
  
  TTI_allmode <- TTT / max(TTT_ff_like, eps)
  
  return(list(
    summary = tibble(
      share_car = share_car_k,
      share_pt = share_pt_k,
      share_walk = share_walk_k,
      share_bike = share_bike_k,
      total_demand_car_peak = total_demand_car_peak,
      ATT = ATT,
      ATT_car = ATT_car,
      TTI = TTI,
      TTI_allmode = TTI_allmode,
      TD = TD,
      PTT = PTT,
      WTT = WTT,
      BTT = BTT,
      TTT = TTT,
      DIST = DIST,
      CO2 = CO2,
      CO2_cong_excess = CO2_cong_excess,
      UE_gap = gap_val,
      UE_iterations = if (!is.null(ue$iteration)) as.integer(ue$iteration) else NA_integer_,
      outer_iterations = outer_iterations,
      outer_converged = outer_converged
    ),
    od = od,
    history = history,
    ue_data = ue_data
  ))
}

observed_traffic_assignment <- function(
    city,
    unit_name,
    asc_car,
    asc_walk,
    asc_bike,
    beta,
    phi_bus = 0.25,
    flow_col = "flow",
    ue_max_gap = 0.001,
    ue_max_it  = 30,
    tol_att = 0.001,
    max_outer_it = 15,
    use_msa = TRUE,
    verbose = TRUE
){
  od <- fread(paste0("data/od_tables/",city,"_full_od_time.csv"))
  unit <- st_read(paste0("data/taz/",city,"_",unit_name,".shp"),quiet=TRUE)
  unit <- unit %>%
    mutate(
      node_id = paste0("zone:", id)
    )
  od$flow[od$o_id == od$d_id] <- 0
  od <- od %>%
    mutate(
      o_node_id = paste0("zone:", o_id),
      d_node_id = paste0("zone:", d_id)
    )
  o <- od %>%
    group_by(o_id) %>%
    summarise(outflow = sum(flow), .groups = "drop") %>%
    rename(id = o_id)
  o <- o %>%
    inner_join(unit %>% st_drop_geometry() %>% select(id, node_id), by = "id")
  d <- od %>%
    group_by(d_id) %>%
    summarise(inflow = sum(flow), .groups = "drop") %>%
    rename(id = d_id)
  d <- d %>%
    inner_join(unit %>% st_drop_geometry() %>% select(id, node_id), by = "id")
  drive_edges <- fread(paste0("data/transport_network/table_data/", city, "/drive_edges.csv"))
  drive_edges[, alpha := 0.6]
  drive_sgr <- makegraph(
    df = drive_edges[, .(u, v, fft)],
    directed = TRUE,
    capacity = drive_edges$capacity,
    alpha = drive_edges$alpha,
    beta = drive_edges$beta
  )
  result <- run_mode_choice_ue_iter_4modes(od,
                                           drive_edges,
                                           drive_sgr,
                                           o,d,
                                           asc_car = asc_car,
                                           asc_walk = asc_walk,
                                           asc_bike = asc_bike,
                                           beta = beta,
                                           phi_bus = phi_bus,
                                           flow_col = flow_col,
                                           ue_max_gap = ue_max_gap,
                                           ue_max_it  = ue_max_it,
                                           tol_att = tol_att,
                                           max_outer_it = max_outer_it,
                                           use_msa = use_msa,
                                           verbose = verbose)
  #write.csv(result$history,paste0("results/traffic_assignment_results/",city,"/observed_traffic_assignment_history.csv"),row.names = FALSE)
  write.csv(result$ue_data,paste0("results/traffic_assignment_results/",city,"/observed_commuting_traffic_assignment.csv"),row.names = FALSE)
  write.csv(result$od,paste0("results/traffic_assignment_results/",city,"/observed_commuting_flows_congested_time.csv"),row.names = FALSE)
}

observed_traffic_assignment("beijing","grid_1k",-2.651,1.661,2.080,-0.128)
observed_traffic_assignment("shanghai","grid_1k",-2.521,1.515,1.933,-0.128)
observed_traffic_assignment("shenzhen","grid_1k",-2.347,1.393,1.807,-0.128)
observed_traffic_assignment("london","msoa",-1.437,3.939,3.163,-0.128)
observed_traffic_assignment("losangeles","tract",0.451,2.343,2.373,-0.062)
observed_traffic_assignment("newyork","tract",-2.099,1.450,0.640,-0.056)