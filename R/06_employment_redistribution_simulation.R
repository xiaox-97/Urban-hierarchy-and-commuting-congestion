# Employment-redistribution simulation
#
# Main design:
# 1. Observed destination commuting inflows represent baseline workplace employment.
# 2. Employment is redistributed proportionally from every eligible source-level unit.
# 3. Receiving units obtain employment in proportion to their existing employment.
# 4. Redistribution proportions are restricted to 0-30%.
# 5. The sequential pathway is L1 -> L2 -> L3 -> L4.
# 6. The stage-specific source-level employment share (prop) is the policy intensity.
# 7. Total input employment is conserved after every redistribution step.
# 8. An origin-constrained gravity model generates model-implied destination employment.
# 9. UE uses max_gap=0.001, max_it=30; modal split-UE uses max 15 outer iterations.
# 10. Final UE gap and outer convergence are retained as diagnostics, not scenario filters.
# 11. Agglomeration is calculated from destination employment; exponential weights decline to 0.01 at 5 km in the main analysis.

library(cppRouting)
library(sf)
library(dplyr)
library(purrr)
library(tidyr)
library(stringr)
library(readr)
library(reshape2)
library(fixest)
library(DescTools)
library(data.table)
RcppParallel::setThreadOptions(numThreads = 12)
setwd('D:/urban_hierarchy_congestion')


prepare_data <- function(city, unit_name) {
  # ---------------------------------------------------------
  # 1. OD data
  # ---------------------------------------------------------
  od <- fread(paste0("data/od_tables/", city, "_full_od_time.csv"))
  setDT(od)
  
  # OD cost table for gravity model
  od_cost <- od[, .(
    dist = first(rdist)
  ), by = .(o_id, d_id)]
  
  od_cost <- od_cost[!is.na(dist) & dist > 0]
  
  # origin totals
  o <- od[, .(Oi = sum(flow, na.rm = TRUE)), by = o_id] %>%
    mutate(node_id = paste0("zone:", o_id))
  
  # destination totals
  d <- od[, .(Dj = sum(flow, na.rm = TRUE)), by = d_id] %>%
    mutate(node_id = paste0("zone:", d_id))
  
  # ---------------------------------------------------------
  # 2. Hierarchy shapefile
  # ---------------------------------------------------------
  unit <- st_read(paste0("results/hierarchy_identification_results/", city, "_hierarchy.shp"), quiet = TRUE)
  unit$id <- as.integer(unit$id)
  
  unit_attr <- unit %>% st_drop_geometry()
  
  # Attach hierarchy attributes to the destination table.
  # The hierarchy shapefile name is retained for compatibility, but POI fields
  # are no longer used in the redistribution model.
  d <- merge(
    d,
    as.data.table(unit_attr[, c("id", "area", "level")]),
    by.x = "d_id",
    by.y = "id",
    all.x = TRUE
  )
  d$id <- as.character(d$d_id)

  # Observed destination commuting inflow is used as the baseline proxy for
  # workplace employment. Employment is the scenario-state attraction variable;
  # Dj_obs is retained as an unchanged baseline diagnostic.
  d$Dj_obs <- as.numeric(d$Dj)
  d$Employment <- as.numeric(d$Dj)
  
  # ---------------------------------------------------------
  # 3. Euclidean centroid distance matrix
  #    (for AGG or diagnostics, not the gravity model)
  # ---------------------------------------------------------
  cent <- st_centroid(unit)
  coords <- st_coordinates(cent)
  
  unit_attr <- unit_attr %>%
    mutate(
      cx = coords[, 1],
      cy = coords[, 2]
    )
  
  coord_mat <- as.matrix(unit_attr[, c("cx", "cy")])
  dist_mat <- as.matrix(dist(coord_mat)) / 1000
  rownames(dist_mat) <- unit_attr$id
  colnames(dist_mat) <- unit_attr$id
  diag(dist_mat) <- sqrt(unit$area / pi)
  
  # ---------------------------------------------------------
  # 4. Road network / assignment objects
  # ---------------------------------------------------------
  drive_edges <- fread(paste0("data/transport_network/table_data/", city, "/drive_edges.csv"))
  drive_edges[, alpha := 0.6]
  
  drive_sgr <- makegraph(
    df = drive_edges[, .(u, v, fft)],
    directed = TRUE,
    capacity = drive_edges$capacity,
    alpha = drive_edges$alpha,
    beta = drive_edges$beta
  )
  
  # ---------------------------------------------------------
  # 5. Destination layer table for sparsification
  # ---------------------------------------------------------
  dest_layer <- d[, .(d_id, level)]
  
  return(list(
    od = od,
    od_cost = od_cost,
    o = o,
    d = d,
    dest_layer = dest_layer,
    unit = unit,
    drive_edges = drive_edges,
    drive_sgr = drive_sgr,
    dist_mat = dist_mat
  ))
}

calc_impedance <- function(dist, lambda, mu) {
  ifelse(
    is.finite(dist) & dist > 0,
    dist^(lambda) * exp(mu * dist),
    0
  )
}

single_constrained_gravity <- function(
    od_cost,
    origins,
    dests,
    lambda,
    mu,
    alpha,
    verbose = FALSE
) {
  od <- as.data.table(copy(od_cost))
  ori <- as.data.table(copy(origins))
  des <- as.data.table(copy(dests))

  stopifnot(all(c("o_id", "d_id", "dist") %in% names(od)))
  stopifnot(all(c("o_id", "Oi") %in% names(ori)))
  stopifnot(all(c("d_id", "Employment") %in% names(des)))

  des[, Employment := as.numeric(Employment)]
  if (any(!is.finite(des$Employment)) || any(des$Employment < 0)) {
    stop("Employment must contain finite non-negative values.")
  }

  # Merge fixed origin totals and scenario workplace employment.
  od <- merge(od, ori, by = "o_id", all.x = TRUE)
  od <- merge(od, des, by = "d_id", all.x = TRUE)

  if (any(is.na(od$Oi))) stop("Some o_id in od_cost not found in origins")
  if (any(is.na(od$Employment))) {
    stop("Some d_id in od_cost not found in destination employment table")
  }

  # Impedance and destination-attraction terms.
  od[, f_ij := calc_impedance(dist, lambda = lambda, mu = mu)]
  od[, attr_j := (Employment)^alpha]
  od[, w_ij := attr_j * f_ij]

  # Remove infeasible or numerically invalid OD pairs.
  od <- od[is.finite(w_ij) & !is.na(w_ij) & w_ij > 0]
  if (nrow(od) == 0) stop("No feasible OD pairs remain after filtering")

  o_missing <- setdiff(ori$o_id, unique(od$o_id))
  if (length(o_missing) > 0) {
    stop(sprintf(
      "Some origins have no feasible destinations: %s",
      paste(head(o_missing, 10), collapse = ", ")
    ))
  }

  # Origin-constrained normalization.
  od[, denom_i := sum(w_ij), by = o_id]
  if (any(od$denom_i <= 0 | is.na(od$denom_i))) {
    stop("Some origins have zero denominator after weighting")
  }
  od[, Tij := Oi * w_ij / denom_i]

  # Binding origin-total check.
  row_chk <- od[, .(pred_Oi = sum(Tij)), by = o_id]
  row_chk <- merge(row_chk, ori, by = "o_id", all.x = TRUE)
  row_chk[, row_rel_err := abs(pred_Oi - Oi) / pmax(Oi, 1e-12)]

  # Destination totals are model outputs, not binding constraints.
  col_chk <- od[, .(pred_Dj = sum(Tij)), by = d_id]
  des_diag_cols <- intersect(
    c("d_id", "Employment", "Dj_obs"),
    names(des)
  )
  col_chk <- merge(
    col_chk,
    unique(des[, ..des_diag_cols]),
    by = "d_id",
    all.x = TRUE
  )
  col_chk[
    ,
    col_rel_err_vs_input_employment :=
      abs(pred_Dj - Employment) / pmax(Employment, 1e-12)
  ]
  if ("Dj_obs" %in% names(col_chk)) {
    col_chk[
      ,
      col_rel_err_vs_observed :=
        abs(pred_Dj - Dj_obs) / pmax(Dj_obs, 1e-12)
    ]
  }

  meta <- list(
    model = "origin_constrained_employment_attraction_gravity",
    lambda = lambda,
    mu = mu,
    alpha = alpha,
    max_row_rel_err = max(row_chk$row_rel_err, na.rm = TRUE),
    n_pairs = nrow(od),
    n_origins = uniqueN(od$o_id),
    n_dests = uniqueN(od$d_id),
    total_observed_O = sum(ori$Oi, na.rm = TRUE),
    total_predicted_T = sum(od$Tij, na.rm = TRUE)
  )

  if (verbose) {
    message(sprintf(
      paste0(
        "Done | pairs = %d | origins = %d | dests = %d | ",
        "max row rel err = %.3e"
      ),
      meta$n_pairs,
      meta$n_origins,
      meta$n_dests,
      meta$max_row_rel_err
    ))
  }

  list(
    flows = od[
      ,
      .(
        o_id,
        d_id,
        dist,
        Oi,
        Employment,
        f_ij,
        attr_j,
        w_ij,
        Tij
      )
    ],
    meta = meta,
    row_check = row_chk,
    col_check = col_chk
  )
}

build_assignment_support <- function(
    flow_dt,
    dest_layer = NULL,
    method = c("top_k", "cum_share"),
    k = 500,
    cum_share = 0.80,
    keep_l12 = FALSE,
    rescale_origin = TRUE
) {
  method <- match.arg(method)
  
  dt <- as.data.table(copy(flow_dt))
  stopifnot(all(c("o_id", "d_id", "Tij", "Oi") %in% names(dt)))
  
  # ---------------------------------------------------------
  # 1. Optional: Merge destination level
  # ---------------------------------------------------------
  if (keep_l12) {
    if (is.null(dest_layer)) {
      stop("dest_layer must be provided when keep_l12 = TRUE")
    }
    dl <- as.data.table(copy(dest_layer))
    stopifnot(all(c("d_id", "level") %in% names(dl)))
    dt <- merge(dt, dl, by = "d_id", all.x = TRUE)
  } else {
    dt[, level := NA_integer_]
  }
  
  # ---------------------------------------------------------
  # 2. Sort by traffic volume within "origin" from highest to lowest
  # ---------------------------------------------------------
  setorder(dt, o_id, -Tij, d_id)
  
  # ---------------------------------------------------------
  # 3. Primary Filtering Rules
  # ---------------------------------------------------------
  if (method == "top_k") {
    dt[, rank_flow := seq_len(.N), by = o_id]
    dt[, keep_main := rank_flow <= k]
  }
  
  if (method == "cum_share") {
    dt[, flow_sum := sum(Tij), by = o_id]
    
    dt[, flow_share := 0]
    dt[flow_sum > 0, flow_share := Tij / flow_sum]
    
    dt[, cum_flow_share := cumsum(flow_share), by = o_id]
    dt[, keep_main := cum_flow_share <= cum_share, by = o_id]
    
    # Keep at least the first edge for each origin
    dt[, row_id_in_o := seq_len(.N), by = o_id]
    dt[row_id_in_o == 1, keep_main := TRUE]
  }
  
  # ---------------------------------------------------------
  # 4. Reserve additional L1/L2
  # ---------------------------------------------------------
  if (keep_l12) {
    dt[, keep_l12_flag := level %in% c(1, 2)]
    dt[, keep := keep_main | keep_l12_flag]
  } else {
    dt[, keep := keep_main]
  }
  
  # ---------------------------------------------------------
  # 5. Thinned-out traffic
  # ---------------------------------------------------------
  dt[, Tij_sparse := 0]
  dt[keep == TRUE, Tij_sparse := Tij]
  
  # ---------------------------------------------------------
  # 6. Relabel by origin（Used for traffic assignment）
  # ---------------------------------------------------------
  dt[, sparse_sum := sum(Tij_sparse), by = o_id]
  
  if (rescale_origin) {
    dt[, scale_factor := 0]
    dt[sparse_sum > 0, scale_factor := Oi / sparse_sum]
    dt[, Tij_assign := Tij_sparse * scale_factor]
  } else {
    dt[, scale_factor := 1]
    dt[, Tij_assign := Tij_sparse]
  }
  
  # ---------------------------------------------------------
  # 7. Origin diagnosis
  # ---------------------------------------------------------
  origin_diag <- dt[, .(
    Oi = first(Oi),
    full_sum = sum(Tij),
    sparse_sum_before_rescale = sum(Tij_sparse),
    sparse_sum_after_rescale = sum(Tij_assign),
    n_total = .N,
    n_kept = sum(keep)
  ), by = o_id]
  
  
  # ---------------------------------------------------------
  # 9. Overall diagnosis
  # ---------------------------------------------------------
  overall_diag <- dt[, .(
    od_share_retained = mean(keep),
    flow_share_retained_before_rescale = sum(Tij_sparse) / sum(Tij),
    flow_total_original = sum(Tij),
    flow_total_before_rescale = sum(Tij_sparse),
    flow_total_after_rescale = sum(Tij_assign)
  )]
  
  # ---------------------------------------------------------
  # 10. Return
  # ---------------------------------------------------------
  list(
    flows = dt[],
    origin_diag = origin_diag[],
    overall_diag = overall_diag[]
  )
}


# =========================================================
# 5. Hierarchical employment redistribution
#    - all source units contribute proportionally
#    - targets receive employment in proportion to current employment
#    - no donor-density filter and no alternative allocation methods
# =========================================================
prepare_destination_state <- function(df) {
  dt <- as.data.table(copy(df))

  needed <- c("id", "d_id", "level", "Employment")
  missing_cols <- setdiff(needed, names(dt))
  if (length(missing_cols) > 0) {
    stop(
      "Destination table is missing columns: ",
      paste(missing_cols, collapse = ", ")
    )
  }

  dt[, id := as.character(id)]
  dt[, Employment := as.numeric(Employment)]

  if (any(!is.finite(dt$Employment)) || any(dt$Employment < 0)) {
    stop("Employment must contain finite non-negative values.")
  }

  dt[]
}

validate_redistribution_prop <- function(prop, max_prop = 0.30) {
  if (length(prop) != 1 || !is.finite(prop)) {
    stop("prop must be one finite numeric value.")
  }
  if (!is.finite(max_prop) || max_prop <= 0 || max_prop > 1) {
    stop("max_prop must be within (0, 1].")
  }
  if (prop < 0 || prop > max_prop + 1e-12) {
    stop(sprintf("prop must be within [0, %.2f].", max_prop))
  }
  invisible(TRUE)
}

# Target allocation weights are based on the target employment distribution
# immediately before each redistribution step.
existing_employment_weights <- function(df, target_idx) {
  if (length(target_idx) == 0) return(numeric())

  w <- as.numeric(df$Employment[target_idx])
  w[!is.finite(w) | w < 0] <- 0

  if (sum(w) <= 0) {
    stop(
      "All eligible target units have zero employment; ",
      "proportional allocation is undefined."
    )
  }

  w / sum(w)
}

redistribute_once <- function(
    df,
    source_filter,
    target_filter,
    prop,
    group_var = NULL,
    max_prop = 0.30,
    strict_group_match = FALSE,
    conservation_tol = 1e-8
) {
  validate_redistribution_prop(prop, max_prop)
  dt <- prepare_destination_state(df)

  total_before <- sum(dt$Employment)

  # The zero-share scenario returns the current employment state.
  if (prop <= 0) return(dt[])

  source_idx_all <- which(source_filter(dt) & dt$Employment > 0)
  target_idx_all <- which(target_filter(dt))

  if (length(source_idx_all) == 0) {
    stop("No source units with positive employment were found.")
  }
  if (length(target_idx_all) == 0) {
    stop("No target units were found for this redistribution scenario.")
  }

  allocate_block <- function(source_idx, target_idx) {
    if (length(source_idx) == 0 || length(target_idx) == 0) {
      return(invisible(NULL))
    }

    moved_out <- dt$Employment[source_idx] * prop
    moved_total <- sum(moved_out)
    if (!is.finite(moved_total) || moved_total <= 0) {
      return(invisible(NULL))
    }

    target_w <- existing_employment_weights(dt, target_idx)
    moved_in <- moved_total * target_w

    set(
      dt,
      i = source_idx,
      j = "Employment",
      value = dt$Employment[source_idx] - moved_out
    )
    set(
      dt,
      i = target_idx,
      j = "Employment",
      value = dt$Employment[target_idx] + moved_in
    )

    invisible(NULL)
  }

  if (is.null(group_var)) {
    allocate_block(source_idx_all, target_idx_all)
  } else {
    if (!group_var %in% names(dt)) {
      stop("group_var not found in destination table: ", group_var)
    }

    source_group_values <- dt[[group_var]][source_idx_all]
    target_group_values <- dt[[group_var]][target_idx_all]

    if (strict_group_match && any(is.na(source_group_values))) {
      stop("Some source units have missing ", group_var, " labels.")
    }

    source_groups <- unique(source_group_values[!is.na(source_group_values)])
    target_groups <- unique(target_group_values[!is.na(target_group_values)])

    unmatched <- setdiff(source_groups, target_groups)
    if (strict_group_match && length(unmatched) > 0) {
      stop(sprintf(
        "Some source groups have no eligible target units in %s: %s",
        group_var,
        paste(head(unmatched, 20), collapse = ", ")
      ))
    }

    # With strict_group_match = FALSE, unmatched source communities remain
    # unchanged instead of causing the whole city simulation to fail.
    for (g in intersect(source_groups, target_groups)) {
      s_idx <- which(
        source_filter(dt) &
          !is.na(dt[[group_var]]) &
          dt[[group_var]] == g &
          dt$Employment > 0
      )
      t_idx <- which(
        target_filter(dt) &
          !is.na(dt[[group_var]]) &
          dt[[group_var]] == g
      )
      allocate_block(s_idx, t_idx)
    }
  }

  total_after <- sum(dt$Employment)
  conservation_error <- total_after - total_before
  scale <- max(abs(total_before), 1)

  if (abs(conservation_error) > conservation_tol * scale) {
    stop(sprintf(
      paste0(
        "Employment conservation failed: before=%.12f, ",
        "after=%.12f, error=%.3e"
      ),
      total_before,
      total_after,
      conservation_error
    ))
  }

  dt[]
}

scenario_l1_to_l2 <- function(df, prop, max_prop = 0.30) {
  redistribute_once(
    df = df,
    source_filter = function(x) x$level == 1,
    target_filter = function(x) x$level == 2,
    prop = prop,
    group_var = NULL,
    max_prop = max_prop
  )
}

scenario_l1_to_l3 <- function(df, prop, max_prop = 0.30) {
  redistribute_once(
    df = df,
    source_filter = function(x) x$level == 1,
    target_filter = function(x) x$level == 3,
    prop = prop,
    group_var = NULL,
    max_prop = max_prop
  )
}

scenario_l1_to_l4 <- function(df, prop, max_prop = 0.30) {
  redistribute_once(
    df = df,
    source_filter = function(x) x$level == 1,
    target_filter = function(x) x$level == 4,
    prop = prop,
    group_var = NULL,
    max_prop = max_prop
  )
}

scenario_l2_to_l3 <- function(df, prop, max_prop = 0.30) {
  redistribute_once(
    df = df,
    source_filter = function(x) x$level == 2,
    target_filter = function(x) x$level == 3,
    prop = prop,
    group_var = "region_1",
    max_prop = max_prop,
    strict_group_match = FALSE
  )
}

scenario_l3_to_l4 <- function(df, prop, max_prop = 0.30) {
  redistribute_once(
    df = df,
    source_filter = function(x) x$level == 3,
    target_filter = function(x) x$level == 4,
    prop = prop,
    group_var = "region_2",
    max_prop = max_prop,
    strict_group_match = FALSE
  )
}

run_scenario <- function(df, scenario_name, prop, max_prop = 0.30) {
  switch(
    scenario_name,
    L1_to_L2 = scenario_l1_to_l2(df, prop, max_prop),
    L1_to_L3 = scenario_l1_to_l3(df, prop, max_prop),
    L1_to_L4 = scenario_l1_to_l4(df, prop, max_prop),
    L2_to_L3 = scenario_l2_to_l3(df, prop, max_prop),
    L3_to_L4 = scenario_l3_to_l4(df, prop, max_prop),
    stop("Unknown scenario: ", scenario_name)
  )
}

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
  
  return(tibble(
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
  ))
}


# =========================================================
# 7. Spatial agglomeration potential with exponential decay
# =========================================================
# The impedance matrix is fixed across redistribution scenarios. In the main
# analysis it is the Euclidean centroid-distance matrix in kilometres.
# interaction_range is defined transparently: the spatial interaction weight
# reaches weight_at_range at d = interaction_range. The default therefore gives
# w(5 km) = 0.01 rather than treating 5 km as the e-folding distance.
calc_agg_exp <- function(
    d,
    impedance_mat,
    interaction_range = 5,
    weight_at_range = 0.01,
    id_col = "id",
    emp_col = "Employment",
    normalize = FALSE,
    chunk_size = 1000
) {
  if (is.null(rownames(impedance_mat)) || is.null(colnames(impedance_mat))) {
    stop("impedance_mat must have row names and column names as unit IDs.")
  }
  if (!is.finite(interaction_range) || interaction_range <= 0) {
    stop("interaction_range must be positive.")
  }
  if (!is.finite(weight_at_range) || weight_at_range <= 0 || weight_at_range >= 1) {
    stop("weight_at_range must be within (0, 1).")
  }
  if (nrow(impedance_mat) != ncol(impedance_mat) ||
      !identical(rownames(impedance_mat), colnames(impedance_mat))) {
    stop("impedance_mat must be a square matrix with matching row/column IDs.")
  }

  dt <- as.data.table(copy(d))
  dt[[id_col]] <- as.character(dt[[id_col]])
  dt[[emp_col]] <- as.numeric(dt[[emp_col]])
  dt[[emp_col]][!is.finite(dt[[emp_col]]) | dt[[emp_col]] < 0] <- 0

  mat_ids <- rownames(impedance_mat)
  E_vec <- setNames(rep(0, length(mat_ids)), mat_ids)
  common_ids <- intersect(dt[[id_col]], mat_ids)
  E_vec[common_ids] <- dt[[emp_col]][match(common_ids, dt[[id_col]])]

  total_E <- sum(E_vec)
  if (!is.finite(total_E) || total_E <= 0) return(NA_real_)

  # beta is calibrated so that exp(-beta * interaction_range) = weight_at_range.
  decay_beta <- -log(weight_at_range) / interaction_range

  n <- length(E_vec)
  agg_i <- numeric(n)
  starts <- seq.int(1L, n, by = max(1L, as.integer(chunk_size)))

  for (s in starts) {
    idx <- s:min(s + chunk_size - 1L, n)
    imp_block <- impedance_mat[idx, , drop = FALSE]
    w_block <- exp(-decay_beta * imp_block)
    w_block[!is.finite(w_block)] <- 0
    agg_i[idx] <- as.vector(w_block %*% E_vec)
  }

  numerator <- sum(E_vec * agg_i)
  if (normalize) {
    # Dimensionless index in [0, 1] when all weights are within [0, 1].
    numerator / (total_E^2)
  } else {
    # Raw employment-weighted agglomeration potential:
    # sum_i E_i * ATEM_i
    numerator
  }
}

calc_agg_multiscale <- function(
    d,
    impedance_mat,
    interaction_ranges = c(3, 5, 10),
    weight_at_range = 0.01,
    id_col = "id",
    emp_col = "Employment",
    normalize = TRUE,
    chunk_size = 1000
) {
  interaction_ranges <- unique(as.numeric(interaction_ranges))
  if (length(interaction_ranges) == 0 ||
      any(!is.finite(interaction_ranges) | interaction_ranges <= 0)) {
    stop("interaction_ranges must contain positive finite values.")
  }
  if (!is.finite(weight_at_range) || weight_at_range <= 0 || weight_at_range >= 1) {
    stop("weight_at_range must be within (0, 1).")
  }

  values <- vapply(interaction_ranges, function(rng) {
    calc_agg_exp(
      d = d,
      impedance_mat = impedance_mat,
      interaction_range = rng,
      weight_at_range = weight_at_range,
      id_col = id_col,
      emp_col = emp_col,
      normalize = normalize,
      chunk_size = chunk_size
    )
  }, numeric(1))

  names(values) <- paste0(
    "AGG_exp_", format(interaction_ranges, trim = TRUE),
    "km_w", format(weight_at_range, scientific = FALSE, trim = TRUE)
  )
  values
}


# IMPORTANT: These coefficients must be re-estimated with destination commuting
# inflow (Employment) as the attraction variable. The values below are retained
# only as editable placeholders so the script structure remains runnable.
gravity_model_parameters <- list(
  beijing = list(-1.646,-0.063,1.034),
  shanghai = list(-1.602,-0.098,1.039),
  shenzhen = list(-1.704,-0.046,0.965),
  london = list(-0.906,-0.121,1.011),
  losangeles = list(-0.821,-0.024,0.999),
  newyork = list(-0.811,-0.042,0.973)
)

mode_choice_parameters <- list(
  beijing = list(-2.651,1.661,2.080,-0.128),
  shanghai = list(-2.521,1.515,1.933,-0.128),
  shenzhen = list(-2.347,1.393,1.807,-0.128),
  london = list(-1.437,3.939,3.163,-0.128),
  losangeles = list(0.451,2.343,2.373,-0.062),
  newyork = list(-2.099,1.450,0.640,-0.056)
)


scenario_evaluation <- function(
    city,
    unit_name,
    flow_col = "Tij_assign",
    cum_share = 1,
    props = seq(0, 0.30, by = 0.05),
    max_prop = 0.30,
    ue_max_gap = 0.001,
    ue_max_it = 30,
    max_outer_it = 15,
    tol_att = 0.001,
    tol_share = 1e-4,
    agg_ranges_km = c(3, 5, 10),
    agg_main_range_km = 5,
    agg_weight_at_range = 0.01
) {
  validate_redistribution_prop(max(props), max_prop)
  if (any(props < 0)) stop("All props must be non-negative.")
  if (!agg_main_range_km %in% agg_ranges_km) {
    agg_ranges_km <- sort(unique(c(agg_ranges_km, agg_main_range_km)))
  }
  if (!is.finite(agg_weight_at_range) ||
      agg_weight_at_range <= 0 || agg_weight_at_range >= 1) {
    stop("agg_weight_at_range must be within (0, 1).")
  }

  dat <- prepare_data(city = city, unit_name = unit_name)
  d0 <- prepare_destination_state(dat$d)

  results_summary <- list()
  counter <- 1L

  evaluate_destination_state <- function(d_state, scenario, pathway, prop, configuration) {
    d_state <- as.data.table(copy(d_state))
    grav_res <- single_constrained_gravity(
      od_cost = dat$od_cost,
      origins = dat$o,
      dests = d_state[, .(d_id, Employment, Dj_obs)],
      lambda = gravity_model_parameters[[city]][[1]],
      mu = gravity_model_parameters[[city]][[2]],
      alpha = gravity_model_parameters[[city]][[3]],
      verbose = FALSE
    )
    
    # ---------------------------------------------------------
    # Model-implied workplace employment by destination
    # Use the full gravity-model flows before OD sparsification.
    # ---------------------------------------------------------
    employment_state <- grav_res$flows[
      ,
      .(
        E_model = sum(Tij, na.rm = TRUE)
      ),
      by = d_id
    ]
    
    # Construct the table used for agglomeration calculation
    agg_state <- d_state[, .(id, d_id)]
    
    agg_state[
      employment_state,
      E_model := i.E_model,
      on = "d_id"
    ]
    
    # Destinations receiving no predicted commuting flow
    agg_state[is.na(E_model), E_model := 0]
    
    # Diagnostic: input employment is an attraction measure, whereas E_model is
    # the destination inflow endogenously generated by the single-constrained model.
    agg_state[
      d_state[, .(d_id, Employment_input = Employment)],
      Employment_input := i.Employment_input,
      on = "d_id"
    ]
    employment_input_output_cor <- suppressWarnings(
      cor(
        log1p(agg_state$Employment_input),
        log1p(agg_state$E_model),
        use = "complete.obs"
      )
    )

    employment_total <- sum(
      agg_state$E_model,
      na.rm = TRUE
    )
    
    origin_total <- sum(
      dat$o$Oi,
      na.rm = TRUE
    )
    
    employment_error <-
      employment_total - origin_total
    
    employment_scale <- max(
      abs(origin_total),
      1
    )
    
    if (
      abs(employment_error) >
      1e-8 * employment_scale
    ) {
      stop(
        sprintf(
          paste0(
            "Model-implied employment is not conserved: ",
            "origin total=%.12f, ",
            "destination total=%.12f, ",
            "error=%.3e"
          ),
          origin_total,
          employment_total,
          employment_error
        )
      )
    }

    assign_res <- build_assignment_support(
      flow_dt = grav_res$flows,
      dest_layer = dat$dest_layer,
      method = "cum_share",
      cum_share = cum_share,
      keep_l12 = FALSE,
      rescale_origin = TRUE
    )

    od_full <- merge(
      dat$od,
      assign_res$flows[, .(o_id, d_id, Tij_assign)],
      by = c("o_id", "d_id"),
      all.x = TRUE
    )
    od_full[is.na(Tij_assign), Tij_assign := 0]

    result <- run_mode_choice_ue_iter_4modes(
      od = od_full,
      drive_edges = dat$drive_edges,
      drive_sgr = dat$drive_sgr,
      o = dat$o,
      d = d_state,
      asc_car = mode_choice_parameters[[city]][[1]],
      asc_walk = mode_choice_parameters[[city]][[2]],
      asc_bike = mode_choice_parameters[[city]][[3]],
      beta = mode_choice_parameters[[city]][[4]],
      flow_col = flow_col,
      max_outer_it = max_outer_it,
      ue_max_gap = ue_max_gap,
      ue_max_it = ue_max_it,
      tol_att = tol_att,
      tol_share = tol_share,
      verbose = TRUE
    ) %>%
      mutate(
        scenario = scenario,
        pathway = pathway,
        configuration = configuration,
        prop = prop,
        input_employment_total = sum(d_state$Employment, na.rm = TRUE),
        modelled_employment_total = employment_total,
        employment_input_output_cor = employment_input_output_cor,
        gravity_max_origin_rel_error = grav_res$meta$max_row_rel_err
      )

    agg_values <- calc_agg_multiscale(
      d = agg_state,
      impedance_mat = dat$dist_mat,
      interaction_ranges = agg_ranges_km,
      weight_at_range = agg_weight_at_range,
      id_col = "id",
      emp_col = "E_model",
      normalize = FALSE
    )
    for (nm in names(agg_values)) result[[nm]] <- unname(agg_values[[nm]])
    main_agg_name <- paste0(
      "AGG_exp_", agg_main_range_km,
      "km_w", format(agg_weight_at_range, scientific = FALSE, trim = TRUE)
    )
    result$AGG <- unname(agg_values[[main_agg_name]])
    result$AGG_interaction_range_km <- agg_main_range_km
    result$AGG_weight_at_range <- agg_weight_at_range

    message(sprintf(
      paste0(
        "%s | %s | p=%.2f | TTI=%.5f | ",
        "final UE gap=%.3e | outer converged=%s | AGG=%.6e"
      ),
      pathway,
      scenario,
      prop,
      result$TTI,
      result$UE_gap,
      as.character(result$outer_converged),
      result$AGG
    ))

    result
  }

  # -------------------------------------------------------
  # Baseline
  # -------------------------------------------------------
  baseline_result <- evaluate_destination_state(
    d_state = d0,
    scenario = "baseline",
    pathway = "baseline",
    prop = 0,
    configuration = "baseline"
  )
  results_summary[[counter]] <- baseline_result
  counter <- counter + 1L

  # -------------------------------------------------------
  # Level-by-level pathway: carry the TTI-optimal state forward
  # -------------------------------------------------------
  current_state <- copy(d0)
  selected_configuration <- character()
  sequential_scenarios <- c("L1_to_L2", "L2_to_L3", "L3_to_L4")

  for (sc in sequential_scenarios) {
    stage_results <- list()
    stage_states <- list()

    for (i in seq_along(props)) {
      p <- props[[i]]
      d_new <- run_scenario(
        df = current_state,
        scenario_name = sc,
        prop = p,
        max_prop = max_prop
      )

      config_label <- paste(
        c(selected_configuration, sprintf("%s=%.2f", sc, p)),
        collapse = ";"
      )
      result <- evaluate_destination_state(
        d_state = d_new,
        scenario = sc,
        pathway = "level_by_level",
        prop = p,
        configuration = config_label
      )

      results_summary[[counter]] <- result
      counter <- counter + 1L
      stage_results[[i]] <- result
      stage_states[[i]] <- d_new
    }

    stage_df <- bind_rows(stage_results)

    # UE_gap and outer_converged are retained as numerical diagnostics only.
    # They are not used as hard filters because all scenarios are evaluated
    # with the same tight UE target and the same iteration limits.
    candidates <- stage_df %>%
      filter(is.finite(TTI), is.finite(TTT))
    if (nrow(candidates) == 0) {
      stop(sprintf("No scenarios produced finite TTI and TTT for %s in %s.", sc, city))
    }

    best <- candidates %>%
      arrange(TTI, TTT, prop) %>%
      slice(1)
    best_idx <- which.min(abs(props - best$prop))
    current_state <- copy(stage_states[[best_idx]])
    selected_configuration <- c(
      selected_configuration,
      sprintf("%s=%.2f", sc, best$prop[[1]])
    )
  }

  # -------------------------------------------------------
  # Direct cross-level pathways: each starts from baseline
  # -------------------------------------------------------
  direct_scenarios <- c("L1_to_L3", "L1_to_L4")

  for (sc in direct_scenarios) {
    for (p in props) {
      d_new <- run_scenario(
        df = d0,
        scenario_name = sc,
        prop = p,
        max_prop = max_prop
      )

      result <- evaluate_destination_state(
        d_state = d_new,
        scenario = sc,
        pathway = "direct_cross_level",
        prop = p,
        configuration = sprintf("%s=%.2f", sc, p)
      )

      results_summary[[counter]] <- result
      counter <- counter + 1L
    }
  }

  summary_df <- bind_rows(results_summary)

  baseline_tti <- baseline_result$TTI[[1]]
  baseline_agg <- baseline_result$AGG[[1]]
  baseline_co2_excess <- baseline_result$CO2_cong_excess[[1]]

  summary_df <- summary_df %>%
    mutate(
      congestion_mitigation_pct = 100 * (
        (baseline_tti - 1) - (TTI - 1)
      ) / pmax(baseline_tti - 1, 1e-12),
      AGG_loss_pct = 100 * (baseline_agg - AGG) / pmax(baseline_agg, 1e-12),
      CO2_excess_reduction_pct = 100 * (
        baseline_co2_excess - CO2_cong_excess
      ) / pmax(abs(baseline_co2_excess), 1e-12)
    )

  out_dir <- file.path("results", "traffic_assignment_results", city)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  write.csv(
    summary_df,
    file.path(out_dir, "scenario_simulation_results_final_beta_0.8.csv"),
    row.names = FALSE
  )

  summary_df
}

# =========================================================
# Optional batch runner
# =========================================================
unit_names <- c(
  beijing = "grid_1k",
  shanghai = "grid_1k",
  shenzhen = "grid_1k",
  london = "msoa",
  losangeles = "tract",
  newyork = "tract"
)

run_employment_redistribution_cities <- function(
    cities = names(unit_names),
    props = seq(0, 0.30, by = 0.05),
    cum_share = 1
) {
  out <- vector("list", length(cities))
  names(out) <- cities

  for (city in cities) {
    out[[city]] <- scenario_evaluation(
      city = city,
      unit_name = unit_names[[city]],
      flow_col = "Tij_assign",
      cum_share = cum_share,
      props = props,
      ue_max_gap = 0.001,
      ue_max_it = 30,
      max_outer_it = 15,
      tol_att = 0.001,
      tol_share = 1e-4,
      agg_ranges_km = c(3, 5, 10),
      agg_main_range_km = 5,
      agg_weight_at_range = 0.01
    )
  }

  out
}

# No simulation is run automatically when this file is sourced.
# Example:
# shenzhen_results <- scenario_evaluation(
#   city = "shenzhen",
#   unit_name = "grid_1k",
#   flow_col = "Tij_assign",
#   cum_share = 1,
#   props = seq(0, 0.30, by = 0.05)
# )

