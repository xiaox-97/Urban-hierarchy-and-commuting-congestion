# Monte Carlo robustness test for the selected employment-redistribution scenario
#
# Design:
# 1. The deterministic employment-redistribution simulation is run first.
# 2. The single best non-baseline configuration is selected from its result table.
# 3. Baseline and optimal configurations are held fixed during Monte Carlo testing.
# 4. For each origin, stochastic OD flows are drawn from a multinomial distribution
#    using the expected shares from the origin-constrained gravity model.
# 5. Baseline and optimal cases use the same seed within each simulation, providing
#    a paired common-random-numbers comparison.
# 6. The test evaluates robustness to stochastic OD realizations; it does not
#    re-optimize the redistribution configuration in each simulation.

main_script <- "D:/urban_hierarchy_congestion/R/06_employment_redistribution_simulation.R"
source(main_script)

library(data.table)
library(dplyr)
library(stringr)

# =========================================================
# 1. Standardize and select fixed Monte Carlo cases
# =========================================================
standardize_scenario_results <- function(res) {
  dt <- as.data.table(copy(res))

  required <- c("scenario", "configuration", "TTI")
  missing_cols <- setdiff(required, names(dt))
  if (length(missing_cols) > 0) {
    stop(
      "Scenario result table is missing columns: ",
      paste(missing_cols, collapse = ", ")
    )
  }

  if (!"pathway" %in% names(dt)) dt[, pathway := NA_character_]
  if (!"prop" %in% names(dt)) dt[, prop := NA_real_]

  dt[, scenario := as.character(scenario)]
  dt[, pathway := as.character(pathway)]
  dt[, configuration := as.character(configuration)]
  dt[, prop := suppressWarnings(as.numeric(prop))]

  dt[]
}

pick_mc_cases <- function(
    city,
    objective_col = "TTI",
    objective_direction = c("min", "max")
) {
  objective_direction <- match.arg(objective_direction)

  path <- file.path(
    "results",
    "traffic_assignment_results",
    city,
    "scenario_simulation_results_final.csv"
  )
  if (!file.exists(path)) stop("Cannot find: ", path)

  res <- standardize_scenario_results(fread(path))
  if (!objective_col %in% names(res)) {
    stop("Cannot find objective_col: ", objective_col)
  }

  baseline <- res[
    scenario == "baseline" |
      configuration == "baseline"
  ]
  if (nrow(baseline) == 0) {
    baseline <- res[is.na(prop) | prop == 0]
  }
  if (nrow(baseline) == 0) {
    stop("No baseline row found for ", city)
  }
  baseline <- baseline[1]

  configuration_has_positive_share <- function(x) {
    vapply(x, function(one) {
      if (is.na(one) || !nzchar(one) || one == "baseline") return(FALSE)
      parts <- strsplit(one, ";", fixed = TRUE)[[1]]
      vals <- suppressWarnings(
        as.numeric(
          vapply(
            strsplit(parts, "=", fixed = TRUE),
            function(kv) if (length(kv) == 2) trimws(kv[2]) else NA_character_,
            character(1)
          )
        )
      )
      any(is.finite(vals) & vals > 0)
    }, logical(1))
  }

  candidates <- res[
    scenario != "baseline" &
      configuration != "baseline" &
      !is.na(configuration) &
      nzchar(configuration) &
      is.finite(get(objective_col))
  ]
  candidates <- candidates[
    configuration_has_positive_share(configuration)
  ]
  if (nrow(candidates) == 0) {
    stop("No finite non-baseline redistribution candidates found for ", city)
  }

  if (objective_direction == "min") {
    setorderv(
      candidates,
      cols = c(objective_col, "TTT", "prop"),
      order = c(1, 1, 1),
      na.last = TRUE
    )
  } else {
    setorderv(
      candidates,
      cols = c(objective_col, "TTT", "prop"),
      order = c(-1, 1, 1),
      na.last = TRUE
    )
  }
  optimal <- candidates[1]

  metric_value <- function(row, nm) {
    if (nm %in% names(row)) row[[nm]][1] else NA_real_
  }

  rbindlist(
    list(
      data.table(
        city = city,
        case_type = "baseline",
        scenario = baseline$scenario[1],
        pathway = baseline$pathway[1],
        configuration = "baseline",
        selected_TTI = metric_value(baseline, "TTI"),
        selected_TTT = metric_value(baseline, "TTT"),
        selected_CO2 = metric_value(baseline, "CO2"),
        selected_CO2_cong_excess = metric_value(
          baseline,
          "CO2_cong_excess"
        ),
        selected_AGG = metric_value(baseline, "AGG")
      ),
      data.table(
        city = city,
        case_type = "optimal",
        scenario = optimal$scenario[1],
        pathway = optimal$pathway[1],
        configuration = optimal$configuration[1],
        selected_TTI = metric_value(optimal, "TTI"),
        selected_TTT = metric_value(optimal, "TTT"),
        selected_CO2 = metric_value(optimal, "CO2"),
        selected_CO2_cong_excess = metric_value(
          optimal,
          "CO2_cong_excess"
        ),
        selected_AGG = metric_value(optimal, "AGG")
      )
    ),
    fill = TRUE
  )
}

# =========================================================
# 2. Reconstruct the fixed employment configuration
# =========================================================
parse_configuration <- function(configuration) {
  configuration <- trimws(as.character(configuration)[1])

  if (
    is.na(configuration) ||
      configuration == "" ||
      configuration == "baseline"
  ) {
    return(data.table(
      scenario = character(),
      prop = numeric()
    ))
  }

  parts <- strsplit(configuration, ";", fixed = TRUE)[[1]]
  parts <- trimws(parts)
  parts <- parts[nzchar(parts)]

  parsed <- lapply(parts, function(part) {
    kv <- strsplit(part, "=", fixed = TRUE)[[1]]
    if (length(kv) != 2) {
      stop("Invalid configuration element: ", part)
    }

    scenario_name <- trimws(kv[1])
    prop_value <- suppressWarnings(as.numeric(trimws(kv[2])))

    if (!nzchar(scenario_name) || !is.finite(prop_value)) {
      stop("Invalid configuration element: ", part)
    }

    data.table(
      scenario = scenario_name,
      prop = prop_value
    )
  })

  rbindlist(parsed)
}

apply_case_configuration <- function(
    base_d,
    case_row,
    max_prop = 0.30
) {
  d_case <- prepare_destination_state(base_d)
  config_steps <- parse_configuration(case_row$configuration[1])

  if (nrow(config_steps) == 0) return(d_case[])

  for (i in seq_len(nrow(config_steps))) {
    d_case <- run_scenario(
      df = d_case,
      scenario_name = config_steps$scenario[i],
      prop = config_steps$prop[i],
      max_prop = max_prop
    )
  }

  prepare_destination_state(d_case)
}

# =========================================================
# 3. Build expected OD for one fixed employment case
# =========================================================
build_expected_od_for_case <- function(
    dat,
    city,
    case_row,
    assignment_method = "cum_share",
    cum_share = 1,
    keep_l12 = FALSE,
    max_prop = 0.30
) {
  d_case <- apply_case_configuration(
    base_d = dat$d,
    case_row = case_row,
    max_prop = max_prop
  )

  gpar <- gravity_model_parameters[[city]]
  if (is.null(gpar) || length(gpar) < 3) {
    stop("No complete gravity_model_parameters found for ", city)
  }

  gravity_result <- single_constrained_gravity(
    od_cost = dat$od_cost,
    origins = dat$o,
    dests = d_case[, .(d_id, Employment, Dj_obs)],
    lambda = gpar[[1]],
    mu = gpar[[2]],
    alpha = gpar[[3]],
    verbose = FALSE
  )

  support <- build_assignment_support(
    flow_dt = gravity_result$flows,
    dest_layer = dat$dest_layer,
    method = assignment_method,
    cum_share = cum_share,
    keep_l12 = keep_l12,
    rescale_origin = TRUE
  )

  flow_assign <- as.data.table(support$flows)[
    ,
    .(
      o_id,
      d_id,
      flow_expected = Tij_assign
    )
  ]

  od_template <- as.data.table(copy(dat$od))
  od_template[, od_order_mc := .I]

  od_expected <- merge(
    od_template,
    flow_assign,
    by = c("o_id", "d_id"),
    all.x = TRUE,
    sort = FALSE
  )

  origin_total <- as.data.table(dat$o)[
    ,
    .(
      o_id,
      Oi_mc = Oi
    )
  ]
  od_expected <- merge(
    od_expected,
    origin_total,
    by = "o_id",
    all.x = TRUE,
    sort = FALSE
  )

  setorder(od_expected, od_order_mc)
  od_expected[, od_order_mc := NULL]

  od_expected[is.na(flow_expected), flow_expected := 0]
  od_expected[!is.finite(flow_expected), flow_expected := 0]
  od_expected[flow_expected < 0, flow_expected := 0]

  if (!"o_node_id" %in% names(od_expected)) {
    od_expected[, o_node_id := paste0("zone:", o_id)]
  }
  if (!"d_node_id" %in% names(od_expected)) {
    od_expected[, d_node_id := paste0("zone:", d_id)]
  }

  attr(od_expected, "gravity_meta") <- gravity_result$meta
  attr(od_expected, "support_diag") <- support$overall_diag
  attr(od_expected, "employment_state") <- d_case

  od_expected[]
}

# =========================================================
# 4. Multinomial OD sampling by origin
# =========================================================
sample_multinomial_od <- function(
    od_expected,
    seed = NULL,
    expected_col = "flow_expected",
    sampled_col = "flow_mc",
    origin_total_col = "Oi_mc"
) {
  if (!is.null(seed)) set.seed(seed)

  dt <- as.data.table(copy(od_expected))
  dt[, mc_order := .I]
  setorder(dt, o_id, d_id)

  dt[!is.finite(get(expected_col)), (expected_col) := 0]
  dt[get(expected_col) < 0, (expected_col) := 0]

  if (!origin_total_col %in% names(dt)) {
    dt[
      ,
      (origin_total_col) := sum(get(expected_col), na.rm = TRUE),
      by = o_id
    ]
  }

  dt[
    ,
    (sampled_col) := {
      mu <- get(expected_col)
      total_mu <- sum(mu, na.rm = TRUE)
      size_i <- as.integer(round(first(get(origin_total_col))))

      if (!is.finite(total_mu) || total_mu <= 0 || size_i <= 0) {
        rep(0, .N)
      } else {
        prob_i <- mu / total_mu
        as.numeric(
          rmultinom(
            n = 1,
            size = size_i,
            prob = prob_i
          )[, 1]
        )
      }
    },
    by = o_id
  ]

  setorder(dt, mc_order)
  dt[, mc_order := NULL]
  dt[]
}

# =========================================================
# 5. Traffic assignment for one sampled OD
# =========================================================
run_assignment_for_sample <- function(
    dat,
    city,
    od_mc,
    phi_bus = 0.25,
    max_outer_it = 15,
    ue_max_gap = 0.001,
    ue_max_it = 30,
    tol_att = 0.001,
    tol_share = 1e-4,
    walk_time_cap = 90,
    bike_time_cap = 60
) {
  mpar <- mode_choice_parameters[[city]]
  if (is.null(mpar) || length(mpar) < 4) {
    stop("No complete mode_choice_parameters found for ", city)
  }

  out <- run_mode_choice_ue_iter_4modes(
    od = od_mc,
    drive_edges = dat$drive_edges,
    drive_sgr = dat$drive_sgr,
    o = dat$o,
    d = dat$d,
    asc_car = mpar[[1]],
    asc_walk = mpar[[2]],
    asc_bike = mpar[[3]],
    beta = mpar[[4]],
    phi_bus = phi_bus,
    flow_col = "flow_mc",
    drive_col = "drive_time",
    pt_col = "pt_time",
    walk_col = "walk_time",
    bike_col = "bike_time",
    from_col = "o_node_id",
    to_col = "d_node_id",
    scale_factor = 1 / 4,
    max_outer_it = max_outer_it,
    ue_max_gap = ue_max_gap,
    ue_max_it = ue_max_it,
    tol_att = tol_att,
    tol_share = tol_share,
    walk_time_cap = walk_time_cap,
    bike_time_cap = bike_time_cap,
    verbose = FALSE
  )

  as.data.table(out)
}

# =========================================================
# 6. Paired Monte Carlo for one city
# =========================================================
run_mc_city <- function(
    city,
    unit_name = unit_names[[city]],
    S = 100,
    seed_base = 20260523,
    objective_col = "TTI",
    objective_direction = "min",
    phi_bus = 0.25,
    assignment_method = "cum_share",
    cum_share = 1,
    keep_l12 = FALSE,
    max_prop = 0.30,
    max_outer_it = 15,
    ue_max_gap = 0.001,
    ue_max_it = 30,
    tol_att = 0.001,
    tol_share = 1e-4
) {
  if (is.null(unit_name) || is.na(unit_name)) {
    stop("No unit_name found for ", city)
  }
  if (!is.finite(S) || S < 1) stop("S must be a positive integer.")
  S <- as.integer(S)

  message("Preparing data for ", city)
  dat <- prepare_data(city, unit_name)

  cases <- pick_mc_cases(
    city = city,
    objective_col = objective_col,
    objective_direction = objective_direction
  )

  out_dir <- file.path(
    "results",
    "traffic_assignment_results",
    city
  )
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  fwrite(
    cases,
    file.path(out_dir, "monte_carlo_selected_cases.csv")
  )

  expected_list <- setNames(
    vector("list", nrow(cases)),
    cases$case_type
  )
  for (r in seq_len(nrow(cases))) {
    message(
      "Building expected OD: ",
      city,
      " | ",
      cases$case_type[r],
      " | ",
      cases$configuration[r]
    )
    expected_list[[cases$case_type[r]]] <- build_expected_od_for_case(
      dat = dat,
      city = city,
      case_row = cases[r],
      assignment_method = assignment_method,
      cum_share = cum_share,
      keep_l12 = keep_l12,
      max_prop = max_prop
    )
  }

  raw_list <- vector("list", S * nrow(cases))
  counter <- 1L

  for (s in seq_len(S)) {
    message("Paired MC run ", s, " / ", S, " | ", city)

    # The same seed is reused for baseline and optimal cases in simulation s.
    # Their OD row order and origin totals are held constant.
    paired_seed <- seed_base + s

    for (r in seq_len(nrow(cases))) {
      case_row <- cases[r]
      case_type_r <- case_row$case_type[1]

      one <- tryCatch({
        od_mc <- sample_multinomial_od(
          od_expected = expected_list[[case_type_r]],
          seed = paired_seed,
          sampled_col = "flow_mc"
        )

        res <- run_assignment_for_sample(
          dat = dat,
          city = city,
          od_mc = od_mc,
          phi_bus = phi_bus,
          max_outer_it = max_outer_it,
          ue_max_gap = ue_max_gap,
          ue_max_it = ue_max_it,
          tol_att = tol_att,
          tol_share = tol_share
        )

        res[
          ,
          `:=`(
            city = city,
            sim = s,
            paired_seed = paired_seed,
            case_type = case_type_r,
            scenario = case_row$scenario[1],
            pathway = case_row$pathway[1],
            configuration = case_row$configuration[1],
            selected_TTI = case_row$selected_TTI[1],
            selected_TTT = case_row$selected_TTT[1],
            selected_CO2 = case_row$selected_CO2[1],
            selected_CO2_cong_excess =
              case_row$selected_CO2_cong_excess[1],
            selected_AGG = case_row$selected_AGG[1],
            total_flow_expected = sum(
              od_mc$flow_expected,
              na.rm = TRUE
            ),
            total_flow_sampled = sum(
              od_mc$flow_mc,
              na.rm = TRUE
            ),
            ok = TRUE,
            error = NA_character_
          )
        ]
        res
      }, error = function(e) {
        data.table(
          city = city,
          sim = s,
          paired_seed = paired_seed,
          case_type = case_type_r,
          scenario = case_row$scenario[1],
          pathway = case_row$pathway[1],
          configuration = case_row$configuration[1],
          selected_TTI = case_row$selected_TTI[1],
          selected_TTT = case_row$selected_TTT[1],
          selected_CO2 = case_row$selected_CO2[1],
          selected_CO2_cong_excess =
            case_row$selected_CO2_cong_excess[1],
          selected_AGG = case_row$selected_AGG[1],
          total_flow_expected = NA_real_,
          total_flow_sampled = NA_real_,
          ok = FALSE,
          error = conditionMessage(e)
        )
      })

      raw_list[[counter]] <- one
      counter <- counter + 1L
    }
  }

  raw <- rbindlist(raw_list, fill = TRUE)
  fwrite(
    raw,
    file.path(out_dir, "monte_carlo_od_metrics_raw.csv")
  )

  # Construct explicit baseline-optimal pairs.
  pair_metrics <- c(
    "TTI",
    "TTI_allmode",
    "TTT",
    "ATT",
    "TD",
    "CO2",
    "CO2_cong_excess",
    "UE_gap",
    "outer_iterations",
    "outer_converged"
  )
  pair_metrics <- intersect(pair_metrics, names(raw))

  base_cols <- c("city", "sim", pair_metrics)
  optimal_cols <- c("city", "sim", pair_metrics)

  baseline_mc <- raw[
    case_type == "baseline" & ok == TRUE,
    ..base_cols
  ]
  optimal_mc <- raw[
    case_type == "optimal" & ok == TRUE,
    ..optimal_cols
  ]

  setnames(
    baseline_mc,
    old = pair_metrics,
    new = paste0(pair_metrics, "_baseline")
  )
  setnames(
    optimal_mc,
    old = pair_metrics,
    new = paste0(pair_metrics, "_optimal")
  )

  paired <- merge(
    baseline_mc,
    optimal_mc,
    by = c("city", "sim"),
    all = FALSE
  )

  if (nrow(paired) > 0) {
    paired[
      ,
      `:=`(
        delta_TTI = TTI_optimal - TTI_baseline,
        congestion_mitigated = TTI_optimal < TTI_baseline,
        congestion_mitigation_pct = 100 * (
          TTI_baseline - TTI_optimal
        ) / pmax(TTI_baseline - 1, 1e-12)
      )
    ]
  } else {
    paired[
      ,
      `:=`(
        delta_TTI = numeric(),
        congestion_mitigated = logical(),
        congestion_mitigation_pct = numeric()
      )
    ]
  }

  if (all(c("TTT_baseline", "TTT_optimal") %in% names(paired))) {
    paired[, delta_TTT := TTT_optimal - TTT_baseline]
  }
  if (
    all(
      c(
        "CO2_cong_excess_baseline",
        "CO2_cong_excess_optimal"
      ) %in% names(paired)
    )
  ) {
    paired[
      ,
      CO2_excess_reduction_pct := 100 * (
        CO2_cong_excess_baseline -
          CO2_cong_excess_optimal
      ) / pmax(abs(CO2_cong_excess_baseline), 1e-12)
    ]
  }

  fwrite(
    paired,
    file.path(out_dir, "monte_carlo_paired_effects.csv")
  )

  q_safe <- function(x, p) {
    if (length(x) == 0 || all(!is.finite(x))) return(NA_real_)
    as.numeric(quantile(x, probs = p, na.rm = TRUE, names = FALSE))
  }

  if (nrow(paired) > 0) {
    summary <- paired[
      ,
      .(
        n_requested = S,
        n_complete_pairs = .N,
        probability_of_congestion_mitigation =
          mean(congestion_mitigated, na.rm = TRUE),
        mean_delta_TTI = mean(delta_TTI, na.rm = TRUE),
        median_delta_TTI = median(delta_TTI, na.rm = TRUE),
        delta_TTI_q025 = q_safe(delta_TTI, 0.025),
        delta_TTI_q975 = q_safe(delta_TTI, 0.975),
        mean_congestion_mitigation_pct =
          mean(congestion_mitigation_pct, na.rm = TRUE),
        median_congestion_mitigation_pct =
          median(congestion_mitigation_pct, na.rm = TRUE),
        mitigation_pct_q025 =
          q_safe(congestion_mitigation_pct, 0.025),
        mitigation_pct_q975 =
          q_safe(congestion_mitigation_pct, 0.975)
      ),
      by = city
    ]
  } else {
    summary <- data.table(
      city = city,
      n_requested = S,
      n_complete_pairs = 0L,
      probability_of_congestion_mitigation = NA_real_,
      mean_delta_TTI = NA_real_,
      median_delta_TTI = NA_real_,
      delta_TTI_q025 = NA_real_,
      delta_TTI_q975 = NA_real_,
      mean_congestion_mitigation_pct = NA_real_,
      median_congestion_mitigation_pct = NA_real_,
      mitigation_pct_q025 = NA_real_,
      mitigation_pct_q975 = NA_real_
    )
  }

  summary[
    ,
    `:=`(
      optimal_scenario = cases[
        case_type == "optimal",
        scenario
      ][1],
      optimal_pathway = cases[
        case_type == "optimal",
        pathway
      ][1],
      optimal_configuration = cases[
        case_type == "optimal",
        configuration
      ][1]
    )
  ]

  fwrite(
    summary,
    file.path(out_dir, "monte_carlo_robustness_summary.csv")
  )

  list(
    selected_cases = cases,
    raw = raw,
    paired = paired,
    summary = summary
  )
}

# =========================================================
# 7. Optional all-city runner
# =========================================================
run_mc_all_cities <- function(
    cities = names(unit_names),
    S = 100,
    seed_base = 20260523,
    objective_col = "TTI",
    objective_direction = "min",
    phi_bus = 0.25
) {
  all_results <- vector("list", length(cities))
  names(all_results) <- cities

  for (city in cities) {
    all_results[[city]] <- run_mc_city(
      city = city,
      unit_name = unit_names[[city]],
      S = S,
      seed_base = seed_base,
      objective_col = objective_col,
      objective_direction = objective_direction,
      phi_bus = phi_bus
    )
  }

  raw_all <- rbindlist(
    lapply(all_results, `[[`, "raw"),
    fill = TRUE
  )
  paired_all <- rbindlist(
    lapply(all_results, `[[`, "paired"),
    fill = TRUE
  )
  summary_all <- rbindlist(
    lapply(all_results, `[[`, "summary"),
    fill = TRUE
  )

  out_dir <- file.path("results", "traffic_assignment_results")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  fwrite(
    raw_all,
    file.path(out_dir, "monte_carlo_od_metrics_raw_all_cities.csv")
  )
  fwrite(
    paired_all,
    file.path(out_dir, "monte_carlo_paired_effects_all_cities.csv")
  )
  fwrite(
    summary_all,
    file.path(out_dir, "monte_carlo_robustness_summary_all_cities.csv")
  )

  all_results
}

# No Monte Carlo simulation is run automatically.
# Example:
mc_london <- run_mc_city(
  city = "london",
  unit_name = unit_names[["london"]],
  seed_base = 20260805,
  S = 50,
  objective_col = "TTI",
  objective_direction = "min"
)
