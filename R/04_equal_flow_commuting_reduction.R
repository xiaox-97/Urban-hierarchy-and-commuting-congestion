# ============================================================
# Compare congestion responses to removing the same absolute
# commuting flow from:
#   1. multi-level upward commuting
#   2. adjacent-level upward commuting
#   3. non-upward commuting
# ============================================================

library(data.table)
library(sf)
library(ggplot2)

project_dir <- "D:/urban_hierarchy_congestion"
main_script <- "D:/urban_hierarchy_congestion/R/06_employment_redistribution_simulation.R"

source(main_script)
setwd(project_dir)

unit_names <- c(
  beijing = "grid_1k",
  shanghai = "grid_1k",
  shenzhen = "grid_1k",
  london = "msoa",
  losangeles = "tract",
  newyork = "tract"
)

cities <- names(unit_names)

# q is the removed flow as a share of the smallest of the three groups.
q_values <- c(0.25, 0.50, 0.75, 1.00)


# ------------------------------------------------------------
# Run the existing mode-choice and UE model
# ------------------------------------------------------------
run_ue <- function(od, dat, city) {

  as.data.table(
    run_mode_choice_ue_iter_4modes(
      od = od,
      drive_edges = dat$drive_edges,
      drive_sgr = dat$drive_sgr,
      o = dat$o,
      d = dat$d,
      asc_car = mode_choice_parameters[[city]][[1]],
      asc_walk = mode_choice_parameters[[city]][[2]],
      asc_bike = mode_choice_parameters[[city]][[3]],
      beta = mode_choice_parameters[[city]][[4]],
      flow_col = "flow",
      max_outer_it = 15,
      ue_max_gap = 0.001,
      ue_max_it = 30,
      tol_att = 0.001,
      tol_share = 1e-4,
      verbose = FALSE
    )
  )
}


# ------------------------------------------------------------
# One city
# ------------------------------------------------------------
run_city <- function(city) {

  message("\n===== ", city, " =====")

  dat <- prepare_data(
    city = city,
    unit_name = unit_names[[city]]
  )

  od0 <- copy(dat$od)

  hierarchy <- as.data.table(
    st_drop_geometry(dat$unit)
  )[
    ,
    .(
      id = as.character(id),
      level
    )
  ]

  # Use match rather than merge so the original OD row order is unchanged.
  od0[
    ,
    origin_level := hierarchy$level[
      match(as.character(o_id), hierarchy$id)
    ]
  ]

  od0[
    ,
    destination_level := hierarchy$level[
      match(as.character(d_id), hierarchy$id)
    ]
  ]

  od0[
    ,
    hierarchy_gap := origin_level - destination_level
  ]

  od0[
    ,
    commuting_group := fcase(
      hierarchy_gap >= 2,
      "multi_level_upward",

      hierarchy_gap == 1,
      "adjacent_level_upward",

      default = "non_upward"
    )
  ]

  group_totals <- od0[
    ,
    .(
      group_flow = sum(flow, na.rm = TRUE),
      n_od = .N
    ),
    by = commuting_group
  ]

  total_flow <- sum(od0$flow)
  common_base_flow <- min(group_totals$group_flow)

  group_totals[
    ,
    `:=`(
      city = city,
      share_of_city_flow = group_flow / total_flow,
      common_base_flow = common_base_flow
    )
  ]

  # Baseline
  baseline <- run_ue(
    od = od0,
    dat = dat,
    city = city
  )

  baseline[
    ,
    `:=`(
      city = city,
      commuting_group = "baseline",
      q = 0,
      removed_flow = 0,
      removed_share_total = 0,
      group_reduction_rate = 0
    )
  ]

  scenario_results <- list(baseline)
  counter <- 2L

  # Equal absolute flow removal
  for (q in q_values) {

    removed_flow <- q * common_base_flow

    for (group_name in c(
      "multi_level_upward",
      "adjacent_level_upward",
      "non_upward"
    )) {

      group_flow <- group_totals[
        commuting_group == group_name,
        group_flow
      ]

      reduction_rate <- removed_flow / group_flow

      od_scenario <- copy(od0)

      od_scenario[
        commuting_group == group_name,
        flow := flow * (1 - reduction_rate)
      ]

      message(
        city,
        " | ",
        group_name,
        " | q=",
        q
      )

      result <- run_ue(
        od = od_scenario,
        dat = dat,
        city = city
      )

      result[
        ,
        `:=`(
          city = city,
          commuting_group = group_name,
          q = q,
          removed_flow = removed_flow,
          removed_share_total = removed_flow / total_flow,
          group_reduction_rate = reduction_rate
        )
      ]

      scenario_results[[counter]] <- result
      counter <- counter + 1L
    }
  }

  results <- rbindlist(
    scenario_results,
    use.names = TRUE,
    fill = TRUE
  )

  baseline_TTI <- baseline$TTI
  baseline_ATT <- baseline$ATT
  baseline_TD <- baseline$TD
  baseline_CO2_excess <- baseline$CO2_cong_excess

  results[
    ,
    `:=`(
      TTI_absolute_reduction = baseline_TTI - TTI,

      TTI_excess_reduction_pct =
        100 * (baseline_TTI - TTI) /
        (baseline_TTI - 1),

      ATT_reduction_pct =
        100 * (baseline_ATT - ATT) /
        baseline_ATT,

      TD_reduction_pct =
        100 * (baseline_TD - TD) /
        baseline_TD,

      CO2_excess_reduction_pct =
        100 * (
          baseline_CO2_excess -
            CO2_cong_excess
        ) /
        baseline_CO2_excess,

      removed_share_total_pct =
        100 * removed_share_total
    )
  ]

  results[
    commuting_group != "baseline",
    response_rank_within_q :=
      frank(
        -TTI_excess_reduction_pct,
        ties.method = "min"
      ),
    by = q
  ]

  output_dir <- paste0(
    "results/traffic_assignment_results/",
    city
  )

  fwrite(
    results,
    paste0(
      output_dir,
      "/equal_flow_upward_commuting_results_5_levels.csv"
    )
  )

  fwrite(
    group_totals,
    paste0(
      output_dir,
      "/equal_flow_upward_commuting_group_totals_5_levels.csv"
    )
  )

  list(
    results = results,
    group_totals = group_totals
  )
}


# ------------------------------------------------------------
# Run all cities
# ------------------------------------------------------------
city_outputs <- lapply(
  cities,
  run_city
)

names(city_outputs) <- cities

all_results <- rbindlist(
  lapply(city_outputs, `[[`, "results"),
  use.names = TRUE,
  fill = TRUE
)

all_group_totals <- rbindlist(
  lapply(city_outputs, `[[`, "group_totals"),
  use.names = TRUE,
  fill = TRUE
)

fwrite(
  all_results,
  "results/traffic_assignment_results/equal_flow_upward_commuting_results_all_cities.csv"
)

fwrite(
  all_group_totals,
  "results/traffic_assignment_results/equal_flow_upward_commuting_group_totals_all_cities.csv"
)
