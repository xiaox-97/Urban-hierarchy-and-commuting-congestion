# ============================================================
# Origin-constrained commuting restructuring experiment
#
# Multi-level upward commuting:
#   L3 -> L1
#   L4 -> L1
#   L4 -> L2
#
# Adjacent-level upward commuting:
#   L2 -> L1
#   L3 -> L2
#   L4 -> L3
#
# A proportion of multi-level upward flows are transferred
# to adjacent-level upward destinations within the same origin.
#
# Origin demand is fixed.
# Full UE results are preserved.
# ============================================================


library(data.table)
library(sf)


project_dir <- "D:/urban_hierarchy_congestion"


source(
  file.path(
    project_dir,
    "R/06_employment_redistribution_simulation.R"
  )
)


setwd(project_dir)



# ------------------------------------------------------------
# Settings
# ------------------------------------------------------------

unit_names <- c(
  beijing = "grid_1k",
  shanghai = "grid_1k",
  shenzhen = "grid_1k",
  london = "msoa",
  osangeles = "tract",
  newyork = "tract"
)


cities <- names(unit_names)


redistribution_share <- c(
  0,
  0.10,
  0.20,
  0.30,
  0.40,
  0.50,
  0.60,
  0.70,
  0.80,
  0.90,
  1.00
)




# ------------------------------------------------------------
# Multi-level upward -> adjacent-level upward
# ------------------------------------------------------------
redistribute_multi_to_adjacent <- function(
    od,
    hierarchy,
    share
){
  
  dt <- copy(od)
  
  # Make sure flow is numeric
  dt[, flow := as.numeric(flow)]
  
  # Hierarchy levels
  dt[, o_level :=
       hierarchy$level[
         match(as.character(o_id), hierarchy$id)
       ]]
  
  dt[, d_level :=
       hierarchy$level[
         match(as.character(d_id), hierarchy$id)
       ]]
  
  dt[, gap := o_level - d_level]
  
  # Keep original flow for redistribution
  dt[, flow0 := flow]
  
  
  # ----------------------------------------------------------
  # Calculate, for each origin:
  # 1. total multi-level upward flow
  # 2. total existing adjacent-level upward flow
  # ----------------------------------------------------------
  
  dt[
    ,
    `:=`(
      multi_total =
        sum(flow0[gap >= 2], na.rm = TRUE),
      
      adjacent_total =
        sum(flow0[gap <= 1 & gap >= 0 & flow0 > 0], na.rm = TRUE)
    ),
    by = o_id
  ]
  
  
  # Only origins with existing adjacent-level upward flow
  # are eligible for restructuring
  dt[
    ,
    eligible :=
      multi_total > 0 &
      adjacent_total > 0
  ]
  
  
  # ----------------------------------------------------------
  # Remove a proportion of multi-level upward flow
  # ----------------------------------------------------------
  
  dt[
    eligible & gap >= 2,
    flow := flow0 * (1 - share)
  ]
  
  
  # ----------------------------------------------------------
  # Add exactly the same amount to adjacent-level upward flows
  # according to their original proportions
  # ----------------------------------------------------------
  
  dt[
    eligible &
      gap <= 1 &
      gap >= 0 &
      flow0 > 0,
    
    flow :=
      flow0 +
      (multi_total * share) *
      flow0 / adjacent_total
  ]
  
  
  # ----------------------------------------------------------
  # Check flow conservation by origin
  # ----------------------------------------------------------
  
  check <- dt[
    ,
    .(
      before = sum(flow0, na.rm = TRUE),
      after  = sum(flow,  na.rm = TRUE)
    ),
    by = o_id
  ]
  
  check[, diff := after - before]
  
  max_diff <- max(abs(check$diff), na.rm = TRUE)
  
  cat(
    "Total flow difference:",
    sum(dt$flow, na.rm = TRUE) -
      sum(dt$flow0, na.rm = TRUE),
    "\n"
  )
  
  cat(
    "Maximum origin flow difference:",
    max_diff,
    "\n"
  )
  
  
  # Remove temporary columns
  dt[
    ,
    c(
      "o_level",
      "d_level",
      "gap",
      "flow0",
      "multi_total",
      "adjacent_total",
      "eligible"
    ) := NULL
  ]
  
  
  dt[]
}



# ------------------------------------------------------------
# Calculate flow hierarchy shares
# ------------------------------------------------------------

calculate_flow_structure <- function(
    od,
    hierarchy
){
  
  
  dt <- copy(od)
  
  
  dt[
    ,
    o_level :=
      hierarchy$level[
        match(
          as.character(o_id),
          hierarchy$id
        )
      ]
  ]
  
  
  dt[
    ,
    d_level :=
      hierarchy$level[
        match(
          as.character(d_id),
          hierarchy$id
        )
      ]
  ]
  
  
  
  total_flow <- sum(
    dt$flow,
    na.rm=TRUE
  )
  
  
  multi_flow <- dt[
    o_level-d_level >=2,
    sum(flow,na.rm=TRUE)
  ]
  
  
  adjacent_flow <- dt[
    o_level-d_level ==1,
    sum(flow,na.rm=TRUE)
  ]
  
  
  data.table(
    
    total_flow =
      total_flow,
    
    multi_level_upward_flow =
      multi_flow,
    
    adjacent_level_upward_flow =
      adjacent_flow,
    
    
    multi_level_upward_share =
      multi_flow/total_flow,
    
    
    adjacent_level_upward_share =
      adjacent_flow/total_flow
  )
  
}







# ------------------------------------------------------------
# Run UE
# ------------------------------------------------------------

run_ue <- function(
    od,
    dat,
    city
){
  
  
  run_mode_choice_ue_iter_4modes(
    
    od = od,
    
    drive_edges =
      dat$drive_edges,
    
    drive_sgr =
      dat$drive_sgr,
    
    o =
      dat$o,
    
    d =
      dat$d,
    
    
    asc_car =
      mode_choice_parameters[[city]][[1]],
    
    asc_walk =
      mode_choice_parameters[[city]][[2]],
    
    asc_bike =
      mode_choice_parameters[[city]][[3]],
    
    beta =
      mode_choice_parameters[[city]][[4]],
    
    
    flow_col="flow",
    
    max_outer_it=15,
    
    ue_max_gap=0.001,
    
    ue_max_it=30,
    
    tol_att=0.001,
    
    tol_share=1e-4,
    
    verbose=FALSE
    
  )
  
}



# ------------------------------------------------------------
# Run one city
# ------------------------------------------------------------

run_city <- function(city){
  
  
  message(
    "\n====== ",
    city,
    " ======"
  )
  
  
  
  dat <- prepare_data(
    city,
    unit_names[[city]]
  )
  
  
  
  od0 <- copy(dat$od)
  
  
  
  # important:
  # avoid integer truncation
  od0[
    ,
    flow :=
      as.numeric(flow)
  ]
  
  
  
  hierarchy <- as.data.table(
    st_drop_geometry(dat$unit)
  )[
    ,
    .(
      id =
        as.character(id),
      
      level
    )
  ]
  
  
  
  results <- list()
  
  
  
  for(p in redistribution_share){
    
    
    
    message(
      city,
      " | redistribution = ",
      p
    )
    
    
    
    if(p==0){
      
      
      od_scenario <- copy(od0)
      
      
    }else{
      
      
      od_scenario <-
        redistribute_multi_to_adjacent(
          od0,
          hierarchy,
          p
        )
      
    }
    
    
    
    
    # flow conservation check
    
    message(
      "flow difference = ",
      sum(od_scenario$flow) -
        sum(od0$flow)
    )
    
    
    
    ue_result <- run_ue(
      od_scenario,
      dat,
      city
    )
    
    
    
    flow_structure <-
      calculate_flow_structure(
        od_scenario,
        hierarchy
      )
    
    
    
    # keep all UE outputs
    
    ue_dt <- as.data.table(
      ue_result
    )
    
    
    
    ue_dt[
      ,
      `:=`(
        
        city = city,
        
        scenario =
          ifelse(
            p==0,
            "baseline",
            "multi_to_adjacent"
          ),
        
        redistribution_share =
          p
        
      )
    ]
    
    
    
    ue_dt <- cbind(
      ue_dt,
      flow_structure
    )
    
    
    
    results[[length(results)+1]] <-
      ue_dt
    
    
    
  }
  
  
  
  out <- rbindlist(
    results,
    fill=TRUE
  )
  
  
  
  baseline_TTI <- out[
    redistribution_share==0,
    TTI
  ][1]
  
  
  baseline_TD <- out[
    redistribution_share==0,
    TD
  ][1]
  
  
  
  out[
    ,
    `:=`(
      
      TTI_excess_reduction_pct =
        100*
        (baseline_TTI-TTI)/
        (baseline_TTI-1),
      
      
      TD_reduction_pct =
        100*
        (baseline_TD-TD)/
        baseline_TD
      
    )
  ]
  
  
  
  fwrite(
    
    out,
    
    file.path(
      
      "results",
      
      "traffic_assignment_results",
      
      city,
      
      "multi_to_adjacent_and_within_commuting_restructuring_results.csv"
      
    )
    
  )
  
  
  out
  
}




# ------------------------------------------------------------
# Run all cities
# ------------------------------------------------------------

all_results <- rbindlist(
  
  lapply(
    cities,
    run_city
  ),
  
  fill=TRUE
  
)



fwrite(
  
  all_results,
  
  file.path(
    
    "results",
    
    "traffic_assignment_results",
    
    "multi_to_adjacent_and_within_commuting_restructuring_all_cities.csv"
    
  )
  
)