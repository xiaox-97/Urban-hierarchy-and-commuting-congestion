library(fixest)
library(dplyr)
library(tidyr)
library(sf)
library(data.table)
setwd('D:/urban_hierarchy_congestion')

city <- "newyork"
od   <- fread(paste0("data/od_tables/", city, "_full_od_time.csv"))
d <- od[, .(Dj = sum(flow, na.rm = TRUE)), by = d_id] %>%
  mutate(node_id = paste0("zone:", d_id))
od2 <- od %>%
  filter(o_id != d_id,
         !is.na(rdist)) %>% 
  inner_join(d %>% select(d_id, Dj), by = "d_id") %>%
  mutate(
    o_id      = factor(o_id),
    d_id      = factor(d_id),
    log_Dj    = log(Dj),
    log_rdist = log(rdist))

model <- fepois(
  flow ~ log_rdist + rdist + log_Dj | o_id,
  data = od2,
  cluster = ~ o_id + d_id
)
summary(model)

od2$flow_pred  <- predict(model, type = "response")
write.csv(od2[,c('o_id','d_id','flow','rdist','flow_pred')], paste0("results/gravity_model_fit/",city,"_full_od_gravity_fit.csv"), row.names = FALSE)
