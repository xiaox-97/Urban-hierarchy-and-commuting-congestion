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

compute_od_time <- function(
    city,
    unit_name
){
  unit <- st_read(paste0("data/taz/",city,"_",unit_name,".shp"),quiet=TRUE)
  unit <- unit %>%
    mutate(
      node_id = paste0("zone:", id)
    )
  od <- fread(paste0("data/od_tables/",city,"_od.csv"))
  full_od <- CJ(
    o_id = unique(od$o_id),
    d_id = unique(od$d_id),
    sorted = TRUE
  )
  
  full_od <- merge(full_od, od, by = c("o_id", "d_id"), all.x = TRUE, sort = TRUE)
  full_od[is.na(flow), flow := 0]
  
  full_od <- full_od %>%
    mutate(
      o_node_id = paste0("zone:", o_id),
      d_node_id = paste0("zone:", d_id)
    )
  o <- full_od %>%
    group_by(o_id) %>%
    summarise(outflow = sum(flow), .groups = "drop") %>%
    rename(id = o_id)
  o <- o %>%
    inner_join(unit %>% select(id, node_id), by = "id")
  d <- full_od %>%
    group_by(d_id) %>%
    summarise(inflow = sum(flow), .groups = "drop") %>%
    rename(id = d_id)
  d <- d %>%
    inner_join(unit %>% select(id, node_id), by = "id")
  
  drive_edges <- read.csv(paste0("data/transport_network/table_data/",city,"/drive_edges.csv"))
  sgr <- makegraph(df = drive_edges[,c("u", "v", "length")],directed = TRUE)
  rij_matrix <- unname(as.matrix(get_distance_matrix(Graph = sgr, from = o$node_id, to = d$node_id)))
  full_od$rdist <- as.vector(t(rij_matrix))
  sgr <- makegraph(df = drive_edges[,c("u", "v", "fft")],directed = TRUE)
  tij_matrix <- unname(as.matrix(get_distance_matrix(Graph = sgr, from = o$node_id, to = d$node_id)))
  full_od$drive_time <- as.vector(t(tij_matrix))
  
  public_transport_edges <- read.csv(paste0("data/transport_network/table_data/",city,"/public_transport_edges.csv"))
  sgr <- makegraph(df = public_transport_edges[,c("u", "v", "time")],directed = TRUE)
  tij_matrix <- unname(as.matrix(get_distance_matrix(Graph = sgr, from = o$node_id, to = d$node_id)))
  full_od$pt_time <- as.vector(t(tij_matrix))
  
  bike_edges <- read.csv(paste0("data/transport_network/table_data/",city,"/bike_edges.csv"))
  sgr <- makegraph(df = bike_edges[,c("u", "v", "time")],directed = TRUE)
  tij_matrix <- unname(as.matrix(get_distance_matrix(Graph = sgr, from = o$node_id, to = d$node_id)))
  full_od$bike_time <- as.vector(t(tij_matrix))
  
  walk_edges <- read.csv(paste0("data/transport_network/table_data/",city,"/walk_edges.csv"))
  sgr <- makegraph(df = walk_edges[,c("u", "v", "time")],directed = TRUE)
  tij_matrix <- unname(as.matrix(get_distance_matrix(Graph = sgr, from = o$node_id, to = d$node_id)))
  full_od$walk_time <- as.vector(t(tij_matrix))
  
  write.csv(full_od, file = paste0("data/od_tables/",city,"_full_od_time.csv"), row.names = FALSE)
}

compute_od_time('beijing','grid_1k')
compute_od_time('shanghai','grid_1k')
compute_od_time('shenzhen','grid_1k')
compute_od_time('london','msoa')
compute_od_time('losangeles','tract')
compute_od_time('newyork','tract')
