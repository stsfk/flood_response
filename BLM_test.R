if (!require("pacman")) {
  install.packages("pacman")
}

pacman::p_load(
  tidyverse,
  sf,
  zeallot,
  tidymodels,
  stringr,
  modelr,
  rethinking
)

# data --------------------------------------------------------------------

load("./data/raw_forcing.Rda")

catchment_info <- sf::read_sf("./data/7_final_events/shapefile/9_final_station.shp") %>%
  select(-area) 

catchment_names <- dir("./data/7_final_events/data/") %>%
  str_extract("[0-9]+")


# processing --------------------------------------------------------------

LABs <- catchment_info %>%
  st_drop_geometry()%>%
  count(LAB)%>%
  arrange(desc(n))

temp <- catchment_info %>%
  st_drop_geometry() %>%
  select(catchment_name = grdc_no, area_new, elevation, Sum_CAP_MC, LAB) %>%
  mutate(catchment_name = as.character(catchment_name),
         elevation = as.numeric(elevation),
         Sum_CAP_MC = as.numeric(Sum_CAP_MC))

data_process <- data_raw %>% 
  left_join(temp, by = "catchment_name") # joining catchment area and elevation



# analysis ----------------------------------------------------------------

data_process$catchment_name %>% unique() %>% length() # n catchment = 2150

data_process$IA %>% unique() %>% length() # IA = 48945


data_process2 <- data_process %>%
  group_by(catchment_name) %>%
  mutate(flood_peak_rank = rank(-flood_peak)) %>%
  ungroup() %>%
  filter(flood_peak_rank == 1) %>%
  mutate(flood_peak_n = flood_peak/area_new,
         sum_pcp_bool = Sum_CAP_MC > 0)


data_process2 %>%
  ggplot(aes(IA, flood_increasing_rate, color = sum_pcp_bool))+
  geom_point() +
  facet_wrap(~LAB, scales = "free")


data_process2 %>%
  ggplot(aes(IA, flood_increasing_rate, color = sum_pcp_bool))+
  geom_point() +
  facet_grid(LAB~sum_pcp_bool, scales = "free")


data_process2 %>%
  filter(LAB %in% c("CEU")) %>%
  ggplot(aes(IA, `runoff depth` - initial_dis/area_new, color = sum_pcp_bool))+
  geom_point() +
  facet_grid(LAB~sum_pcp_bool, scales = "free")



data_process %>%
  filter(LAB %in% c("CEU"),
         runoff_coe!=0) %>%
  ggplot(aes(IA, ))+
  geom_point()



data_process2 %>%
  filter(LAB %in% c("CEU")) %>%
  ggplot(aes(IA, flood_lag_time, color = area_new))+
  geom_point()+
  facet_wrap(~sum_pcp_bool)


data_process %>%
  filter(LAB %in% c("CEU")) %>%
  ggplot(aes(IA, flood_increasing_rate, color = max_rain))+
  geom_point()


data_process %>%
  filter(LAB %in% c("CEU")) %>%
  mutate(sum_pcp_bool = Sum_CAP_MC > 0.1)  %>%
  ggplot(aes(IA, flood_peak/area_new, color = Sum_CAP_MC/area_new))+
  geom_point()+
  scale_colour_gradientn(colours = terrain.colors(10))+
  facet_wrap(~sum_pcp_bool)

data_process %>%
  filter(LAB %in% c("CEU"),
         year > 2000) %>%
  mutate(sum_pcp_bool = Sum_CAP_MC > 0.1)  %>%
  ggplot(aes(IA, flood_increasing/area_new, color = Sum_CAP_MC/area_new))+
  geom_point()+
  scale_colour_gradientn(colours = terrain.colors(10))+
  facet_wrap(~sum_pcp_bool)


data_process %>%
  filter(LAB %in% c("CEU"),
         year >= 1985,
         `runoff depth` > 0,
         runoff_coe < 5,
         catchment_name == "6123150") %>%
  group_by(catchment_name, year) %>%
  ungroup() %>%
  ggplot(aes(year, IA, group =catchment_name, color = runoff_coe)) +
  geom_point()+
  scale_colour_gradientn(colours = terrain.colors(10))+
  labs(title = "IA of each catchment change with time")









data_process %>%
  group_by(catchment_name, year) %>%
  summarise(IA = mean(IA),
            LAB = LAB[[1]]) %>%
  ungroup() %>%
  ggplot(aes(year, IA, group =catchment_name)) +
  geom_line()+
  facet_wrap(~LAB)+
  labs(title = "IA of each catchment change with time")



data_process %>%
  filter(year >= 1990) %>%
  mutate(flood_peak_n = flood_peak)




data_process %>%
  filter(LAB %in% c("CEU"), year >= 2015) %>%
  mutate(sum_pcp_bool = Sum_CAP_MC > 0.1)  %>%
  ggplot(aes(IA, flood_increasing/area_new, color = Sum_CAP_MC/area_new))+
  geom_point()+
  scale_colour_gradientn(colours = terrain.colors(10))+
  facet_wrap(~sum_pcp_bool)



data_process %>%
  ggplot(aes(sum_rain, sum_pcp))+
  geom_point()
  






data_process <- data_raw %>%
  filter(catchment_name %in% catchment_name_analysis) %>%
  select(catchment_name:day_e, initial_dis, flood_peak, rain_duration_time:barren_land)

temp <- catchment_info %>%
  st_drop_geometry() %>%
  select(catchment_name = grdc_no, area_new, elevation, Sum_CAP_MC) %>%
  mutate(catchment_name = as.character(catchment_name),
         elevation = as.numeric(elevation),
         Sum_CAP_MC = as.numeric(Sum_CAP_MC))

data_process <- data_process %>% left_join(temp, by = "catchment_name") # joining catchment area and elevation





