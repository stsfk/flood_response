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


# Quality checks ----------------------------------------------------------

# filter complete data
1 - (data_process %>% complete.cases() %>% sum())/nrow(data_process) 
data_process <- data_process[data_process %>% complete.cases(),]

# remove 0 runoff depth
data_process <- data_process %>%
  filter(`runoff depth` > 0)

# land use sum
data_process %>%
  transmute(land_use_sum = IA + cropland + forest + grassland + shrubland + tundra + barren_land,
            land_use_interval = cut_interval(land_use_sum, 10)) %>%
  count(land_use_interval)

data_process %>%
  transmute(land_use_sum = IA + cropland + forest + grassland + shrubland + tundra + barren_land) %>%
  filter(land_use_sum < 0.95 |
         land_use_sum > 1.05)

data_process <- data_process %>%
  mutate(land_use_sum = IA + cropland + forest + grassland + shrubland + tundra + barren_land) %>%
  filter(land_use_sum >= 0.95, land_use_sum <= 1.05) %>%
  select(-land_use_sum)

# filter year >= 1985
data_process <- data_process %>%
  filter(year >= 1985) 

# filter flood increase
data_process <- data_process %>%
  filter(flood_increasing > 0,
         flood_increasing_rate > 0,
         flood_recession_rate > 0)


# summarize ---------------------------------------------------------------

# number of flood each catchment
data_process %>% count(catchment_name)%>% arrange(n)
data_process %>% count(catchment_name)%>% arrange(desc(n))

data_process %>% count(catchment_name) %>% pull(n) %>% hist(main = "# events per catchment")

# number of flood each catchment
data_process %>% count(LAB)%>% arrange(n)
data_process %>% count(LAB)%>% arrange(desc(n))

data_process %>% count(LAB) %>% pull(n) %>% hist()

log10(data_process$`runoff depth`) %>% hist(main = "Histogram of log 10 runoff depth")

data_anual_max = data_process %>%
  group_by(year, catchment_name) %>%
  mutate(flood_peak_rank = rank(-flood_peak)) %>%
  ungroup() %>%
  filter(flood_peak_rank == 1) %>%
  select(-flood_peak_rank)

log10(data_anual_max$`runoff depth`) %>% hist(main = "Histogram of log 10 runoff depth")

# Model floods from one region --------------------------------------------

# remove small events
data_process = data_process %>%
  group_by(year, catchment_name) %>%
  mutate(flood_peak_rank = rank(-flood_peak)) %>%
  ungroup() %>%
  filter(flood_peak_rank == 1) %>%
  select(-flood_peak_rank)

ID_variables <- data_process %>%
  select(catchment_name:day_e)

catchment_id <- ID_variables %>%
  mutate(catchment_name = factor(catchment_name) %>% as.numeric()) %>%
  select(catchment_name)

flood_charater <- data_process %>%
  rename(runoff_depth="runoff depth") %>%
  mutate(flood_peak = flood_peak/area_new,
         flood_peak = log10(flood_peak),
         flood_increasing_rate = log10(flood_increasing_rate),
         flood_recession_rate = log10(flood_recession_rate),
         runoff_depth = log10(runoff_depth)) %>%
  select(starts_with("flood"), runoff_depth, sm_ave) %>%
  select(-flood_increasing)

forcing <- data_process %>%
  select(
    initial_dis,
    rain_duration_time,
    rain_pre_1day,
    sum_rain,
    max_rain,
    rain_intensity,
    sm_pre_1day,
    tem_pre_1day,
    tem_ave:sum_snowfall
  ) %>%
  mutate(
    sum_rain = log10(sum_rain),
    sum_pcp = log10(sum_pcp),
    rain_intensity = log10(rain_intensity)
  )

catchment_charater <- data_process %>%
  select(IA:LAB) %>%
  mutate(LAB = factor(LAB) %>% as.numeric(),
         Sum_CAP_MC = Sum_CAP_MC/area_new,
         area_new = log10(area_new))
  

out <- catchment_id %>%
  cbind(catchment_charater) %>%
  cbind(forcing) %>%
  cbind(flood_charater)


write_csv(out, file = "./data/annual_max_flood.csv")











data_process2 <- data_process %>%
  filter(LAB == "CNA")

ggplot(data_process, aes(sum_pcp, `runoff depth`, color = log10(area_new)))+
  geom_point()+
  scale_colour_gradientn(colours = terrain.colors(10))+
  facet_wrap(~LAB, scales = "free") +
  labs(x = "Total precipitation [mm]",
       y = "Runoff depth [mm]",
       color = "Log10 of catchment area")







data_recipe <- recipe(x = data_process) %>%
  update_role(rain_duration_time:flood_type, new_role = "predictor") %>%
  update_role('runoff depth', new_role = "outcome") %>%
  update_role(id, new_role = "id variable") %>%
  step_dummy(all_nominal_predictors(), one_hot = T) %>%
  prep()








# explore -----------------------------------------------------------------

data_process %>%
  filter()


data_process %>%
  filter() %>%
  summarise(across(c(initial_dis:Sum_CAP_MC), min, na.rm = T)) %>%
  view()



data_process <- data_process %>%
  filter(year >= 1985,
         `runoff depth` > 0) # IA is available after 1985, some runoff depth is 0

data_process %>%
  transmute(total_land = IA + cropland + forest + grassland + shrubland + tundra + barren_land) %>%
  pull(total_land) %>%
  table()


  
  
  
  
  
  
  
  
  


