library( tidyverse )
library( data.table )
library(glue)



allFiles <- list.files("reference_data/")

aspenFiles <- allFiles[ grepl( pattern = "CO2-CH4", allFiles)]

aspenData <- lapply( aspenFiles, function( filename ){
  
  dataChunk <- fread( paste0("reference_data/", filename), skip = 1)
  
  return( dataChunk )
  
}) %>%
  bind_rows()


aspenData <- aspenData  %>% 
  mutate(
    yco = round( `Mole fraction CO2`, 1),
    t_c = round( `Temperature (K)` - 273.15, 1 ),
    p_bar = round(`Pressure (N/sqm)`*10**(-5), 1),
    rhoco = yco*`Density (kg/cum)`/44.01,
    rhome = (1-yco)*`Density (kg/cum)`/16.04
  )



aspenData %>%
  group_by( t_c ) %>%
  summarise(
    count = n(),
    p_min = min(p_bar),
    p_max = max(p_bar)
  ) %>%
  arrange( desc(count))

aspenData %>%
  mutate(
    yco = `Mole fraction CO2`,
    t_c = round( `Temperature (K)` - 273.15, 1 ),
    p_bar = `Pressure (N/sqm)`*10**(-5)
  )


aspenData %>%
  distinct(t_c, p_bar, yco) %>%
  mutate(
    exp = 1:n(),
    s_box = 34,
    n_moves = 1000,
    n_equil = as.integer(1000000),
    n_prod =  as.integer( 100000)
  ) %>%
  select( exp, everything() ) %>%
  write.csv("aspen_comparison.csv",row.names = FALSE)
