library( tidyverse )
library( data.table )
library(glue)

#### Collect the data ####
homeDir <-"./single_component_individual_runs/lower_p/"
dataDir <- paste0( homeDir, "only_co2_45/")
expList <- fread( paste0( homeDir, "only_co2_45.csv"))

dataFiles <- list.files(dataDir)

# all the raw data
rawDataTable <- lapply( dataFiles, function(file){
  
  exp <- str_remove( file, pattern = ".csv")
  
  expInputs <- expList %>%
    filter( exp == !!exp )
  
  expOutputs <- fread( paste0( dataDir, file ) )
  
  
  expData <- bind_cols(expInputs, expOutputs) %>%
    mutate( iter = 1:n() )
  
  return( expData )
  
}) %>% bind_rows()


nistData <- fread( "reference_data/c02.txt" ) %>%
  mutate(
    p_bar = `Pressure (bar)`,
    rhocov = `Density (mol/l)`*6.02*10**(-4)
  )


#### Process Data ####

const <- expList %>% distinct(yco, t_c, s_box, n_moves, n_equil, n_prod )

# Single Pressure to Zoom on
targetP <- 100

# all data at pressure
targetRawData <- rawDataTable %>%
  filter( p_bar == targetP )

# summarized data at pressure over iterations
targetMeanRunData <- rawDataTable %>%
  filter( p_bar == targetP ) %>%
  group_by( iter, p_bar ) %>%
  summarize(
    Pvm = mean( Pv ),
    Envm = mean(Env),
    rhocovm = mean(rhocov),
    pdiff = Pvm - p_bar
  )


# summarized data at all pressures over iterations
meanRunData <- rawDataTable %>%
  group_by( iter, p_bar ) %>%
  summarise( 
    Pvm = mean(Pv),
    Envm = mean(Env),
    rhocovm = mean(rhocov),
    pdiff = Pvm - p_bar
  ) %>%
  ungroup()

# Iteration to call "equilibrated"
targetIter <- 250

# All reps Data at all pressures only result
resultData <- rawDataTable %>%
  filter( iter > targetIter ) %>%
  group_by( exp, p_bar ) %>%
  summarise(
    Pvm = mean(Pv),
    Envm = mean(Env),
    rhocovm = mean(rhocov),
    pdiff = Pvm - p_bar,
    pdiff_frac = pdiff/p_bar
  )

# Summarized Data at all pressures only result
meanResultData <- rawDataTable %>%
  filter( iter > targetIter ) %>%
  group_by( p_bar ) %>%
  summarise(
    Pvm = mean(Pv),
    Envm = mean(Env),
    rhocovm = mean(rhocov),
    pdiff = Pvm - p_bar,
    pdiff_frac = pdiff/p_bar
  )


#### Plot Pressures ####

ggplot() +
  stat_smooth(
    data = targetRawData,
    aes(x = iter, y = Pv ),
    alpha = 1,
    color = "blue"
  ) + 
  geom_point(
    data = targetRawData,
    aes(x = iter, y = Pv ),
    alpha = 0.5,
    color = "black"
  ) +
  # geom_line(
  #   data = targetMeanRunData,
  #   aes( x = iter, y = Pvm ),
  #   color = "red",
  #   size = 1
  # ) + 
  geom_point(
    data = targetMeanRunData,
    aes( x = iter, y = p_bar ),
    color = "red",
    size = 2
  ) +
  ylab( "Pressure of Simulation [MPa]") + 
  xlab( glue("Cycles [{const$n_moves} steps]") ) +
  labs(
    title = glue("Simulation Pressure of Pure CO2 @ {targetP} MPa Over Simulation"),
    subtitle = glue("Pure CO2 @ {const$t_c} C"),
    caption = glue("s: {const$s_box} A | Ntotal: {const$n_prod}")
  )
  

ggplot() + 
  geom_smooth(
    data = meanRunData,
    aes( x = iter, y = Pvm, group = p_bar, col = p_bar),
  ) +
  ylab( "Pressure of Simulation [MPa]") + 
  xlab( glue("Cycles [{const$n_moves} steps]") ) +
  labs(
    title = glue("Simulation Pressure of Pure CO2 @ Various P Over Simulation"),
    subtitle = glue("Pure CO2 @ {const$t_c} C"),
    caption = glue("s: {const$s_box} A | Ntotal: {const$n_prod}"),
    col = "Res. P [MPa]"
  )


ggplot() + 
  stat_smooth(
    data = meanRunData,
    aes( x = iter, y = Envm, group = p_bar, col = p_bar)
  ) +
  ylab( "Internal Energy [J/K]") + 
  xlab( glue("Cycles [{const$n_moves} steps]") ) +
  labs(
    title = glue("Internal Energy of Pure CO2 @ Various P Over Simulation"),
    subtitle = glue("Pure CO2 @ {const$t_c} C"),
    caption = glue("s: {const$s_box} A | Ntotal: {const$n_prod}"),
    col = "Res. P [MPa]"
  )


ggplot() + 
  stat_smooth(
    data = meanRunData,
    aes( x = iter, y = rhocovm, group = p_bar, col = p_bar)
  ) +
  ylab( "Density [A^-3]") + 
  xlab( glue("Cycles [{const$n_moves} steps]") ) +
  labs(
    title = glue("Density of Pure CO2 @ Various P Over Simulation"),
    subtitle = glue("Pure CO2 @ {const$t_c} C"),
    caption = glue("s: {const$s_box} A | Ntotal: {const$n_prod}"),
    col = "Res. P [MPa]"
  )


ggplot() + 
  stat_smooth(
    data = meanRunData,
    aes( x = iter, y = pdiff, group = p_bar, col = p_bar)
  ) +
  xlim( targetIter , 1000)



ggplot() + 
  geom_point(
    data = resultData,
    aes( x = p_bar, y = Pvm )
  ) +
  geom_line(
    data = resultData,
    aes( x = p_bar, y = p_bar )
  ) +
  stat_smooth(
    data = resultData,
    aes( x = p_bar, y = Pvm)
  ) +
  ylab( "Pressure of Simulation [MPa]") + 
  xlab( glue("Pressure of Reservoir [MPa]") ) +
  labs(
    title = glue("Pressure Consistency of Pure CO2 @ Various P"),
    subtitle = glue("Pure CO2 @ {const$t_c} C"),
    caption = glue("s: {const$s_box} A | Ntotal: {const$n_prod}")
  )

ggplot() + 
  geom_point(
    data = resultData,
    aes( x = p_mpa, pdiff )
  ) +
  geom_line(
    data = resultData,
    aes( x = p_mpa, y = 0)
  )   +
  ylab( "Pressure Difference [MPa]") + 
  xlab( glue("Pressure of Reservoir [MPa]") ) +
  labs(
    title = glue("Pressure Consistency of Pure CO2 @ Various P"),
    subtitle = glue("Pure CO2 @ {const$t_c} C"),
    caption = glue("s: {const$s_box} A | Ntotal: {const$n_prod}")
  )


ggplot() + 
  geom_point(
    data = resultData,
    aes( x = p_mpa, pdiff_frac )
  ) +
  geom_line(
    data = resultData,
    aes( x = p_mpa, y = 0)
  )  +
  stat_smooth(
    data = resultData,
    aes( x = p_mpa, pdiff_frac)
  ) +
  ylab( "Pressure Difference Fraction [MPa]") + 
  xlab( glue("Pressure of Reservoir [MPa]") ) +
  labs(
    title = glue("Pressure Consistency of Pure CO2 @ Various P"),
    subtitle = glue("Pure CO2 @ {const$t_c} C"),
    caption = glue("s: {const$s_box} A | Ntotal: {const$n_prod}")
  )



ggplot() + 
  geom_point(
    data = resultData,
    aes( x = p_mpa, Envm )
  ) +
  ylab( "Internal Energy [J/K]") + 
  xlab( glue("Pressure of Reservoir [MPa]") ) +
  labs(
    title = glue("Internal Energy of Pure CO2 @ Various P"),
    subtitle = glue("Pure CO2 @ {const$t_c} C"),
    caption = glue("s: {const$s_box} A | Ntotal: {const$n_prod}")
  )


ggplot() + 
  geom_point(
    data = resultData,
    aes( x = p_bar, rhocovm )
  ) +
  ylab( "Density [A^-3]") + 
  xlab( glue("Pressure of Reservoir [MPa]") ) +
  labs(
    title = glue("Density of Pure CO2 @ Various P"),
    subtitle = glue("Pure CO2 @ {const$t_c} C"),
    caption = glue("s: {const$s_box} A | Ntotal: {const$n_prod}")
  ) +
  geom_line(
    data = nistData,
    aes( x = p_bar, y = rhocov)
  )



