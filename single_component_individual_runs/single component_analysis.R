library( tidyverse )
library( data.table )
library( glue )



#### Collect the data ####
homeDir <-"./single_component_individual_runs/lower_p/"

combineRunData <- function( dataDir, expListFile ){
  
  dataFiles <- list.files(dataDir)
  expList <- fread(  expListFile )
  
  # all the raw data
  rawDataTable <- lapply( dataFiles, function(file){
    
    exp <- str_remove( file, pattern = ".csv")
    
    expInputs <- expList %>%
      filter( exp == !!exp )
    
    expOutputs <- fread( paste0( dataDir, file ) )
    
    
    expData <- bind_cols(expInputs, expOutputs) %>%
      mutate( iter = 1:n() )
    
    return( expData )
    
  }) %>% bind_rows
  
  return( rawDataTable )
}

co2_30 <- combineRunData(
  paste0(homeDir, "only_co2_30/"),
  paste0( homeDir, "only_co2_30.csv" )
)

co2_45 <- combineRunData(
  paste0(homeDir, "only_co2_45/"),
  paste0( homeDir, "only_co2_45.csv" )
)

co2_60 <- combineRunData(
  paste0(homeDir, "only_co2_60/"),
  paste0( homeDir, "only_co2_60.csv" )
)


nist_co2_30 <- fread('reference_data/c02_30.txt')
nist_co2_45 <- fread('reference_data/c02_45.txt')
nist_co2_60 <- fread('reference_data/c02_60.txt')


allRunDataCO2 <- bind_rows(
  co2_30, co2_45, co2_60
) %>%
  mutate(
    rhoco = rhocov*10**4/6.02
  )

allRunDataCO2 %>%
  glimpse()

allNistDataCO2 <- bind_rows(
  nist_co2_30, nist_co2_45, nist_co2_60
) %>%
  mutate(
    t_c = as.factor(`Temperature (C)`),
    p_bar = `Pressure (bar)`,
    rhoco = `Density (mol/l)`
  )


middle <- allRunDataCO2 %>%
  filter( p_bar == 100 ) %>%
  group_by( t_c, iter ) %>%
  summarise(
    Pvm = mean( Pv ),
    Envm = mean( Env ),
    rhocom = mean( rhoco )
  ) %>%
  ungroup() %>%
  mutate(
    t_c = as.factor(t_c)
  )

ggplot()+
  geom_smooth(
    data = middle,
    mapping = aes(x = iter, y = Pvm, col = t_c, group = t_c)
  )

ggplot()+ 
  geom_smooth(
    data = middle,
    mapping = aes(x = iter, y = rhocom, col = t_c, group = t_c)
  ) 

allNistDataCO2 %>%
  filter( p_bar == 100 )


at_60 <- allRunDataCO2 %>%
  filter( t_c == 60 ) %>%
  group_by( p_bar, iter ) %>%
  summarise(
    Pvm = mean( Pv ),
    Envm = mean( Env ),
    rhocom = mean( rhoco )
  ) %>%
  ungroup()

ggplot()+
  geom_smooth(
    data = at_60,
    mapping = aes(x = iter, y = Pvm, col = p_bar, group = p_bar)
  )

at_30 <- allRunDataCO2 %>%
  filter( t_c ==30 ) %>%
  group_by( p_bar, iter ) %>%
  summarise(
    Pvm = mean( Pv ),
    Envm = mean( Env ),
    rhocom = mean( rhoco )
  ) %>%
  ungroup()

ggplot()+
  geom_smooth(
    data = at_30,
    mapping = aes(x = iter, y = Pvm, col = p_bar, group = p_bar)
  )




# Iteration to call "equilibrated"
targetIter <- 250

# All reps Data at all pressures only result
resultData <- allRunDataCO2 %>%
  filter( iter > targetIter ) %>%
  group_by( exp, p_bar, t_c ) %>%
  summarise(
    Pvm = mean(Pv),
    Envm = mean(Env),
    rhocom = mean(rhoco),
    pdiff = Pvm - p_bar,
    pdiff_frac = pdiff/p_bar
  ) %>%
  ungroup() %>%
  mutate( t_c = as.factor( t_c ))

ggplot() + 
  geom_point(
    data = resultData,
    mapping = aes(x = p_bar, y = rhocom, color = t_c, group = t_c)
  ) + 
  geom_line(
    data = allNistDataCO2,
    aes( x = p_bar, y = rhoco, color = t_c, group = t_c)
  )





me_30 <- combineRunData(
  paste0(homeDir, "only_me_30/"),
  paste0( homeDir, "only_me_30.csv" )
)

me_45 <- combineRunData(
  paste0(homeDir, "only_me_45/"),
  paste0( homeDir, "only_me_45.csv" )
)

me_60 <- combineRunData(
  paste0(homeDir, "only_me_60/"),
  paste0( homeDir, "only_me_60.csv" )
)


nist_co2_30 <- fread('reference_data/c02_30.txt')
nist_co2_45 <- fread('reference_data/c02_45.txt')
nist_co2_60 <- fread('reference_data/c02_60.txt')


allRunDataCO2 <- bind_rows(
  co2_30, co2_45, co2_60
) %>%
  mutate(
    rhoco = rhocov*10**4/6.02
  )

allRunDataCO2 %>%
  glimpse()

allNistDataCO2 <- bind_rows(
  nist_co2_30, nist_co2_45, nist_co2_60
) %>%
  mutate(
    t_c = as.factor(`Temperature (C)`),
    p_bar = `Pressure (bar)`,
    rhoco = `Density (mol/l)`
  )


middle <- allRunDataCO2 %>%
  filter( p_bar == 100 ) %>%
  group_by( t_c, iter ) %>%
  summarise(
    Pvm = mean( Pv ),
    Envm = mean( Env ),
    rhocom = mean( rhoco )
  ) %>%
  ungroup() %>%
  mutate(
    t_c = as.factor(t_c)
  )

ggplot()+
  geom_smooth(
    data = middle,
    mapping = aes(x = iter, y = Pvm, col = t_c, group = t_c)
  )

ggplot()+ 
  geom_smooth(
    data = middle,
    mapping = aes(x = iter, y = rhocom, col = t_c, group = t_c)
  ) 

allNistDataCO2 %>%
  filter( p_bar == 100 )


at_60 <- allRunDataCO2 %>%
  filter( t_c == 60 ) %>%
  group_by( p_bar, iter ) %>%
  summarise(
    Pvm = mean( Pv ),
    Envm = mean( Env ),
    rhocom = mean( rhoco )
  ) %>%
  ungroup()

ggplot()+
  geom_smooth(
    data = at_60,
    mapping = aes(x = iter, y = Pvm, col = p_bar, group = p_bar)
  )

at_30 <- allRunDataCO2 %>%
  filter( t_c ==30 ) %>%
  group_by( p_bar, iter ) %>%
  summarise(
    Pvm = mean( Pv ),
    Envm = mean( Env ),
    rhocom = mean( rhoco )
  ) %>%
  ungroup()

ggplot()+
  geom_smooth(
    data = at_30,
    mapping = aes(x = iter, y = Pvm, col = p_bar, group = p_bar)
  )




# Iteration to call "equilibrated"
targetIter <- 250

# All reps Data at all pressures only result
resultData <- allRunDataCO2 %>%
  filter( iter > targetIter ) %>%
  group_by( exp, p_bar, t_c ) %>%
  summarise(
    Pvm = mean(Pv),
    Envm = mean(Env),
    rhocom = mean(rhoco),
    pdiff = Pvm - p_bar,
    pdiff_frac = pdiff/p_bar
  ) %>%
  ungroup() %>%
  mutate( t_c = as.factor( t_c ))

ggplot() + 
  geom_point(
    data = resultData,
    mapping = aes(x = p_bar, y = rhocom, color = t_c, group = t_c)
  ) + 
  geom_line(
    data = allNistDataCO2,
    aes( x = p_bar, y = rhoco, color = t_c, group = t_c)
  )









