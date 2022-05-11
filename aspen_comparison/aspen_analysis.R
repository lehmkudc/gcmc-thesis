

library( tidyverse )
library( data.table )
library(glue)

source( "reference_data/pr.R")


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
    rhome = (1-yco)*`Density (kg/cum)`/16.04,
    source = "Aspen"
  )


#### Collect the data ####
homeDir <-"./aspen_comparison/"

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


aspenRunData <- combineRunData(
  paste0(homeDir, "data/"),
  paste0(homeDir, "aspen_comparison.csv")
) %>%
  mutate(
    rhocov = rhocov*10000/6.02,
    rhomev = rhomev*10000/6.02
  )

comparisonData <- aspenData %>%
  #filter( t_c < 80) %>%
  #filter( p_bar < 250 ) %>%
  select( yco, t_c, p_bar, rhoco, rhome) %>%
  left_join(
    aspenRunData %>%
      select( yco, t_c, p_bar, Pv, Env, rhocov, rhomev ),
    by = c("yco", "t_c", "p_bar")
  )

comparisonData %>%
  mutate(
    z_pr = PR_Zmix(p_bar, t_c + 273.15, yco ),
    rho_pr = p/z/R/(t+273.15),
    rhoco_pr = yco*rho_pr,
    rhome_pr = 1
  )

ggplot( comparisonData ) +
  geom_point(
    aes( x = yco, y = p_bar )
  )



ggplot( comparisonData ) +
  geom_point(
    aes( x = yco, y = rhocov/(rhocov + rhomev) ),
    alpha = 0.2
  )




ggplot( comparisonData ) +
  geom_point(
    aes( x = p_bar, y = Pv, col = yco ),
    alpha = 1
  ) +
  geom_line(
    aes( x = p_bar, y = p_bar)
  )



ggplot( comparisonData ) +
  geom_point(
    aes( x = p_bar, y = Pv, col = t_c ),
    alpha = 1
  ) +
  geom_line(
    aes( x = p_bar, y = p_bar)
  )


ggplot( comparisonData ) +
  geom_point(
    aes( x = rhoco, y = rhocov, col = t_c ),
    alpha = 1
  ) +
  geom_line(
    aes( x = rhoco, y = rhoco)
  )


ggplot( comparisonData ) +
  geom_point(
    aes( x = rhoco, y = rhocov, col = p_bar ),
    alpha = 1
  ) +
  geom_line(
    aes( x = rhoco, y = rhoco)
  )


ggplot( comparisonData ) +
  geom_point(
    aes( x = rhoco, y = rhocov, col = yco ),
    alpha = 1
  ) +
  geom_line(
    aes( x = rhoco, y = rhoco)
  )





ggplot( comparisonData ) +
  geom_point(
    aes( x = rhome, y = rhomev, col = t_c ),
    alpha = 1
  ) +
  geom_line(
    aes( x = rhome, y = rhome)
  )


ggplot( comparisonData ) +
  geom_point(
    aes( x = rhome, y = rhomev, col = p_bar ),
    alpha = 1
  ) +
  geom_line(
    aes( x = rhome, y = rhome)
  )



ggplot( comparisonData ) +
  geom_point(
    aes( x = rhome, y = rhomev, col = yco ),
    alpha = 1
  ) +
  geom_line(
    aes( x = rhome, y = rhome)
  )

