

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


combineRunData(
  paste0(homeDir, "data/"),
  paste0(homeDir, "aspen_comparison.csv")
)
