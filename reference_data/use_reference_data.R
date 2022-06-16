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
  write.csv("./aspen_comparison/aspen_comparison.csv",row.names = FALSE)





table_pieces <- readLines("reference_data/mondejar_et_al_raw.txt",encoding="Latin-1") %>%
  strsplit( "SPACE") %>% unlist() %>% map( function(x){
    table_piece <- strsplit( x, " ") %>% 
      unlist() %>%
      .[. != "ï»¿"] %>%
      .[. != ""] %>%  
      str_replace(string = ., pattern = "âˆ’", replacement = "-") %>% 
      as.numeric() %>%
      matrix(., ncol = 5,byrow = TRUE) %>%
      as.data.frame()
    
    colnames( table_piece ) <- c(
      "Tk", "P", "density", "rho_diff", "uncertainty"
    )
    
    table_piece <- table_piece %>%
      select( Tk, P, density )
    
    return(
      table_piece 
    )
  })

bind_rows(
  table_pieces[[1]] %>% mutate(yco = 0.2),
  table_pieces[[2]] %>% mutate(yco = 0.2),
  table_pieces[[3]] %>% mutate(yco = 0.2),
  table_pieces[[4]] %>% mutate(yco = 0.2),
  table_pieces[[5]] %>% mutate(yco = 0.4),
  table_pieces[[6]] %>% mutate(yco = 0.4),
  table_pieces[[7]] %>% mutate(yco = 0.6),
  table_pieces[[8]] %>% mutate(yco = 0.6),
) %>%
  mutate(
    t_c = Tk - 273.15, #[*C]
    p_bar = P*10, #[bar]
    mfrac_co = ( 44.01*yco )/( 44.01*yco + 16.04*(1-yco)), 
    rhoco = mfrac_co*density/44.01, #[mol/L]
    rhome = (1-mfrac_co)*density/16.04
  )



bind_rows(
  table_pieces[[1]] %>% mutate(yco = 0.2),
  table_pieces[[2]] %>% mutate(yco = 0.2),
  table_pieces[[3]] %>% mutate(yco = 0.2),
  table_pieces[[4]] %>% mutate(yco = 0.2),
  table_pieces[[5]] %>% mutate(yco = 0.4),
  table_pieces[[6]] %>% mutate(yco = 0.4),
  table_pieces[[7]] %>% mutate(yco = 0.6),
  table_pieces[[8]] %>% mutate(yco = 0.6),
) %>% 
  mutate( 
    exp = 1:n(),
    t_c = Tk - 273.15,
    p_bar = P*10,
    s_box = 34,
    n_moves = 1000,
    n_equil = as.integer(1000000),
    n_prod  = as.integer(100000)
  ) %>%
  select( exp, yco, p_bar, t_c, s_box, n_moves, n_equil, n_prod ) %>%
  write_csv( "mixture_grid/mondejar_et_al.csv")






