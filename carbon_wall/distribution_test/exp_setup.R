suppressPackageStartupMessages({
  library(tidyverse)
})

yco <- 0.5
p_bar <- 10 #bar
t_c <- 45 #C
s_box <- 100 #A
n_moves <- 1000
n_equil <- 50000
n_prod <- 10000
reps <- 10

sf <- c("False", "True", "True", "True")
W <- c(0, 20, 30, 40)



data.frame(
  yco = yco,
  p_bar = p_bar,
  t_c = t_c,
  s_box = s_box,
  n_moves = as.integer(n_moves),
  n_equil = as.integer(n_equil),
  n_prod = as.integer(n_prod),
  sf = rep( sf, each = reps),
  W = rep( W, each = reps)
) |>
  mutate( exp = 1:n()) |>
  write_csv("carbon_wall/distribution_test/distribution_test_exp.csv")
