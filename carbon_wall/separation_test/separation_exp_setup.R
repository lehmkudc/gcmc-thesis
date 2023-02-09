suppressPackageStartupMessages({
  library(tidyverse)
})
{
yco <- round(seq(0.1, 0.9, by = 0.1),1)
p_bar <- c(25, 50, 100) #bar
t_c <- c(30,45,60) #C
s_box <- 100 #A
n_moves <- 1000
n_equil <- 500000
n_prod <- 100000
reps <- 3

sf <- 1
W <- c(20, 30, 40)
}

const_p_t <- data.frame(
  yco = rep( yco, each = reps ),
  p_bar = 100,
  t_c = 45,
  s_box = s_box,
  n_moves = as.integer(n_moves),
  n_equil = as.integer(n_equil),
  n_prod = as.integer(n_prod),
  sf = 1,
  W = rep( W, each = reps*length(yco))
)

const_w_t <- data.frame(
  yco = rep( yco, each = reps ),
  p_bar = rep(p_bar, each = reps*length(yco)),
  t_c = 45,
  s_box = s_box,
  n_moves = as.integer(n_moves),
  n_equil = as.integer(n_equil),
  n_prod = as.integer(n_prod),
  sf = 1,
  W = 30
)

const_p_w <- data.frame(
  yco = rep( yco, each = reps ),
  p_bar = 100,
  t_c = rep( t_c, each = reps*length(yco)),
  s_box = s_box,
  n_moves = as.integer(n_moves),
  n_equil = as.integer(n_equil),
  n_prod = as.integer(n_prod),
  sf = 1,
  W = 30
)

bind_rows(
  const_p_t,
  const_w_t,
  const_p_w
) |> 
  mutate(
    exp = 1:n()
  ) |> 
  filter( exp > 201 ) |>
  write_csv("./carbon_wall/separation_test/separation_exp_list2.csv")



