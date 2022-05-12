library(tidyverse)


kappa <- function(ohm){
  return(  0.37464 + 1.5422*ohm - 0.26992*ohm**2 )
}

alpha <- function(kappa, Tc, Tr){
  return( ( 1 + kappa*(1 - sqrt(Tr/Tc) ) )**2 )
}

a <- function( Tc, Pc, alpha, P_res, T ){
  R = .08314 #L*bar/mol/K
  return( 0.45724*R**2*Tc**2*alpha/Pc )
}

A <- function(a, P_res, Tr ){
  R = .08314 #L*bar/mol/K
  return(a*P_res/R**2/Tr**2)
}

b <- function( Tc, Pc ){
  R = .08314 #L*bar/mol/K
  return( 0.07780*R*Tc/Pc )
}

B <- function( b, P_res, Tr){
  R = .08314 #L*bar/mol/K
  return( b*P_res/R/Tr)
}

solveZ <- function( A, B ){
  c1 = B-1
  c2 = A - 3*B**2 - 2*B
  c3 = B**3 + B**2 - A*B
  
  roots = polyroot( c(c3, c2, c1, 1))
  #roots = polyroot( c(1, c1, c2, c3) )
  
  for (r in roots){
    if(  round( Im(r), 7) == 0 ){
      return( Re(r) )
    }
  }
}

Fi <- function( Aii, Ajj, Aij, A_mix, Bi, B_mix, Z_mix, yi, P_res ){
  
  term1 = Bi*(Z_mix - 1)/B_mix
  term2 = log( Z_mix - B_mix )
  term3 = A_mix/( 2*sqrt(2)*B_mix )
  term4 = 2*( yi*Aii + (1-yi)*Aij )/A_mix - Bi/B_mix
  term5 = (Z_mix + (1 + sqrt(2))*B_mix )/(Z_mix + (1 - sqrt(2))*B_mix )
  term6 = term1 - term2 - term3*term4*log( term5 )
  return(yi*P_res*np.exp( term6 ) )
  
}

Fone <- function( Z, A, B, P_res){
  term1 = Z - 1 - log( Z - B)
  term2 = A/( 2*sqrt(2)*B )
  term3 = (Z + (1 + sqrt(2))*B )/(Z + (1 - sqrt(2))*B )
  term4 = term1 - term2*log( term3 )
  return( exp( term4 )*P_res )
}

PR_Fugacity <- function(P_res, Tr, yco){
  Tc_co = 304.2 #K
  Tc_me = 190.6 #K
  Pc_co = 73.76 #bar
  Pc_me = 46 #bar
  ohm_co = 0.225
  ohm_me = 0.008
  
  kappa_me = kappa( ohm_me )
  alpha_me = alpha( kappa_me, Tc_me, Tr)
  a_me = a( Tc_me, Pc_me, alpha_me, P_res, Tr )
  A_me = A( a_me, P_res, Tr )
  b_me = b( Tc_me, Pc_me )
  B_me = B( b_me, P_res, Tr)
  
  kappa_co = kappa( ohm_co )
  alpha_co = alpha( kappa_co, Tc_co, Tr)
  a_co = a( Tc_co, Pc_co, alpha_co, P_res, Tr )
  A_co = A( a_co, P_res, Tr )
  b_co = b( Tc_co, Pc_co )
  B_co = B( b_co, P_res, Tr)
  
  a_cm = sqrt(a_me*a_co )*(1-0.0919)
  A_cm = A( a_cm, P_res, Tr)
  
  a_mix = yco*yco*a_co + 2*yco*(1-yco)*a_cm + (1-yco)*(1-yco)*a_me
  b_mix = yco*b_co + (1-yco)*b_me
  
  A_mix = A( a_mix, P_res, Tr )
  B_mix = B( b_mix, P_res, Tr )
  Z_mix = solveZ( A_mix, B_mix )
  
  F_me = Fi( A_me, A_co, A_cm, A_mix, B_me, B_mix, Z_mix, 1-yco, P_res)
  F_co = Fi( A_co, A_me, A_cm, A_mix, B_co, B_mix, Z_mix, yco, P_res)
  
  return( c(F_co, F_me) )

}

PR_Fugacity_Single <- function( P_res, Tr, species ){
  Tc_co = 304.2 #K
  Tc_me = 190.6 #K
  Pc_co = 73.76 #bar
  Pc_me = 46 #bar
  ohm_co = 0.225
  ohm_me = 0.008
  
  if (species == "co2"){
    kappa_co = kappa( ohm_co )
    alpha_co = alpha( kappa_co, Tc_co, Tr)
    a_co = a( Tc_co, Pc_co, alpha_co, P_res, Tr )
    A_co = A( a_co, P_res, Tr )
    b_co = b( Tc_co, Pc_co )
    B_co = B( b_co, P_res, Tr)
    Z_co = solveZ( A_co, B_co )
    return( Fone( Z_co, A_co, B_co, P_res ) )
  } else if (species == "me"){
    kappa_me = kappa( ohm_me )
    alpha_me = alpha( kappa_me, Tc_me, Tr)
    a_me = a( Tc_me, Pc_me, alpha_me, P_res, Tr )
    A_me = A( a_me, P_res, Tr )
    b_me = b( Tc_me, Pc_me )
    B_me = B( b_me, P_res, Tr)
    Z_me = solveZ( A_me, B_me )
    return( Fone( Z_me, A_me, B_me, P_res ) )
  } else {
    stop("Unknown Species")
  }

}

PR_Zmix <- function(P_res, Tr, yco){
  
  Tc_co = 304.2 #K
  Tc_me = 190.6 #K
  Pc_co = 73.76 #bar
  Pc_me = 46 #bar
  ohm_co = 0.225
  ohm_me = 0.008
  
  kappa_me = kappa( ohm_me )
  alpha_me = alpha( kappa_me, Tc_me, Tr)
  a_me = a( Tc_me, Pc_me, alpha_me, P_res, Tr )
  A_me = A( a_me, P_res, Tr )
  b_me = b( Tc_me, Pc_me )
  B_me = B( b_me, P_res, Tr)
  
  kappa_co = kappa( ohm_co )
  alpha_co = alpha( kappa_co, Tc_co, Tr)
  a_co = a( Tc_co, Pc_co, alpha_co, P_res, Tr )
  A_co = A( a_co, P_res, Tr )
  b_co = b( Tc_co, Pc_co )
  B_co = B( b_co, P_res, Tr)
  
  a_cm = sqrt(a_me*a_co )*(1-0.0919)
  A_cm = A( a_cm, P_res, Tr)
  
  a_mix = yco*yco*a_co + 2*yco*(1-yco)*a_cm + (1-yco)*(1-yco)*a_me
  b_mix = yco*b_co + (1-yco)*b_me
  
  A_mix = A( a_mix, P_res, Tr )
  B_mix = B( b_mix, P_res, Tr )
  Z_mix = solveZ( A_mix, B_mix )
  
  return( Z_mix )
}
