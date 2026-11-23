# Runs an SDP optimisation for the bushfire example

#' @param mortality_w_fire Daily mortality in the wild if a fire occurs
#' @param mortality_w_fire_sd Variation in fire mortality
#' @param phi_wc Wild to captive translocation success rate (surviving transition + first month in captivity)
#' @param phi_wc_sd Variation in wild to captive translocation success rate
#' @param phi_cw Captive to wild translocation success rate (surviving transition + first month in the wild)
#' @param phi_cw_sd Variation in captive to wild translocation success rate 
#' @param cw_translocate_cost Cost of translocating 100 individuals captive->wild
#' @param wc_translocate_cost Cost of translocating 100 individuals wild->captive
#' @param captive_cost Cost of maintaining 100 individuals in captivity
#' @param establish_cost Cost of establishing a new captive population if none currently exists
#' @param cost_weight Proportion of weight to apply to cost versus to population
#' @param fire_transitions (4*4 matrix of fire transitions)
#' @param Tmax Time horizon in days (defaults to 60 days)
#' @param discount Discount rate (defaults to 0.99)
#' @param nsamples Number of draws from joint probability distribution
#' @param seed Seed for random draws of joint probability distribution
#' @returns Array of optimal decisions given current state.


run_bushfire_sdp <- function(
  mortality_w_fire, mortality_w_fire_sd, 
  phi_wc, phi_wc_sd, 
  phi_cw, phi_cw_sd, 
  cw_translocate_cost,
  wc_translocate_cost,
  captive_cost,
  establish_cost,
  fire_transitions,
  discount = 0.999,
  Tmax = 360,
  cost_weight,
  nsamples = 10000,
  seed = 1512
){






}