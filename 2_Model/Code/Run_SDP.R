# Function to run stochastic dynamic programming for one species ----------
#'
#' @param fecundity_w Fecundity of wild population
#' @param fecundity_c Fecundity of captive population
#' @param breeding_month Month where species breeds
#' @param mortality_w_base Monthly mortality in the wild (uninfected)
#' @param mortality_w_flu Monthly mortality in the wild (infected)
#' @param mortality_c Monthly mortality in captivity
#' @param phi_wc Wild to captive translocation success rate (surviving transition + first month in captivity)
#' @param phi_cw Captive to wild translocation success rate (surviving transition + first month in the wild)
#' @param disease_transitions (6*6 matrix of disease transitions)
#' @param discount Discount rate (defaults to 0.99)
#' @param Tmax Time horizon in months (defaults to 30 years)
#' @param start_month First month (as number)
#' @param cw_translocate_cost Cost of translocating 100 individuals captive->wild
#' @param wc_translocate_cost Cost of translocating 100 individuals wild->captive
#' @param captive_cost Cost of maintaining 100 individuals in captivity
#' @param establish_cost Cost of establishing a new captive population if none currently exists
#' @param cost_weight Proportion of weight to apply to cost versus to population
#'
#' @returns Array of optimal decisions given current state.

run_SDP <- function(fecundity_w, fecundity_c, breeding_month,
                    mortality_w_base, mortality_w_flu, mortality_c,
                    phi_wc, phi_cw,
                    disease_transitions,
                    discount = 0.999,
                    Tmax = 360,
                    start_month,
                    cw_translocate_cost,
                    wc_translocate_cost,
                    captive_cost,
                    establish_cost,
                    cost_weight){

  # Define states -----------------------------------------------------------
  # Define the states in the model. There are four dimensions to the state space:
  # Wild population, captive pop, H5N1 exposure, and month.
  
  # Wild pop: 0-2000, increments of 50
  Nstates_nWild <- 41
  states_nWild <- 1:41
  
  # Captive pop: 0-1000, increments of 50
  Nstates_nCaptive <- 21
  states_nCaptive <- 1:21
  
  # Infected: Exposure level; 6 is 'infected'
  Nstates_exposure <- 6
  states_exposure <- 1:6
  
  # Month: 1:12
  Nstates_month <- 12
  states_month <- 1:12
  
  # Flag non-breeding months
  nonbreeding_months <- seq(1,12)[1:12 %notin% breeding_month]
  
  ## Also include binary states for breeding/nonbreeding and infected/noninfected 
  # Breeding: Binary
  Nstates_breeding <- 2
  states_breeding <- c(FALSE, TRUE)
  
  # Infected: binary
  Nstates_infected <- 2
  states_infected <- c(FALSE, TRUE)
  
  ## Count number of states
  # Number of states (only population) 
  nStates_pop <- Nstates_nWild*Nstates_nCaptive
  
  # Total number of states
  nStates_full <- Nstates_nWild*Nstates_nCaptive*Nstates_exposure*Nstates_month


  # Define possible actions -------------------------------------------------
  # Positive actions are translocation captive to wild, negative actions are 
  # translocation wild to captive
  actions <- seq(-100, 100, by = 5)
  
  # Total number of actions
  nActions <- length(actions)
  

  # Define population model -------------------------------------------------
  # This function determines the next population of the wild and captive 
  # population given an action.
  dynamic <- function(actualWild, actualCaptive, action,
                      disease, breeding, 
                      mortality_w_base, mortality_w_base_sd = 0.01, 
                      mortality_w_flu, mortality_w_flu_sd = 0.1, 
                      mortality_c, mortality_c_sd = 0.01, 
                      fecundity_w, fecundity_w_sd = 0.1, 
                      fecundity_c, fecundity_c_sd = 0.1, 
                      phi_wc, phi_wc_sd = 0.01, 
                      phi_cw, phi_cw_sd = 0.01, 
                      nsamples = 10000){
    
    nextWild <- rep(NA, nsamples)
    nextCaptive <- rep(NA, nsamples)
    
    tau_cw <- ifelse(action > 0, action, 0)/100
    tau_wc <- ifelse(action < 0, abs(action), 0)/100
    
    
    mortality_w_flu_sample <- truncnorm::rtruncnorm(nsamples, 0, 1, 
                                                      mean = mortality_w_flu,
                                                      sd = mortality_w_flu_sd)
    mortality_w_base_sample <- truncnorm::rtruncnorm(nsamples, 0, 1, 
                                                     mean = mortality_w_base,
                                                     sd = mortality_w_base_sd)
    mortality_c_sample <- truncnorm::rtruncnorm(nsamples, 0, 1, 
                                                mean = mortality_c,
                                                sd = mortality_c_sd)
    fecundity_c_sample <- rnorm(nsamples,
                                mean = fecundity_c, 
                                sd = fecundity_c_sd)
    fecundity_w_sample <- rnorm(nsamples,
                                mean = fecundity_w, 
                                sd = fecundity_w_sd)
    phi_cw_sample <- truncnorm::rtruncnorm(nsamples, 0, 1, 
                                           mean = phi_cw,
                                           sd = phi_cw_sd)
    phi_wc_sample <- truncnorm::rtruncnorm(nsamples, 0, 1, 
                                           mean = phi_wc,
                                           sd = phi_wc_sd)

    mortality_w_sample <- ifelse(disease == 1, mortality_w_flu_sample, mortality_w_base_sample)
    
    nextWild <- 
      ((actualWild-(actualWild*tau_wc))*(1-mortality_w_sample)) + 
      ((actualWild-(actualWild*tau_wc))*(1-mortality_w_sample)*breeding*fecundity_w_sample) +
      (actualCaptive*tau_cw*phi_cw_sample*(1-mortality_c_sample))
    
    nextCaptive <- 
      ((actualCaptive-(actualCaptive*tau_cw))*(1-mortality_c_sample)) + 
      ((actualCaptive-(actualCaptive*tau_cw))*(1-mortality_c_sample)*breeding*fecundity_c_sample) +
      (actualWild*tau_wc*phi_wc_sample*(1-mortality_w_sample))
    
    nextWild    <- pmin(pmax(round(nextWild    / 50), 0), 40)
    nextCaptive <- pmin(pmax(round(nextCaptive / 50), 0), 20)
    
    # Get joint distribution of nextwild and nextcaptive
    pair_idx <- (nextWild + 1) + 41 * nextCaptive
    matrix(tabulate(pair_idx, nbins = 41 * 21), nrow = 41, ncol = 21) / nsamples

  }
  

  # Define utility function -------------------------------------------------
  # Define the value of an action given the next wild/captive populations. Here,
  # wild individuals are worth 2x captive individuals.
  state_util <- matrix(0,
                       nrow = 41,
                       ncol = 21)
  for(w in 1:41){
    for(x in 1:21){
      state_util[w,x] <- ((w-1)*2)+((x-1))
    }
  }
  
  get_utility <- function(state_util,
                          lastWildPop, 
                          lastCaptivePop, 
                          transitions_pop,
                          infected,
                          breeding,
                          action) {
    

    maxPopUtility <- max(state_util)
    
    
    popUtility <- state_util*transitions_pop[lastWildPop,
                                             lastCaptivePop,
                                             ,
                                             ,
                                             infected+1,
                                             breeding+1,
                                             which(actions==action)]
    
    maxCost <- max(
      c((establish_cost + (40*wc_translocate_cost) + (20*captive_cost)),
        (20*captive_cost),
        (20*cw_translocate_cost))
      )
    
    tau_cw <- ifelse(action > 0, action, 0)/100
    tau_wc <- ifelse(action < 0, abs(action), 0)/100
    
    if (lastCaptivePop == 1 & tau_wc > 0){
      setupCost <- establish_cost
    } else {
      setupCost <- 0
    }
    
  
    costUtility <- (((lastCaptivePop-1)-((lastCaptivePop-1)*(tau_cw)) +
                       ((lastWildPop-1)*tau_wc))*captive_cost) +
      (((lastCaptivePop-1)*tau_cw)*cw_translocate_cost) +
      (((lastWildPop-1)*tau_wc)*wc_translocate_cost) +
      (setupCost)
    
    totalUtility <- ((sum(popUtility)/maxPopUtility)*(1-cost_weight)) - 
      ((costUtility/maxCost)*cost_weight)
    
    return(totalUtility)
  }

  
  # Fill population transitions array --------------------------------------
  # Population transitions: Prob of going from one state to another. Determines
  # only future population, not month/exposure states.
  transitions_pop <- array(0, 
                           dim = c(Nstates_nWild, # W pop at t
                                   Nstates_nCaptive, # C pop at t
                                   Nstates_nWild, # W pop at t+1
                                   Nstates_nCaptive, # C pop at t+1
                                   2, # uninfected/infected at t
                                   2, # nonbreeding/breeding at t
                                   nActions)) # action
  
  # Utility array: Utility of an action in a given state
  utility <- array(0, 
                   dim = c(Nstates_nWild,
                           Nstates_nCaptive,
                           Nstates_exposure,
                           Nstates_month,
                           nActions))
  
  # Iterate over all possible population states
  for (w in 1:Nstates_nWild){
    for (x in 1:Nstates_nCaptive){
      for (y in 1:Nstates_infected){
        for (z in 1:Nstates_breeding){
          # And possible actions
          for (i in 1:nActions){
            transitions_pop[w, x, , , y, z, i] <- dynamic(actualWild = (states_nWild[w]-1)*50, 
                                                          actualCaptive = (states_nCaptive[x]-1)*50,
                                                          disease = states_infected[y],
                                                          action = actions[i],
                                                          mortality_w_base = mortality_w_base, 
                                                          mortality_w_flu = mortality_w_flu, 
                                                          mortality_c = mortality_c, 
                                                          breeding = states_breeding[z], 
                                                          fecundity_w = fecundity_w, 
                                                          fecundity_c = fecundity_c, 
                                                          phi_wc = phi_wc, 
                                                          phi_cw = phi_cw)
          }
        }
      }
    }
  }
  
  # Fill utility matrix
  for (w in 1:Nstates_nWild){
    for (x in 1:Nstates_nCaptive){
      for (y in 1:Nstates_infected){
        for (z in 1:Nstates_breeding){
          # And possible actions
          for (i in 1:nActions){
            # Compute utility of action - utility uses exposure/month states rather
            # than infected/breeding states, hence the below logic to fill it fully.
            if (y == 1 & z == 1){
              utility[w, x, 1:5, nonbreeding_months, i] <- get_utility(w,
                                                                       x,
                                                                       transitions_pop,
                                                                       FALSE,
                                                                       FALSE,
                                                                       actions[i])
            }
            if (y == 1 & z == 2){
              utility[w, x, 1:5, breeding_month, i] <- get_utility(w,
                                                                   x,
                                                                   transitions_pop,
                                                                   FALSE,
                                                                   TRUE,
                                                                   actions[i])
            }
            if (y == 2 & z == 1){
              utility[w, x, 6, nonbreeding_months, i] <- get_utility(w,
                                                                     x,
                                                                     transitions_pop,
                                                                     TRUE,
                                                                     FALSE,
                                                                     actions[i])
            }
            if (y == 2 & z == 2){
              utility[w, x, 6, breeding_month, i] <- get_utility(w,
                                                                 x,
                                                                 transitions_pop,
                                                                 TRUE,
                                                                 TRUE,
                                                                 actions[i])
            }
          }
        }
      }
    }
  }
  
          

  # Assign months to each timestep based on breeding month ------------------
  Tmonths <- rep_len(c(start_month:12,1:(start_month-1)),
                     length.out = Tmax)
  

  # Reshape arrays for matrix algebra ---------------------------------------
  # Number of population states
  nPop <- Nstates_nwild*Nstates_nCaptive 
  # Number of state/action combinations
  nRow <- nPop*nActions

  # Transitions: one (state x action) by next-state matrix per
  # infected (1/2) x breeding (1/2) combination
  PT <- vector("list", 4)
  dim(PT) <- c(2, 2)
  for (inf in 1:2) {
    for (br in 1:2) {
      slice <- transitions_pop[, , , , inf, br, ]   # 41 x 21 x 41 x 21 x nActions
      dim(slice) <- c(nPop, nPop, nActions)        # s x s' x i
      PT[[inf, br]] <- matrix(aperm(slice, c(1, 3, 2)), nrow = nRow)  # (s,i) x s'
    }
  }
  rm(slice)
  
  # Utility: one (state x action) by exposure matrix per month
  U_month <- lapply(1:Nstates_month, function(m) {
    u <- utility[, , , m, ]                        # 41 x 21 x 6 x nActions
    dim(u) <- c(nPop, Nstates_exposure, nActions)
    matrix(aperm(u, c(1, 3, 2)), nrow = nRow)      # (s,i) x y
  })
  
  # Tie-break order: smallest |action| first, then the positive one.
  # max.col(ties.method = "first") then reproduces the original rule.
  pref <- order(abs(actions), -actions)
  
  # Solve Bellman via backwards iteration -----------------------------------
  # Terminal value: (w-1)*2 + (x-1), identical across exposure states
  pop_value <- as.vector(outer((states_nWild - 1) * 2, states_nCaptive - 1, "+"))
  Vnext <- matrix(pop_value, nrow = nPop, ncol = Nstates_exposure)
  
  D <- array(0, dim = c(Nstates_nWild, Nstates_nCaptive,
                        Nstates_exposure, Nstates_month))
  Vt          <- matrix(0, nPop, Nstates_exposure)
  best_action <- matrix(0, nPop, Nstates_exposure)
  
  for (t in (Tmax - 1):1) {
    m  <- Tmonths[t]
    br <- (m == breeding_month) + 1
    
    # Disease step: expected next value for each current exposure y
    # M[s', y] = sum_z disease_transitions[y, z] * Vnext[s', z]
    M <- Vnext %*% t(disease_transitions)          # nPop x 6
    
    # Population step: all states and actions in two matrix products
    Q <- U_month[[m]]
    Q[, 1:5] <- Q[, 1:5] + discount * (PT[[1, br]] %*% M[, 1:5])  # not infected
    Q[, 6]   <- Q[, 6]   + discount * (PT[[2, br]] %*% M[, 6])    # infected
    
    # Best action and its value for each (state, exposure)
    dim(Q) <- c(nPop, nActions, Nstates_exposure)
    for (y in 1:Nstates_exposure) {
      Qy  <- Q[, pref, y]                          # columns in tie-break order
      col <- max.col(Qy, ties.method = "first")
      Vt[, y]          <- Qy[cbind(seq_len(nPop), col)]
      best_action[, y] <- actions[pref][col]
    }
    
    Vnext <- Vt
    print(t)
    
    # Record the policy for the first year
    if (t <= 12) {
      D[, , , m] <- array(best_action,
                          dim = c(Nstates_nWild, Nstates_nCaptive, Nstates_exposure))
    }
  }
 
  sdp_Output <- list(D = D,
                     transitions_pop = transitions_pop,
                     arguments = list(fecundity_w = fecundity_w, 
                                      fecundity_c = fecundity_c, 
                                      breeding_month = breeding_month,
                                      mortality_w_base = mortality_w_base, 
                                      mortality_w_flu = mortality_w_flu, 
                                      mortality_c = mortality_c,
                                      phi_wc = phi_wc, 
                                      phi_cw = phi_cw,
                                      disease_transitions = disease_transitions,
                                      discount = discount,
                                      Tmax = Tmax,
                                      captive_cost = captive_cost,
                                      establish_cost = establish_cost,
                                      cw_translocate_cost = cw_translocate_cost,
                                      wc_translocate_cost = wc_translocate_cost,
                                      cost_weight = cost_weight))
  
  return(sdp_Output)
}