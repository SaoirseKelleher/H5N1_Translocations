# Runs an SDP optimisation for the bushfire example

#' @param mortality_fire Daily mortality in the wild if a fire occurs
#' @param mortality_fire_sd Variation in fire mortality
#' @param phi_wc Wild to captive translocation success rate (surviving transition + first month in captivity)
#' @param phi_wc_sd Variation in wild to captive translocation success rate
#' @param fire_transitions (4*4 matrix of fire transitions)
#' @param max_action Maximum number of individuals to translocate in a day
#' @param Tmax Time horizon in days (defaults to 60 days)
#' @param discount Discount rate (defaults to 0.99)
#' @param nsamples Number of draws from joint probability distribution
#' @param seed Seed for random draws of joint probability distribution
#' @returns Array of optimal decisions given current state.


run_bushfire_sdp <- function(
  mortality_fire, mortality_fire_sd, 
  phi_wc, phi_wc_sd, 
  fire_transitions,
  max_action,
  discount = 0.999,
  Tmax = 360,
  nsamples = 10000,
  seed = 1512
){


  # Define states ---------------------------------------------------
  # Define the states in the model. There are three dimensions: the wild
  # population, the captive population, and fire exposure.

  # Wild pop: 0-2000, increments of 50
  Nstates_nWild <- 41
  states_nWild <- 1:41
  
  # Captive pop: 0-1000, increments of 50
  Nstates_nCaptive <- 21
  states_nCaptive <- 1:21

  # Exposure level to fire: 4 is 'burning habitat'
  Nstates_exposure <- 4
  states_exposure <- 1:4

  # Binary 'burning' states
  Nstates_burning <- 2
  states_burning <- c(FALSE, TRUE)

  # Total number of states:
  # Population only:
  nStates_pop <- Nstates_nWild*Nstates_nCaptive
  # Total:
  nStates_full <- Nstates_nWild*Nstates_nCaptive*Nstates_exposure


  # Define possible actions -------------------------------------------
  # Positive actions are translocation captive to wild, negative actions are 
  # translocation wild to captive
  actions <- seq(-100, 0, by = 5)
  
  # Total number of actions
  nActions <- length(actions)
  
  # Define population model --------------------------------------------
  # Draw random parameter values

  withr::with_seed(seed, {
    mortality_fire_sample <- truncnorm::rtruncnorm(nsamples, 0, 1, 
      mean = mortality_fire,
      sd = mortality_fire_sd)
    
    phi_wc_sample <- truncnorm::rtruncnorm(nsamples, 0, 1, 
      mean = phi_wc,
      sd = phi_wc_sd)
  })

  # This function determines the next population of the wild and captive 
  # population given an action.
  dynamic <- function(actualWild, actualCaptive, action, fire){
    
    nextWild <- rep(NA, nsamples)
    nextCaptive <- rep(NA, nsamples)
    
    tau_wc <- ifelse(action < 0, abs(action), 0)/100

    mortality_w_sample <- ifelse(fire == 1, mortality_fire_sample, 0)
    
    nextWild <- 
      ((actualWild-(actualWild*tau_wc))*(1-mortality_w_sample))
    
    nextCaptive <- 
      actualCaptive + 
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
                          fire,
                          action){
    

    maxPopUtility <- max(state_util)
    
    
    popUtility <- state_util*transitions_pop[lastWildPop,
                                             lastCaptivePop,
                                             ,
                                             ,
                                             fire+1,
                                             which(actions==action)]
      
    totalUtility <- (sum(popUtility)/maxPopUtility)
    
    return(totalUtility)
  }


  # Fill population transitions array --------------------------------------
  # Population transitions: Prob of going from one state to another. Determines
  # only future population, not exposure states.
  transitions_pop <- array(0, 
                           dim = c(Nstates_nWild, # W pop at t
                                   Nstates_nCaptive, # C pop at t
                                   Nstates_nWild, # W pop at t+1
                                   Nstates_nCaptive, # C pop at t+1
                                   2, # burning/not-burning at t
                                   nActions)) # action
  
  # Utility array: Utility of an action in a given state
  utility <- array(0, 
                   dim = c(Nstates_nWild,
                           Nstates_nCaptive,
                           Nstates_exposure,
                           nActions))
  
  # Iterate over all possible population states
  for (w in 1:Nstates_nWild){
    for (x in 1:Nstates_nCaptive){
      for (y in 1:Nstates_burning){
        # And possible actions
        for (i in 1:nActions){
          transitions_pop[w, x, , , y, i] <- dynamic(
            actualWild = (states_nWild[w]-1)*50, 
            actualCaptive = (states_nCaptive[x]-1)*50,
            fire = states_burning[y],
            action = actions[i]
          )
        }
      }
    }
  }
  
  # Fill utility matrix
  for (w in 1:Nstates_nWild){
    for (x in 1:Nstates_nCaptive){
      for (y in 1:Nstates_burning){
        # And possible actions
        for (i in 1:nActions){
          # Compute utility of action - utility uses exposure states rather
          # burning states, hence the below logic to fill it fully.
          if (y == 1){
            utility[w, x, 1:3, i] <- get_utility(state_util,
              w,
              x,
              transitions_pop,
              FALSE,
              actions[i])
          }
          if (y == 2){
            utility[w, x, 4, i] <- get_utility(state_util,
            w,
            x,
            transitions_pop,
            TRUE,
            actions[i])
          }
        }
      }
    }
  }
  
  # Reshape arrays for matrix algebra ---------------------------------------
  # Number of population states
  nPop <- Nstates_nWild*Nstates_nCaptive 
  # Number of state/action combinations
  nRow <- nPop*nActions

  # Transitions: one (state x action) by next-state matrix per
  # breeding state (1/2) 
  PT <- vector("list", 2)
  dim(PT) <- c(2, 1)
  for (burn in 1:2) {
    slice <- transitions_pop[, , , , burn, ]   # 41 x 21 x 41 x 21 x nActions
    dim(slice) <- c(nPop, nPop, nActions)        # s x s' x i
    PT[[burn, 1]] <- matrix(aperm(slice, c(1, 3, 2)), nrow = nRow)  # (s,i) x s'
  }
  rm(slice)

  # Utility: one (state x action) by exposure matrix per month
  u <- utility[, , , ]                        # 41 x 21 x 6 x nActions
  dim(u) <- c(nPop, Nstates_exposure, nActions)
  u <- matrix(aperm(u, c(1, 3, 2)), nrow = nRow)
  
  # Tie-break order: smallest |action| first, then the positive one.
  # max.col(ties.method = "first") then reproduces the original rule.
  pref <- order(abs(actions), -actions)
  

  # Solve Bellman via backwards iteration -----------------------------------
  # Terminal value: (w-1)*2 + (x-1), identical across exposure states
  pop_value <- as.vector(outer((states_nWild - 1) * 2, states_nCaptive - 1, "+"))/100
  Vnext <- matrix(pop_value, nrow = nPop, ncol = Nstates_exposure)
  
  D <- array(0, 
    dim = c(Nstates_nWild,
    Nstates_nCaptive,
    Nstates_exposure))
  Vt          <- matrix(0, nPop, Nstates_exposure)
  best_action <- matrix(0, nPop, Nstates_exposure)
  
  for (t in (Tmax - 1):1) {
    
    # Disease step: expected next value for each current exposure y
    # M[s', y] = sum_z disease_transitions[y, z] * Vnext[s', z]
    M <- Vnext %*% t(fire_transitions)          # nPop x 6
    
    # Population step: all states and actions in two matrix products
    Q <- u
    Q[, 1:3] <- Q[, 1:3] + discount * (PT[[1, 1]] %*% M[, 1:3])  # burning
    Q[, 4]   <- Q[, 4]   + discount * (PT[[2, 1]] %*% M[, 4])    # not burning
    
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
    
    
    D[, , ] <- array(best_action,
      dim = c(Nstates_nWild, Nstates_nCaptive, Nstates_exposure))
    }
  }
 