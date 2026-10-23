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
                      disease, mortality_w_base, 
                      mortality_w_flu, mortality_c, 
                      breeding, fecundity_w, fecundity_c, 
                      phi_wc, phi_cw){
    
    tau_cw <- ifelse(action > 0, action, 0)/100
    tau_wc <- ifelse(action < 0, abs(action), 0)/100
    mortality_w <- ifelse(disease == 1, mortality_w_flu, mortality_w_base)
    
    nextWild <- 
      ((actualWild-(actualWild*tau_wc))*(1-mortality_w)) + 
      ((actualWild-(actualWild*tau_wc))*breeding*fecundity_w) +
      (actualCaptive*tau_cw*phi_cw*(1-mortality_c))
    
    nextCaptive <- 
      ((actualCaptive-(actualCaptive*tau_cw))*(1-mortality_c)) + 
      ((actualCaptive-(actualCaptive*tau_cw))*breeding*fecundity_c) +
      (actualWild*tau_wc*phi_wc*(1-mortality_w))
    
    nextpop <- list(Wild = nextWild,
                    Captive = nextCaptive)
    
    return(nextpop)
  }
  

  # Define utility function -------------------------------------------------
  # Define the value of an action given the next wild/captive populations. Here,
  # wild individuals are worth 2x captive individuals.
  get_utility <- function(lastWildPop, 
                          lastCaptivePop, 
                          transitions_pop,
                          infected,
                          breeding,
                          action) {
    
    
    state_util <- matrix(0,
                         nrow = 41,
                         ncol = 21)
    for(w in 1:41){
      for(x in 1:21){
        state_util[w,x] <- ((w-1)*2)+((x-1))
      }
    }
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
      (((lastCaptivePop-1)*tau_cw)*cw_translocate_cost)+
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
            
            # Calculate next population given the state & action
            nextpop <- dynamic(actualWild = (states_nWild[w]-1)*50, 
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
            
            # Round population projections into states. Add 1 so zero population 
            # goes into state 1.
            nextWildState <- round(nextpop$Wild/50)+1
            nextWildState <- ifelse(nextWildState>41, 41, nextWildState)
            nextCaptiveState <- round(nextpop$Captive/50)+1
            nextCaptiveState <- ifelse(nextCaptiveState>21, 21, nextCaptiveState)
            
            # Calculate probabilities of next states
            probDensity <- function(x, n){
              sum(dpois((((x*50)+1):((x*50)+50)), n))
            }
            
            wildMatrix <- matrix(c(dpois(0, nextpop$Wild), 
                                   sapply(0:38, probDensity, nextpop$Wild),
                                   sum(dpois(3950:10000, nextpop$Wild))), 
                                 nrow = 41, 
                                 ncol = 21)
            
            captiveMatrix <- matrix(c(dpois(0, nextpop$Captive), 
                                      sapply(0:18, probDensity, nextpop$Captive),
                                      sum(dpois(1950:10000, nextpop$Captive))), 
                                    nrow = 41, 
                                    ncol = 21,
                                    byrow = TRUE)
            
            nextMatrix <- wildMatrix*captiveMatrix

            transitions_pop[w, x, , , y, z, i] <- nextMatrix
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
  
  # Solve Bellman via backwards iteration -----------------------------------
  # Calculate value of being in each state at TMax
  Vtmax <- array(0,
                 dim = c(Nstates_nWild, 
                         Nstates_nCaptive, 
                         Nstates_exposure, 
                         Nstates_month))
  
  # Value in the final state does not consider cost.
  for (i in 1:Nstates_nWild){
    for (j in 1:Nstates_nCaptive){
      Vtmax[i,j,,] <- ((states_nWild[i]-1)*2) + (states_nCaptive[j]-1)
    }
  }

  # Create empty arrays for action values at t and t+1
  Vt <-  array(0,
               dim = c(Nstates_nWild,
                       Nstates_nCaptive, 
                       Nstates_exposure, 
                       Nstates_month))
  Vtplus <-  Vtmax 
  
  # Create empty optimal policy array for each state
  D <- array(0,
             dim = c(Nstates_nWild, 
                     Nstates_nCaptive, 
                     Nstates_exposure, 
                     Nstates_month))
  
  # Starting at the maximum timestep and iterate backwards until the first time-step.
  for (t in (Tmax-1):1){
    
    # Define array Q, stores the updated action values for all states and actions
    Q <- array(0, 
               dim = c(Nstates_nWild,
                       Nstates_nCaptive,
                       Nstates_exposure,
                       Nstates_month,
                       nActions))
    # Define array W, the probability of each future state given conditions at t
    W <-  array(0, 
                dim = c(Nstates_nWild,
                        Nstates_nCaptive,
                        Nstates_exposure))
    
    # For each action, loop through all possible states (excluding month, which
    # is deterministic and unrelated to actions/populations)
    for (i in 1:length(actions)) { # Action taken at t
      for (w in 1:Nstates_nWild){ # W pop at t
        for (x in 1:Nstates_nCaptive){ # C pop at t
          for (y in 1:Nstates_exposure){ # Exposure at t
            # Get probability of each state at t+1 given state at t
            for (z in 1:Nstates_exposure){
              # Population transitions & disease transition are independent, can just multiply here
              W[,,z]  <- transitions_pop[w,x,,,(y==6)+1,
                                         (z==breeding_month)+1,i]*disease_transitions[y,z]
            }
            
            # Fill Q with utility - first fill each month with just the utility... 
            for (t2 in 1:12){
              Q[w,x,y,t2,i] <- utility[w,x,y,Tmonths[t],i]
            } 
            # ... but then replace the actual next month with the utility plus the future state potential.
            Q[w,x,y,Tmonths[t],i] <- utility[w,x,y,Tmonths[t],i] + sum(discount*W[,,]*Vtplus[,,,Tmonths[t+1]])
          }
        }
      } 
    }
    
    # Find the value given the optimal action for each state at this timestep
    Vt <- apply(Q, c(1,2,3,4), max)
    # Set the next timestep to be the value of the current timestep
    Vtplus <- Vt
    # Print timestep to track progress
    print(t)
    # For the last 12 timesteps (i.e., the first year from the present),
    # indicate the optimal action for each state
    if (t <= 12){
      for (w in 1:Nstates_nWild){
        for (x in 1:Nstates_nCaptive){
          for (y in 1:Nstates_exposure){
            # Get the action that provides the best future value. 
            # If tied, pick the one that involves the least translocation
            D[w,x,y,Tmonths[t]] <- max(actions[(which(Q[w,x,y,Tmonths[t],] == Vt[w,x,y,Tmonths[t]]))][which(abs(actions[(which(Q[w,x,y,Tmonths[t],] == Vt[w,x,y,Tmonths[t]]))])==min(abs(actions[(which(Q[w,x,y,Tmonths[t],] == Vt[w,x,y,Tmonths[t]]))])))])
          }
        }
      }
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