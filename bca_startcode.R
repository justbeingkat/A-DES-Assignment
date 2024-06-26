## Advanced Discrete Event Simulation
## Assignment Part I: Health Economic Model
##
## Good luck!

## Section 1: Initialization ----

# Clear the workspace
rm(list=ls()); gc();

# Load the required packages, make sure the required packages are installed. For the installation of the packages use line 13 (uncomment shift+Ctrl+C, or remove #)

# install.packages(c("simmer", "simmer.plot", "fitdistrplus")) # Install packages

library(simmer);
library(simmer.plot);
library(fitdistrplus);
library(dplyr); 
library(ggplot2); 

# Set the working directory
#setwd("~/A-DES-Assignment");

# Load funtions for extracting monitored attributes
source("getSingleAttribute.R", echo=T);
source("getMultipleAttributes.R", echo=T);


## Section 2: Data analysis ----

# Load the dataset
load("trial_dataset.RData");

# Define parameters
# Patient characteristics 
p.male <- 0.34;                               # probability of being male
p.women <- 1 - 0.34; 
p.poor <- 0.1986242;                               # probability of a poor condition
m.age.men <- 57.25063;                        # mean age for men
sd.age.men <- 8.05432;                        # sd age for men
m.age.women <- 65.63482;                      # mean age for women
sd.age.women <- 13.27633;                     # sd age for women

x.w <-seq(20,100, by= 1)
x.m <-seq(20,100, by= 1)

age.woman.d <- dnorm(x.w, m.age.women, sd.age.women); 
age.men.d <- dnorm(x.m, m.age.men, sd.age.men);

# Tx general
#t.cycle <- 30;    # average time of a normal treatment cycle
#t.major <- 6;     # average time of a cycle in which major complications occur
#t.death <- 15;    # average time of a cycle in which the patient dies

t.normal <- function(){return(rweibull(n = 1, shape = 16, scale = 31))} #Time to event step 1
t.minor <- function(){return(rweibull(n = 1, shape = 13, scale = 31))} #Time to event step 2
t.major <- function(){return(rweibull(n = 1, shape = 5, scale = 6))} #Time to event step 1
t.death <- function(){return(rweibull(n = 1, shape = 2, scale = 18))} #Time to event step 2

#p.minor <- 0.10;  # probability of minor complications in a cycle
#p.major <- 0.04;  # probability of major complications in a cycle
#p.death <- 0.03;  # probability of death in a cycle

#responsive probabillities 
p.minor.responsive.poor.Tx1 <-0.066666667
p.minor.responsive.good.Tx1 <- 0.038934426

p.major.responsive.poor.Tx1 <- 0.077777778
p.major.responsive.good.Tx1 <-0.024590164

p.death.responsive.poor.Tx1 <- 0.033333333
p.death.responsive.good.Tx1 <- 0.028688525

#not responsive probabillities
p.minor.not.responsive.poor.Tx1 <-0.226950355
p.minor.not.responsive.good.Tx1 <- 0.162162162

p.major.not.responsive.poor.Tx1 <- 0.141843972
p.major.not.responsive.good.Tx1 <-0.11261261

p.death.not.responsive.poor.Tx1 <- 0.021276596
p.death.not.responsive.good.Tx1 <- 0.031531532

c.minor <- 487; 
c.major <- 14897;
u.minor <- 0.051;
u.major <- 0.089

# Tx1 specific
c.Tx1.cycle <- 351;   # costs of a cycle of Tx1
c.Tx1.day <- 8.50;       # additional daily costs when on treatment Tx1
u.Tx1 <- 0.55/36;    # utility per day when on treatment Tx1
u.Tx1.Response <- function(){return(rnorm(n = 1, mean = 0.4695, sd = 0.1214))}; 
u.Tx1.NoResponse <- function(){return(rnorm(n = 1, mean = 0.46423, sd = 0.14172))}; 

p.Tx1.poor <- mean(data$Tx1.C1.Dx.Pet[data$Poor==1]==1, na.rm=T);                               # probability of effective Tx1 treatment when in poor condition
p.Tx1.good <- mean(data$Tx1.C1.Dx.Pet[data$Poor==0]==1, na.rm=T);                               # probability of effective Tx1 treatment when in good condition

# Tx2 specific
c.Tx2.cycle <- 4141;  # costs of a cycle of Tx2
c.Tx2.day <- 19;      # additional dayly costs when on treatment Tx2
u.Tx2 <- 0.5/365;     # utility per day when on treatment Tx2

u.Tx2.Response <- function(){return(rnorm(n = 1, mean = 0.53106, sd = 0.1316209))};
u.Tx2.NoResponse <- function(){return(rnorm(n = 1, mean = 0.54499, sd = 0.113477))};

p.Tx2.respondedTx1 <- 0.3493
p.Tx2.notrespondedTx1 <- 0.6771

# p.Tx2.yes.exp and p.Tx2.no.exp are provided for step 2.4
p.Tx2.yes.exp <- seq(from=0.39, to=0.47, by=0.02);                                               # probability of effective Tx2 treatment when responded to Tx1, dependent on cycles Tx1 
p.Tx2.no.exp <- seq(from=0.87, to=0.31, by=-0.14);                                               # probability of effective Tx2 treatment when not responded to Tx1, dependent on cycles Tx1 

# FU1 en FU2 specific
#t.fu1.full <- 63          # average time spent in the first follow up if the patient survives during follow up
#t.fu1.death <- 42         # average time spent in the first follow up if the patient dies during follow up
#t.fu2 <- 100              # average time spent in the second follow up after Tx2
p.death.followup <- 0.05; # probability of dying during first follow up

t.fu1.normal <- function(){rexp(1, rate = 0.0167)}
t.fu1.dead <- function(){rexp(1, rate = 0.0273)}
t.fu2 <- function(){rexp(1, rate=0.0072)}

## Section 3: Supportive functions ----

# Functions for determining the event to happen

# Function for defining the event during a cycle of treatment
Tx.event.cycle <- function(response, poor) {
  poor <- ifelse(runif(1) < p.poor, 1, 0)
  if (response == 1){
    if (poor == 1){
      p.complete <- 1 - p.minor.responsive.poor.Tx1 - p.major.responsive.poor.Tx1 - p.death.responsive.poor.Tx1  
      event <- sample(1:4, 1, prob = c(p.complete, p.minor.responsive.poor.Tx1, p.major.responsive.poor.Tx1, p.death.responsive.poor.Tx1))
      return(event)
    } else if (poor == 0) {
      p.complete <- 1 - p.minor.responsive.good.Tx1 - p.major.responsive.good.Tx1 - p.death.responsive.good.Tx1  
      event <- sample(1:4, 1, prob = c(p.complete, p.minor.responsive.good.Tx1, p.major.responsive.good.Tx1, p.death.responsive.good.Tx1))
      return(event)
    }
  } else if (response == 0){
    if (poor == 1){
      p.complete <- 1 - p.minor.not.responsive.poor.Tx1 - p.major.not.responsive.poor.Tx1 - p.death.not.responsive.poor.Tx1  
      event <- sample(1:4, 1, prob = c(p.complete, p.minor.not.responsive.poor.Tx1, p.major.not.responsive.poor.Tx1, p.death.not.responsive.poor.Tx1))
      return(event)
    } else if (poor == 0) {
      p.complete <- 1 - p.minor.not.responsive.good.Tx1 - p.major.not.responsive.good.Tx1 - p.death.not.responsive.good.Tx1  
      event <- sample(1:4, 1, prob = c(p.complete, p.minor.not.responsive.good.Tx1, p.major.not.responsive.good.Tx1, p.death.not.responsive.good.Tx1))
      return(event)
    }
  }
}
#Function for defining the event during Follow up 1
Tx1.event.fu1 <- function() {
  
  #select event during follow up
  event <- ifelse(runif(1) < p.death.followup, 2, 1)
  
  return(event)
}

Tx1.utility <- function(response) {
  if (response == 1){
    return (u.Tx1.Response()/365)
  } else if (response == 0){
    return(u.Tx1.NoResponse()/365)
  }
}

Tx2.utility <- function(response) {
  if (response == 1){
    return (u.Tx2.Response()/365)
  } else if (response == 0){
    return(u.Tx2.NoResponse()/365)
  }
}
# Functions for determining the time-to-events

# Function for defining the time spent on a cycle
Tx.time.cycle <- function(event) {
  if (event == 1) {
    
    return(t.normal)  # Duration for completing the cycle without any complication
  } else if (event == 2) {
    return(t.)  # Duration for minor complications
  } else if (event == 3) {
    return(t.major)  # Duration for major complications during treatment
  } else if (event == 4) {
    return(t.death)  # Duration for death
  }
}

Tx1.time.fu1 <- function(event) {
  if (event == 1) {
    return(t.fu1.full)  # Duration for completing the followup
  } else if (event == 2) {
    return(t.fu1.death)  # Duration for death during followup
  }
}
Tx1.Response <- function(poor) {
  poor <- ifelse(runif(1) < p.poor, 1, 0)
  if (poor == 0){
    return(ifelse(runif(1) < p.Tx1.good, 1, 0))
  } else if (poor == 1){
    return(ifelse(runif(1) < p.Tx1.poor, 1, 0))
  }
}

Tx2.Response <- function(response) {
  if (response == 1){
    return(ifelse(runif(1) < p.Tx2.respondedTx1, 1, 0))
  } else if (response == 0){
    return(ifelse(runif(1) < p.Tx2.notrespondedTx1, 1, 0))
  }
}
## Section 4: Discrete event simulation model ----
# Define the model structure for the current practice, i.e., best standard care (BSC)
bsc.model <- trajectory() %>%
  
  # Initialization
  set_attribute(key = "Alive", value = 1) %>%                                                                           # define an attribute to check whether the patient is alive
  set_attribute(key = "Tx1.Response", value = function() Tx1.Response()) %>%                                                          # check whether the patient responded to Tx1
  set_attribute(key = "Tx2.Response", value = function() Tx2.Response(Tx1.Response())) %>%       # check whether the patient responded to Tx2
  set_attribute(key = "Total.Costs", value = 0) %>%
  set_attribute(key = "Total.Utility", value = 0) %>%
  set_attribute(key = "Cycle Count Tx1", value = 0) %>%
  set_attribute(key = "Cycle Count Tx2", value = 0) %>%
  # First-line treatment
  set_attribute(key = "Tx.event.cycle", value = function() Tx.event.cycle(get_attribute(bsc.sim,"Tx1.Response"))) %>%                                         # select the event to happen in this treatment cycle
  branch(option = function() get_attribute(bsc.sim, "Tx.event.cycle"), continue = c(T, T, T, F),
         
         # Event 1: Full cycle
         trajectory() %>%
           set_attribute(keys = "Cycle Count Tx1", mod = "+", values = 1) %>%
           set_attribute(key = "Tx.time.cycle", value = function() t.normal()) %>%  # determine how long the cycle will last
           seize(resource = "Tx1", amount = 1) %>%         
           timeout_from_attribute(key = "Tx.time.cycle") %>%                                                          # stay in first-line treatment for the determined time
           release(resource = "Tx1", amount = 1) %>%                                                                  # leave first-line treatment
           set_attribute(keys = "cycle.costs", value = c.Tx1.day) %>%
           set_attribute(keys = "cycle.costs", mod = "*", values = function() get_attribute(bsc.sim, "Tx.time.cycle")) %>% #determine the costs
           set_attribute(keys = "cycle.costs", mod = "+", values = c.Tx1.cycle) %>% 
           set_attribute(keys = "Total.Costs", mod = "+", values = function() get_attribute(bsc.sim, "cycle.costs")) %>%
           set_attribute(key = "cycle.utility", value = function() Tx1.utility(get_attribute(bsc.sim, "Tx1.Response"))) %>%
           set_attribute(keys = "cycle.utility", mod = "*", values = function() get_attribute(bsc.sim, "Tx.time.cycle")) %>%
           set_attribute(keys = "Total.Utility", mod = "+", values = function() get_attribute(bsc.sim, "cycle.utility")) %>%
           rollback(target = 14, times=5, check = function() get_attribute(bsc.sim, "Cycle Count Tx1") < 5),                                                                          # go back for another cycle
         
         # Event 2: Minor complication
         trajectory() %>%
           set_attribute(keys = "Cycle Count Tx1", mod = "+", values = 1) %>%
           set_attribute(key = "Tx.time.cycle", value = function() t.minor()) %>%  # determine how long the cycle will last
           seize(resource = "Tx1", amount = 1) %>%                                                                    # occupy a place in first-line treatment
           timeout_from_attribute(key = "Tx.time.cycle") %>%                                                          # stay in first-line treatment for the determined time
           release(resource = "Tx1", amount = 1) %>%                                                                  # leave first-line treatment
           set_attribute(keys = "cycle.costs", value = c.Tx1.day) %>%
           set_attribute(keys = "cycle.costs", mod = "*", values = function() get_attribute(bsc.sim, "Tx.time.cycle")) %>% #determine the costs
           set_attribute(keys = c("cycle.costs", "cycle.costs"), mod = "+", values = c(c.Tx1.cycle, c.minor)) %>% 
           set_attribute(keys = "Total.Costs", mod = "+", values = function() get_attribute(bsc.sim, "cycle.costs")) %>%
           set_attribute(key = "cycle.utility", value = function() Tx1.utility(get_attribute(bsc.sim, "Tx1.Response"))) %>%           
           set_attribute(keys = "cycle.utility", mod = "*", values = function() get_attribute(bsc.sim, "Tx.time.cycle")) %>%
           set_attribute(keys = "cycle.utility", mod = "+", values = -u.minor) %>%
           set_attribute(keys = "Total.Utility", mod = "+", values = function() get_attribute(bsc.sim, "cycle.utility")) %>%
           rollback(target = 15, times=5, check = function() get_attribute(bsc.sim, "Cycle Count Tx1") < 5),    
         
         # Event 3: Major complication
         trajectory() %>%
           set_attribute(keys = "Cycle Count Tx1", mod = "+", values = 1) %>%
           set_attribute(key = "Tx.time.cycle", value = function() t.major()) %>%  # determine how long the cycle will last
           seize(resource = "Tx1", amount = 1) %>%                                                                    # occupy a place in first-line treatment
           timeout_from_attribute(key = "Tx.time.cycle") %>%                                                          # stay in first-line treatment for the determined time
           release(resource = "Tx1", amount = 1) %>%                                                                  # leave first-line treatment
           set_attribute(keys = "cycle.costs", value = c.Tx1.day) %>%
           set_attribute(keys = "cycle.costs", mod = "*", values = function() get_attribute(bsc.sim, "Tx.time.cycle")) %>% #determine the costs
           set_attribute(keys = c("cycle.costs", "cycle.costs"), mod = "+", values = c(c.Tx1.cycle, c.major)) %>% 
           set_attribute(keys = "Total.Costs", mod = "+", values = function() get_attribute(bsc.sim, "cycle.costs")) %>%
           set_attribute(key = "cycle.utility", value = function() Tx1.utility(get_attribute(bsc.sim, "Tx1.Response"))) %>%         
           set_attribute(keys = "cycle.utility", mod = "*", values = function() get_attribute(bsc.sim, "Tx.time.cycle")) %>%
           set_attribute(keys = "cycle.utility", mod = "+", values = -u.major) %>%
           set_attribute(keys = "Total.Utility", mod = "+", values = function() get_attribute(bsc.sim, "cycle.utility")),
         
         # Event 4: Death
         trajectory() %>%
           set_attribute(keys = "Cycle Count Tx1", mod = "+", values = 1) %>%
           set_attribute(key = "Tx.time.cycle", value = function() t.death()) %>%  # determine how long the cycle will last
           seize(resource = "Tx1", amount = 1) %>%                                                                    # occupy a place in first-line treatment
           timeout_from_attribute(key = "Tx.time.cycle") %>%                                                          # stay in first-line treatment for the determined time
           release(resource = "Tx1", amount = 1) %>%                                                                  # leave first-line treatment
           set_attribute(keys = "cycle.costs", value = c.Tx1.day) %>%
           set_attribute(keys = "cycle.costs", mod = "*", values = function() get_attribute(bsc.sim, "Tx.time.cycle")) %>% #determine the costs
           set_attribute(keys = c("cycle.costs"), mod = "+", values = c(c.Tx1.cycle)) %>% 
           set_attribute(keys = "Total.Costs", mod = "+", values = function() get_attribute(bsc.sim, "cycle.costs")) %>%
           set_attribute(key = "cycle.utility", value = function() Tx1.utility(get_attribute(bsc.sim, "Tx1.Response"))) %>%           
           set_attribute(keys = "cycle.utility", mod = "*", values = function() get_attribute(bsc.sim, "Tx.time.cycle")) %>%
           set_attribute(keys = "Total.Utility", mod = "+", values = function() get_attribute(bsc.sim, "cycle.utility")) %>%
           set_attribute(key = "Alive", value = 0)                                                                     # update that the patient has died
  ) %>%
  
  # Follow up 1
  set_attribute(key = "Tx1.event.fu1", value = function() Tx1.event.fu1()) %>%                                            # select the event to happen in this treatment cycle
  branch(option = function() get_attribute(bsc.sim, "Tx1.event.fu1"), continue = c(T, F),
         
         # Event 1: Survives follow-up
         trajectory() %>%
           set_attribute(key = "Tx.time.fu1", value = function() t.fu1.normal()) %>%  # determine how long the follow-up will last
           timeout_from_attribute(key = "Tx.time.fu1") %>%
           set_attribute(keys = "Total.Utility", mod = "*", values = 1.1),                                                                # stay in follow-up treatment for the determined time
         
         # Event 2: Dies in follow-up
         trajectory() %>%
           set_attribute(key = "Tx.time.fu1", value = function() t.fu1.dead()) %>%  # determine how long the follow-up will last
           timeout_from_attribute(key = "Tx.time.fu1") %>%                                                             # stay in follow-up treatment for the determined time
           set_attribute(keys = "Total.Utility", mod = "*", values = 1.1) %>%                                                              # stay in follow-up treatment for the determined time
           set_attribute(key = "Alive", value = 0)
  ) %>%
  
  #Second-line treatment
  set_attribute(key = "Tx.event.cycle", value = function() Tx.event.cycle(get_attribute(bsc.sim,"Tx2.Response"))) %>%                                         # select the event to happen in this treatment cycle
  branch(option = function() get_attribute(bsc.sim, "Tx.event.cycle"), continue = c(T, T, T, F),
         
         # Event 1: Full cycle
         trajectory() %>%
           set_attribute(key = "Cycle Count Tx2", value = 0) %>%
           set_attribute(key = "Tx.time.cycle", value = function() t.normal()) %>%  # determine how long the cycle will last
           seize(resource = "Tx2", amount = 1) %>%                                                                    # occupy a place in second-line treatment
           timeout_from_attribute(key = "Tx.time.cycle") %>%                                                          # stay in second-line treatment for the determined time
           release(resource = "Tx2", amount = 1) %>%                                                                  # leave second-line treatment
           set_attribute(keys = "cycle.costs", value = c.Tx2.day) %>%
           set_attribute(keys = "cycle.costs", mod = "*", values = function() get_attribute(bsc.sim, "Tx.time.cycle")) %>% #determine the costs
           set_attribute(keys = c("cycle.costs"), mod = "+", values = c(c.Tx2.cycle)) %>% 
           set_attribute(keys = "Total.Costs", mod = "+", values = function() get_attribute(bsc.sim, "cycle.costs")) %>%
           set_attribute(key = "cycle.utility", value = function() Tx2.utility(get_attribute(bsc.sim, "Tx2.Response"))) %>%           
           set_attribute(keys = "cycle.utility", mod = "*", values = function() get_attribute(bsc.sim, "Tx.time.cycle")) %>%
           set_attribute(keys = "Total.Utility", mod = "+", values = function() get_attribute(bsc.sim, "cycle.utility")) %>%
           rollback(target = 14, times=5, check = function() get_attribute(bsc.sim, "Cycle Count Tx2") < 5),                                                                          # go back for another cycle
         
         # Event 2: Minor complication
         trajectory() %>%
           set_attribute(key = "Cycle Count Tx2", value = 0) %>%
           set_attribute(key = "Tx.time.cycle", value = function() t.minor()) %>%  # determine how long the cycle will last
           seize(resource = "Tx2", amount = 1) %>%                                                                    # occupy a place in second-line treatment
           timeout_from_attribute(key = "Tx.time.cycle") %>%                                                          # stay in second-line treatment for the determined time
           release(resource = "Tx2", amount = 1) %>%                                                                  # leave second-line treatment
           set_attribute(keys = "cycle.costs", value = c.Tx2.day) %>%
           set_attribute(keys = "cycle.costs", mod = "*", values = function() get_attribute(bsc.sim, "Tx.time.cycle")) %>% #determine the costs
           set_attribute(keys = c("cycle.costs", "cycle.costs"), mod = "+", values = c(c.Tx2.cycle, c.minor)) %>% 
           set_attribute(keys = "Total.Costs", mod = "+", values = function() get_attribute(bsc.sim, "cycle.costs")) %>%
           set_attribute(key = "cycle.utility", value = function() Tx2.utility(get_attribute(bsc.sim, "Tx2.Response"))) %>%           
           set_attribute(keys = "cycle.utility", mod = "*", values = function() get_attribute(bsc.sim, "Tx.time.cycle")) %>%
           set_attribute(keys = "cycle.utility", mod = "+", values = -u.minor) %>%
           set_attribute(keys = "Total.Utility", mod = "+", values = function() get_attribute(bsc.sim, "cycle.utility")) %>%
           rollback(target = 15, times=5, check = function() get_attribute(bsc.sim, "Cycle Count Tx2") < 5),   
         
         # Event 3: Major complication
         trajectory() %>%
           set_attribute(key = "Cycle Count Tx2", value = 0) %>%
           set_attribute(key = "Tx.time.cycle", value = function() t.major()) %>%  # determine how long the cycle will last
           seize(resource = "Tx2", amount = 1) %>%                                                                    # occupy a place in second-line treatment
           timeout_from_attribute(key = "Tx.time.cycle") %>%                                                          # stay in second-line treatment for the determined time
           release(resource = "Tx2", amount = 1) %>%                                                                  # leave second-line treatment
           set_attribute(keys = "cycle.costs", value = c.Tx2.day) %>%
           set_attribute(keys = "cycle.costs", mod = "*", values = function() get_attribute(bsc.sim, "Tx.time.cycle")) %>% #determine the costs
           set_attribute(keys = c("cycle.costs", "cycle.costs"), mod = "+", values = c(c.Tx2.cycle, c.major)) %>% 
           set_attribute(keys = "Total.Costs", mod = "+", values = function() get_attribute(bsc.sim, "cycle.costs")) %>%
           set_attribute(key = "cycle.utility", value = function() Tx2.utility(get_attribute(bsc.sim, "Tx2.Response"))) %>%
           set_attribute(keys = "cycle.utility", mod = "*", values = function() get_attribute(bsc.sim, "Tx.time.cycle")) %>%
           set_attribute(keys = "cycle.utility", mod = "+", values = -u.major) %>%
           set_attribute(keys = "Total.Utility", mod = "+", values = function() get_attribute(bsc.sim, "cycle.utility")),
         
         # Event 4: Death
         trajectory() %>%
           set_attribute(key = "Cycle Count Tx2", value = 0) %>%
           set_attribute(key = "Tx.time.cycle", value = function() t.death()) %>%  # determine how long the cycle will last
           seize(resource = "Tx2", amount = 1) %>%                                                                    # occupy a place in second-line treatment
           timeout_from_attribute(key = "Tx.time.cycle") %>%                                                          # stay in second-line treatment for the determined time
           release(resource = "Tx2", amount = 1) %>%                                                                  # leave second-line treatment
           set_attribute(keys = "cycle.costs", value = c.Tx2.day) %>%
           set_attribute(keys = "cycle.costs", mod = "*", values = function() get_attribute(bsc.sim, "Tx.time.cycle")) %>% #determine the costs
           set_attribute(keys = c("cycle.costs"), mod = "+", values = c(c.Tx2.cycle)) %>% 
           set_attribute(keys = "Total.Costs", mod = "+", values = function() get_attribute(bsc.sim, "cycle.costs")) %>%
           set_attribute(key = "cycle.utility", value = function() Tx2.utility(get_attribute(bsc.sim, "Tx2.Response"))) %>% 
           set_attribute(keys = "cycle.utility", mod = "*", values = function() get_attribute(bsc.sim, "Tx.time.cycle")) %>%
           set_attribute(keys = "Total.Utility", mod = "+", values = function() get_attribute(bsc.sim, "cycle.utility")) %>%
           set_attribute(key = "Alive", value = 0)                                                                     # update that the patient has died
  ) %>%
  
  # Follow up 2
  set_attribute(key = "Tx.time.fu2", value = function() t.fu2()) %>%   
  timeout_from_attribute(key = "Tx.time.fu2") %>%   
  set_attribute(keys = "Total.Utility", mod = "*", values = 1.1) %>%                                                              # stay in follow-up treatment for the determined time
  set_attribute(key = "Alive", value = 0)

# Visualize to check whether the defined model structure is ok
plot(bsc.model);


## Section 5: Simulation ----

# Simulation settings
set.seed(5678);       # random number seed for reproducibility
n.patients <- 10000;    # number of patients to simulate 
mon.patients <- 2;    # level of monitoring (see add_generator)

# Define simulation for the best standard care (bsc)
bsc.sim <- simmer() %>%
  add_resource(name="Tx1", capacity=Inf, mon=F) %>%
  add_resource(name="Tx2", capacity=Inf, mon=F) %>%
  add_generator(name_prefix="Patient", trajectory=bsc.model, distribution=at(rep(x=0, times=n.patients)), mon=mon.patients)

# Run the BSC simulation
bsc.sim %>% 
  run()

# Get the outcomes for the monitored attributes
bsc.out <- get_mon_attributes(bsc.sim);             # retrieve the monitor object
bsc.costs <- getMultipleAttributes(c("Total.Costs", "Total.Utility"), bsc.out);   # get outcomes for multiple outcomes at the same

