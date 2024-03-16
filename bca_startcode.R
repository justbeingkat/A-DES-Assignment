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

# Set the working directory
setwd("C:/path/to/your/DES files/");

# Load funtions for extracting monitored attributes
source("getSingleAttribute.R", echo=T);
source("getMultipleAttributes.R", echo=T);


## Section 2: Data analysis ----

# Load the dataset
load("trial_dataset.RData");

# Define parameters
# Patient characteristics 
p.male <- 0.33;                               # probability of being male
p.poor <- 0.20;                               # probability of a poor condition
m.age.men <- 57.25063;                        # mean age for men
sd.age.men <- 8.05432;                        # sd age for men
m.age.women <- 65.63482;                      # mean age for women
sd.age.women <- 13.27633;                     # sd age for women
# Tx general
t.cycle <- 30;    # average time of a normal treatment cycle
t.major <- 6;     # average time of a cycle in which major complications occur
t.death <- 15;    # average time of a cycle in which the patient dies

p.minor <- 0.10;  # probability of minor complications in a cycle
p.major <- 0.04;  # probability of major complications in a cycle
p.death <- 0.03;  # probability of death in a cycle

# Tx1 specific
c.Tx1.cycle <- 375;   # costs of a cycle of Tx1
c.Tx1.day <- 8;       # additional dayly costs when on treatment Tx1
u.Tx1 <- 0.55/365;    # utility per day when on treatment Tx1

p.Tx1.poor <- mean(data$Tx1.C1.Dx.Pet[data$Poor==1]==1, na.rm=T);                               # probability of effective Tx1 treatment when in poor condition
p.Tx1.good <- mean(data$Tx1.C1.Dx.Pet[data$Poor==0]==1, na.rm=T);                               # probability of effective Tx1 treatment when in good condition

# Tx2 specific
c.Tx2.cycle <- 4000;  # costs of a cycle of Tx2
c.Tx2.day <- 15;      # additional dayly costs when on treatment Tx2
u.Tx2 <- 0.5/365;     # utility per day when on treatment Tx1

# p.Tx2.yes.exp and p.Tx2.no.exp are provided for step 2.4
p.Tx2.yes.exp <- seq(from=0.39, to=0.47, by=0.02);                                               # probability of effective Tx2 treatment when responded to Tx1, dependent on cycles Tx1 
p.Tx2.no.exp <- seq(from=0.87, to=0.31, by=-0.14);                                               # probability of effective Tx2 treatment when not responded to Tx1, dependent on cycles Tx1 

# FU1 en FU2 specific
t.fu1.full <- 63          # average time spent in the first follow up if the patient survives during follow up
t.fu1.death <- 42         # average time spent in the first follow up if the patient dies during follow up
t.fu2 <- 100              # average time spent in the second follow up after Tx2
p.death.followup <- 0.05; # probability of dying during first follow up

## Section 3: Supportive functions ----

# Function for determining the event to happen
Tx1.event <- function() {
  
  #Randomly select whether the patient dies with a 10% probability or not
  event <- ifelse(runif(1) < 0.10, 2, 1);
  
  return(event);                                                                                                  # A return value equal to 0 skips the branch and continues to the next activity.
  
} # Function for defining the event during a cycle of Tx1

# Functions for determining the time-to-events
Tx1.time <- function(Tx1.Event) {
  
  return(30);
  
} # Function for defining the time spent on a cycle of Tx1 




## Section 4: Discrete event simulation model ----

# Define the model structure for the current practice, i.e. best standard care (BSC)
bsc.model <- trajectory() %>%
  
  # Initialization
  set_attribute(key="Alive", value=1) %>%                                                                          # define an attribute to check whether the patient is alive
  
  # First-line treatment
  set_attribute(key="Tx1.Event", value=function() Tx1.event()) %>%                                                 # select the event to happen in this treatment cycle          
  branch(option=function() get_attribute(bsc.sim, "Tx1.Event"), continue=c(T, F),
         
         # Event 1: Full cycle
         trajectory() %>%
           set_attribute(key="Tx1.Time", value=function() Tx1.time(get_attribute(bsc.sim, "Tx1.Event"))) %>%       # determine how long the cycle will last
           seize(resource="Tx1", amount=1) %>%                                                                     # occupy a place in first-line treatment
           timeout_from_attribute(key="Tx1.Time") %>%                                                              # stay in first-line treatment for the determined time
           release(resource="Tx1", amount=1) %>%                                                                   # leave first-line treatment
           rollback(amount=6, times=Inf),                                                                          # go back for another cycle (Hint: look at plot trajectory)
         
         # Event 2: Death
         trajectory() %>%
           set_attribute(key="Tx1.Time", value=function() Tx1.time(get_attribute(bsc.sim, "Tx1.Event"))) %>%       # determine how long the cycle will last
           seize(resource="Tx1", amount=1) %>%                                                                     # occupy a place in first-line treatment
           timeout_from_attribute(key="Tx1.Time") %>%                                                              # stay in first-line treatment for the determined time
           release(resource="Tx1", amount=1) %>%                                                                   # leave first-line treatment
           set_attribute(key="Alive", value=0)                                                                     # update that the patient has died
           
         
  ) # branch first-line treatment

# Visualize to check whether the defined model structure is ok
plot(bsc.model);
  




## Section 5: Simulation ----

# Simulation settings
set.seed(5678);       # random number seed for reproducibility
n.patients <- 100;    # number of patients to simulate 
mon.patients <- 2;    # level of monitoring (see add_generator)

# Define simulation for the best standard care (bsc)
bsc.sim <- simmer() %>%
  add_resource(name="Tx1", capacity=Inf, mon=F) %>%
  add_generator(name_prefix="Patient", trajectory=bsc.model, distribution=at(rep(x=0, times=n.patients)), mon=mon.patients)

# Run the BSC simulation
bsc.sim %>% 
  run()

# Get the outcomes for the monitored attributes
bsc.out <- get_mon_attributes(bsc.sim);             # retrieve the monitor object
getSingleAttribute("Alive", bsc.out);               # get patient-level outcomes for the attribute of interest
View(getMultipleAttributes(c("Alive"), bsc.out));   # get outcomes for multiple outcomes at the same time


