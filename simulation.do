// Workshop "Reading the "Diff-in-diff Practitioner's Guide" together" 
// University of Manchester, HOPE, 03/10/2025
// Author: Igor Francetic
// File: Simulating the dataset

// 1. Set up environment and parameteres
clear all
local N 6000 // Number of GP practices (units)
local T 10 // Number of time periods (t=1 to 5)
set obs 60000 // Number of observations (6000 units * 10 periods)
set seed 333 // for reproducibility

// 2. Generate panel structure
gen id = ceil(_n/`T')
gen time = mod(_n-1, `T') + 1
sort id time
xtset id time

// 3. Create confounding variable x (correlated with y and treatment)
gen x_it = rnormal(6, 6) + runiform(0,1)*time
gen ps=invlogit(-3 + 0.9 * x_it + rnormal(0, 1))
gen treatment=(ps>0.3) if time==1
bysort id: replace treatment = treatment[1] // Treatment is a time-invariant group indicator

// 4. Generate practice-level unobserved ability (a_i), observable noise (e.g. case-mix, n_it) and idiosyncratic error (epsilon_it)
gen a_i = rnormal(0, 5) if time == 1
by id: replace a_i = a_i[1]
gen epsilon_it = rnormal(0, 10) 
gen n_it = rnormal(3,3)

// 5. Create four treatment cohorts with one never treated group

// Create a variable to define the cohort based on 'id' with different treatment times
egen countertreat=group(id) if treatment==1
qui: sum countertreat
local TN `r(max)'
gen treattime = .
replace treattime = 0 if treatment==0          												// Never Treated (Control)
replace treattime = 4 if treatment==1&countertreat <= (`TN'/3)  							// Treated at t=4
replace treattime = 6 if treatment==1&countertreat > (`TN'/3) & countertreat <= (2*`TN'/3) 	// Treated at t=6
replace treattime = 9 if treatment==1&countertreat > (2*`TN'/3)                 			// Treated at t=9
egen group=group(treattime)
label define grouplabs 1 "Never treated" 2 "Treated at t=4" 3 "Treated at t=6" 4 "Treated at t=9"
label values group grouplabs

// Create relative time-to-treatment (timetotreat) & instantaneous treatment dummy (treated)
gen timetotreat= time - treattime
replace timetotreat = -100 if treattime == 0 // Tgroup gets a distinct code
gen treated = (timetotreat >= 0) if timetotreat != -100
replace treated = 0 if treattime == 0

// 6. Simulate outcome with heterogeneous (true) treatment effect

// Define true treatment effect for the different cohorts§
gen te=.
replace te=0 if group==1
replace te=10 if group==2 
replace te=-5 if group==3 
replace te=1 if group==4 

local beta_x 0.3 // True effect of X on Y
local tau=0.3 // Define the true learning effect over time

// True model: Y_it = beta_0 + tau*time + beta_x*x + te * treated + n_it + a_i + epsilon_it
gen y = 10 + `tau'*time + `beta_x' * x_it + te * treated + n_it + a_i + epsilon_it