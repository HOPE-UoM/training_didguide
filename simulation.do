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
bys id: g time=_n
sort id time
xtset id time

// 3. Create four treatment cohorts with one never treated group
gen treatment=(id>1999)
egen countertreat=group(id) if treatment==1
qui: sum countertreat
local TN `r(max)'
gen treattime = .
replace treattime = 0 if treatment==0          												// Never Treated (Control)
replace treattime = 4 if treatment==1&countertreat <= (`TN'/3)  							// Treated at t=4
replace treattime = 6 if treatment==1&countertreat > (`TN'/3) & countertreat <= (2*`TN'/3) 	// Treated at t=6
replace treattime = 8 if treatment==1&countertreat > (2*`TN'/3)                 			// Treated at t=8
egen group=group(treattime)
label define grouplabs 1 "Never treated" 2 "Treated at t=4" 3 "Treated at t=6" 4 "Treated at t=8"
label values group grouplabs
gen timetotreat= time - treattime
replace timetotreat = -100 if treattime == 0 // Tgroup gets a distinct code
gen treated = (timetotreat >= 0) if timetotreat != -100
replace treated = 0 if treattime == 0

// 4. Create unobservable time-varying confounding (correlated with y and treatment)
generate effort_it = rnormal(50, 2)+5*treatment if time==1
replace effort_it = l.effort_it+0.1*time*treatment if missing(effort_it)

// 5. Generate practice-level unobserved ability (a_i), observable confounders and noise (e.g. case-mix, staff) and idiosyncratic error (epsilon_it)
gen casemix_it=rnormal(10, 20)
gen a_i = rnormal(0, 5) + treatment*5 if time == 1
by id: replace a_i = a_i[1]
gen epsilon_it = rnormal(0, 10) 

// 6. Simulate outcome with heterogeneous (true) treatment effects across treatment cohorts and different selection mechanisms

// Define true treatment effect for the different cohorts
gen te=.
replace te=0 if group==1
replace te=10 if group==2 
replace te=-10 if group==3 
replace te=3 if group==4 

// The true model with selection on observables: Y_it = beta_0 + tau*time + beta_casemix*casemix_it + te * treated + staff_it + a_i + epsilon_it
gen y_pt = 10 + 0.3*time + te * treated + casemix_it + a_i + epsilon_it
gen y_npt = 10 + 0.3*time + te * treated + casemix_it + 0.5 * effort_it + a_i + epsilon_it

keep y* id group time casemix_it treated treatment treattime timetotreat
order y* id group time casemix_it treated treatment treattime timetotreat
