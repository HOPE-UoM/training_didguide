// Workshop "Reading the "Diff-in-diff Practitioner's Guide" together" 
// University of Manchester, HOPE, 03/10/2025
// Author: Igor Francetic
// File: Step-by-step analysis

// Install two packages
ssc install csdid, replace
ssc install drdid, replace
ssc install estout, replace
ssc install bacondecomp, replace
ssc install coefplot, replace

// Set up the link to do file on GitHub
global simulation "https://raw.githubusercontent.com/HOPE-UoM/training_didguide/refs/heads/main/simulation.do"

// Set up the outcome we look at, we start with the dgp with conditional parallel trends
global outcome y_pt

///////////////////////////////////////////
// Step 1: Canonical 2x2 design without covariates
///////////////////////////////////////////

clear all
qui: do $simulation
xtset id time

// See percentage of treated and controls over time
tab time treated, row
tab treatment if time==1

// Let's start by only looking at controls and those treated in the first wave
// Let's also focus on time 3 (pre) and 4 (post)
keep if inlist(group, 1, 2)
keep if time>2&time<5
generate post=(time>3)
scalar att=10

// Method 1: Estimate DID from 4 averages
sum $outcome if post==1&treatment==1 // Post average for treated
scalar treatpost = `r(mean)'
sum $outcome if post==0&treatment==1 // Pre average for treated
scalar treatpre = `r(mean)'
sum $outcome if post==1&treatment==0 // Post average for controls
scalar contrpost = `r(mean)'
sum $outcome if post==0&treatment==0 // Pre average for controls
scalar contrpre = `r(mean)'
scalar did=(`=treatpost'-`=treatpre')-(`=contrpost'-`=contrpre')
scalar did1=abs((did-`=att')/`=att'*100)
* Knowing that the true treatment effect is 10, one can compute the bias
di "----------------------> This gives a " `=did1' "% bias"

// Method 2: Estimate DID with regression
reg $outcome c.treatment##c.post, cluster(id) 
scalar did2=abs((_b[c.treatment#c.post]-`=att')/`=att'*100)
di "----------------------> This gives a " `=did2' "% bias, the same as before"

///////////////////////////////////////////
// Step 2: Canonical 2x2 design with covariates
///////////////////////////////////////////

// Balance in covariates at baseline
// Note that staff_it is perfectly balanced but on casemix_it there is a small differences
global covs "casemix_it"
eststo treat: estpost summarize $covs if post==0&treatment==1
eststo contr: estpost summarize $covs if post==0&treatment==0
eststo diff: estpost ttest $covs if post==0, by(treatment)
estout contr treat diff, cells("mean(pattern(1 1 0) fmt(5)) sd(pattern(1 1 0)) b(star pattern(0 0 1) fmt(5)) t(pattern(0 0 1) par fmt(5))") label

// Check if normalised difference is more than rule of thumb value of 0.25 for casemix
forvalues i=0/1 {
sum $covs if post==0&treatment==`i'
scalar mean`i'=`r(mean)'
scalar sd`i'=`r(sd)'
}
scalar nd=(mean1-mean0)/((sd1^2+sd0^2)/2)^2
di `=nd'

// Estimate TWFE adjusting for time-varying covariates
xtreg $outcome c.treated c.post $covs, fe cluster(id)
scalar did3=abs((_b[c.treated]-`=att')/`=att'*100)
di "----------------------> This gives a " `=did3' "% bias, so we have not improved over " `=did2' "%!"

// Doubly robust estimation with conditional parallel trends recommended by The Guide
drdid $outcome $covs, t(post) treatment(treatment) dripw // Doubly robust
scalar did4=abs((`e(att1)'-`=att')/`=att'*100) 
di "----------------------> This gives a " `=did4' "% bias, some more improvement over " `=did3' "%!"
	   
///////////////////////////////////////////
// Step 3: DID with multiple time periods but only one treatment timing
///////////////////////////////////////////

clear all
qui: do $simulation
xtset id time

// Let's only keep those treated in year 4
keep if inlist(group, 1,2)
tab time treated, col
scalar att=10

// Method 1: Estimate DID from 4 averages for the first year of adoption
sum $outcome if time==4&treatment==1 // Post average for treated
scalar treatpost = `r(mean)'
sum $outcome if time==3&treatment==1 // Pre average for treated
scalar treatpre = `r(mean)'
sum $outcome if time==4&treatment==0 // Post average for controls
scalar contrpost = `r(mean)'
sum $outcome if time==3&treatment==0 // Pre average for controls
scalar contrpre = `r(mean)'
scalar did=(`=treatpost'-`=treatpre')-(`=contrpost'-`=contrpre')
scalar did5=abs((did-`=att')/`=att'*100)
* Knowing that the true treatment effect is 10, one can compute the bias
di "We get an estimate of " `=did', "which corresponds to a " `=did5' "% bias"

// Do we get the same reading the TWFE event study coefficient for the first year of implementation?
xtreg $outcome treatment##ib3.time, fe cluster(id)
coefplot, keep(1.treatment#*) vertical yline(0, lcol(red))
di "The answer is YES, of course. For the first adoption period the coefficient is again " `=_b[1.treatment#4.time]'

// Method 2: Estimate dynamic DID with regression adjustment and SEs adjusted for Multiple Hypothesis Testing, without covariates
csdid $outcome, time(time) gvar(treattime) method(dripw) event
scalar did=`=_b[t_3_4]'
scalar did6=abs((did-`=att')/`=att'*100)
estat event
csdid_plot
di "Note that we get exactly the same estimate for the first year of adoption, but confidence intervals are wider (MHT)"

// Estimate dynamic DID imposing conditional PT and compare aggregate estimates between CSDID and TWFE
csdid $outcome $covs, ivar(id) time(time) gvar(treattime) method(dripw) event 
estat event
csdid_plot
estat simple, estore(cs_att)
eststo twfe_att: xtreg $outcome treated i.time $covs, fe cluster(id)
esttab cs_att twfe_att, keep(ATT treated)

///////////////////////////////////////////
// Step 4: Staggered DID 
///////////////////////////////////////////

clear all
qui: do $simulation
tab time treated, col
xtset id time
scalar att = 10/3-10/3+3/3

di "Here the true overall ATT is " `=att'

// Wrong way to compute average effect in staggered settings using TWFE (with the problems identified by Goodman-Bacon 2021)
xtreg $outcome treated i.time, fe cluster(id) 
scalar did=`=_b[treated]'
scalar did9=abs((did-`=att')/`=att'*100)
di "The coefficient has a bias of " `=abs(did9)' "%"
bacondecomp $outcome treated, ddetail // This tells you the weight of different comparisons in the overall estimate

// Wrong way to compute average effect in staggered settings using TWFE with CPT (with the problems identified by Goodman-Bacon 2021)
xtreg $outcome treated i.time $covs, fe cluster(id)
scalar did=`=_b[treated]'
scalar did10=abs((did-`=att')/`=att'*100)
di "The coefficient has a bias of " `=abs(did10)' "%"

// Optimal method (Callaway and Sant'Anna  2021): Estimate staggered DID with DR estimator, conditional PT, and adjusting for Multiple Hypothesis Testing, and using never treated as control
csdid $outcome, ivar(id) time(time) gvar(treattime) method(dripw) event
estat group

// Optimal method (Callaway and Sant'Anna  2021): Estimate staggered DID with DR estimator, conditional PT, and adjusting for Multiple Hypothesis Testing, and using never treated as control
csdid $outcome $covs, ivar(id) time(time) gvar(treattime) method(dripw) event
estat group

// Optimal method (Callaway and Sant'Anna  2021): Estimate staggered DID with DR estimator, conditional PT, and adjusting for Multiple Hypothesis Testing, and using not yet treated as control
csdid $outcome $covs, ivar(id) time(time) gvar(treattime) method(dripw) event notyet
estat group