**********************************************************************
***** Program: moments.do
***** By: Stephie Fried and David Lagakos
***** Date: Winter 2022
***** Purpose: Computes targets from the WES data for model calibration
**********************************************************************

clear
set type double

**Set own directory
*cd "C:/Users/sdfried/Dropbox (ASU)/Research/David/IGC_energy/replication_package"

**Set local options
local countries = "Nigeria Ghana Tanzania Uganda" // Set countries to run for.
			 // Valid options include: Nigeria, Ghana, Tanzania, and/or Uganda

foreach country in `countries' {

if "`country'" == "Nigeria" {
	loc year = "2014"
}
else {
	loc year = "2013"
}
use "./data/raw/`country'-`year'-full-data.dta", clear

*Create indicator == 1 if interviewee has generator
gen has_gen = 1 if c10 == 1
replace has_gen = 0 if c10 == 2
label variable has_gen "indicator =1 if has generator, zero otherwise"

*Create variable for fraction of own electricity from generator
gen self_gen_frac = c11 if c11>0
label variable self_gen_frac "Fraction of electricity firms with generators generate themselves"

*Create log size variable
gen size = l1 if l1>0
gen ln_size = log(size)
label variable size "Firm size (number of employees)"
label variable ln_size "Log of firm size"

*Drop if size < 10
drop if size < 10

*Create variable for fraction of no outages
*C6: Over last FY, did this establishment experience outages
gen has_outages = 1 if c6 ==1 //experienced outages
replace has_outages =0 if c6 ==2  // did not experience outages
label variable has_outages "indicator =1 if has outages, zero otherwise"

egen total_firms = total(wstrict) if ~missing(has_outages)
egen no_outages_firms = total(wstrict) if has_outages == 0
label variable total_firms "total number of firms"
label variable no_outages_firms "Number of firms that don't experience outages" 

gen frac_no_outages = no_outages_firms/total_firms
label variable frac_no_outages "Fraction of firms that don't experience outages"

*Get means of variables
sum has_gen [aw = wstrict] if has_outages ==1 //fraction of modern firms that experience outages that have a generator
gen has_gen_mean = r(mean)
label variable has_gen_mean "Average fraction of modern firms that experience outages that have generators"
sum self_gen_frac [aw=wstrict] if has_gen  ==1 //fraction of electricity that firms with a generator generate themselves
gen self_gen_frac_mean = r(mean)
label variable self_gen_frac_mean "Average fraction of of electricity that firms with generators generate themselves"

*Linear probability model
reg has_gen ln_size [aw=wstrict], robust

*Export results to spreadsheet
matrix results =  r(table)
svmat results

scalar est = results[1,1]
scalar uci = results[6,1]
scalar lci = results[5,1]

gen beta = est
gen upperb = uci
gen lowerb = lci
label variable beta "Estimated semi-elasticity of generator ownership with respect to firm size"
label variable upperb "Upper bound on the 95 percent confidence interval"
label variable lowerb "Lower bound on the 95 percent confidence interval"

collapse (max) beta upperb lowerb has_gen_mean self_gen_frac_mean frac_no_outages

*export has_gen, self_gen_frac, and frac_no_outages to spreadsheet
export excel has_gen_mean self_gen_frac_mean frac_no_outages using "./data/output/results_micro.xlsx", ///
sheet(variables_`country') firstrow(variables) sheetreplace keepcellfmt

*export regression results
export excel beta upperb lowerb using "./data/output/results_micro.xlsx", ///
sheet(estimates_`country') firstrow(variables) sheetreplace keepcellfmt

}

