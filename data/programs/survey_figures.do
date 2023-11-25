clear all

set more off


*Set to your own directory
cd "C:\Users\l1sxf17\Dropbox (FRB SF)\Research\David\IGC_energy\replication_package\"
/*-----------------------------------
	Appending and saving Ghana files
-------------------------------------*/
* append and save all rounds of the main survey data from Ghana
use "./data/temp/ghana_main_survey_round1.dta",replace
append using "./data/temp/ghana_main_survey_round2.dta"
append using "./data/temp/ghana_main_survey_round3.dta"
append using "./data/temp/ghana_main_survey_round4.dta"

* Create a dummy variable "placebo" and give it a value of 0 (1 == placebo, 0==otherwise)
gen placebo = 0

save "./data/temp/ghana_main_survey_data.dta",replace

* append and save placebo data
use "./data/temp/ghana_placebo_round1.dta", clear
append using "./data/temp/ghana_placebo_round2.dta" 
append using "./data/temp/ghana_placebo_round3.dta" 
append using "./data/temp/ghana_placebo_round4.dta"

* Create a dummy variable "placebo" and give it a value of 1 (1 == placebo, 0==otherwise)
gen placebo = 1 

save "./data/temp/ghana_placebo_data.dta", replace

* Append Main and Placebo surveys

use "./data/temp/ghana_main_survey_data.dta", clear
append using "./data/temp/ghana_placebo_data.dta"

* Create string variable 'country' and give it the value "Ghana"
gen country = "Ghana"

save "./data/output/survey_data_ghana.dta", replace

/*-----------------------------------
	Appending and saving Nigeria files
-------------------------------------*/
* append and save all rounds of main survey data from Nigeria
use "./data/temp/nigeria_main_survey_round1.dta",replace
append using "./data/temp/nigeria_main_survey_round2.dta"
append using "./data/temp/nigeria_main_survey_round3.dta"
append using "./data/temp/nigeria_main_survey_round4.dta"

* Create a dummy variable "placebo" and give it a value of 0 (1 == placebo, 0==otherwise)
gen placebo = 0

save "./data/temp/nigeria_main_survey_data.dta",replace

* append and save all rounds of 'placebo' survey data from Nigeria
use "./data/temp/nigeria_placebo_round1.dta", clear
append using "./data/temp/nigeria_placebo_round2.dta" 
append using "./data/temp/nigeria_placebo_round3.dta" 
append using "./data/temp/nigeria_placebo_round4.dta"

* Create a dummy variable 'placebo' and give it a value of 1 (1 == placebo, 0==otherwise)
gen placebo = 1 

save "./data/temp/nigeria_placebo_data.dta", replace
	
* Append Main and Placebo surveys from Nigeria
use "./data/temp/nigeria_main_survey_data.dta", clear
append using "./data/temp/nigeria_placebo_data.dta"

* Create string variable 'country' and give it the value "Nigeria"
gen country = "Nigeria"

save "./data/output/survey_data_nigeria.dta", replace

* Now append the Ghana and Nigeria Survey files and save it as one file
use "./data/output/survey_data_ghana.dta", clear
append using "./data/output/survey_data_nigeria.dta"

save "./data/output/google_surveys.dta", replace

/*------------------------------------------------------------------------------
		This code reproduce Figure F.1: Freguency of power outages 
--------------------------------------------------------------------------------*/
use "./data/output/google_surveys.dta", clear
keep if placebo == 0 // keep only responses from main survey
gen yval1 = (outagefreq == 1) 
gen yval2 = (outagefreq == 2)
gen yval3 = (outagefreq == 3)
gen yval4 = (outagefreq == 4) 
gen yval5 = (outagefreq == 5) 
gen yval6 = (outagefreq == 6)
collapse (mean) yval*, by(country)
reshape long yval, i(country) j(xval)
replace yval = 100*yval
reshape wide yval, i(xval) j(country) string
graph bar yval*, over(xval, relabel(1 "Never" 2 "1 to 2 times" 3 "3 to 5 times" 4 "6 to 10" 5 `" "More than 10" "times" "' 6 "I don't know") label(labsize(small))) ///
bargap(3) blabel(total, format(%5.1f)) bar(1, fcolor(red) fint(100) lcolor(black)) bar(2, fcolor(gold) fint(100) lcolor(black)) graphregion(color(white)) ///
legend(order(1 "Ghana" 2 "Nigeria") symxsize(3) symysize(3) rows(1) colgap(30) region(lcolor(white))) name(outagesfreq, replace)
graph export ./figures/outagesfreq.eps, replace

/*------------------------------------------------------------------------------
 This code reproduce Figure 7: Effects of Eliminating Power Outages from the 
 Firms’ Perspective 
--------------------------------------------------------------------------------*/
use "./data/output/survey_data_ghana.dta", clear
append using "./data/output/survey_data_nigeria.dta"
foreach v of varlist raiseownprofit otherbusenter makenewinv hiremoreworkers{
bysort country: egen ms`v' = mean(`v'*100) if placebo == 0 // average results from main survey
bysort country: egen pt`v' = mean(`v'*100) if placebo == 1 // average results from placebo survey
bysort country: fillmissing ms`v' pt`v', with(any) 
gen pte`v' = . // placebo treatment effects
qui: reg `v' placebo if country == "Ghana"
replace pte`v' = -_b[placebo]*100 if country == "Ghana"
qui: reg `v' placebo if country == "Nigeria"
replace pte`v' = -_b[placebo]*100 if country == "Nigeria"
format ms`v' pt`v' pte`v' %5.1f
}
keep country ms* pt* pte*
duplicates drop country, force
reshape long ms pt pte, i(country) j(variable) string
rename (ms pt pte) (yval#), addnumber
reshape long yval, i(country variable) j(xvalue) 
replace xval = xval + 4 if variable == "makenewinv"
replace xval = xval + 8 if variable == "hiremoreworkers"
replace xval = xval + 12 if variable == "otherbusenter"
tostring yval, gen(yvaltodisplay) force format(%5.1f)
replace yvaltodisplay = yvaltodisplay + "***" if  (xval == 3 | xval == 7 |xval == 11 | xval == 15) & yvaltodisplay != "10.6"
replace yvaltodisplay = yvaltodisplay + "**" if xval == 11 & yvaltodisplay == "10.6"
sort country x
levelsof country, local(country)
foreach c of local country{
tw 	(bar yval xval if (xval == 1 | xval == 5 | xval == 9 | xval == 13) & country == "`c'", fcolor(red) fint(100) lcolor(black) barwidth(.99)) ///
	(bar yval xval if (xval == 2 | xval == 6 | xval == 10 | xval == 14)  & country == "`c'", fcolor(orange) fint(100) lcolor(black) barwidth(.99)) ///
	(bar yval xval if (xval == 3 | xval == 7 | xval == 11 | xval == 15)  & country == "`c'", fcolor(gold) fint(100) lcolor(black) barwidth(.99)) ///
	(scatteri 69 1 "Raise own profits", msymb(none) mlabcol(black) mlabsize(medsmall))(scatteri 69 5 "Expand own investment", msymb(none) mlabcol(black) mlabsize(medsmall)) ///
	(scatteri 69 9 "Expand own hiring", msym(none) mlabcol(black) mlabsize(medsmall))(scatteri 69 13 "Increase entry", msymb(none) mlabcol(black) mlabsize(medsmall)) ///  
	(scatter yval xval if country == "`c'", msymbol(none) mlab(yvaltodisplay) mlabpos(12) mlabcolor(black)), ///
ysc(titlegap(2)) yla(0(10)70, format(%9.0g) nogrid) yline(0, lcolor(black) lwidth(.2)) ytitle("Percent") xsc(off) xsize(2) ysize(1) graphregion(color(white)) ///
legend(rows(1) order(1 "Eliminate power outages" 2 "Placebo" 3 "Difference") symxsize(3) symysize(3) colgap(30) region(lcolor(white))) name(fig8_`c', replace)
graph export ./figures/fig8_`c'.eps, replace
}

//################################ Do file ENDs #################################
