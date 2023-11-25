clear all
set more off


*Set to your own directory
cd "C:\Users\l1sxf17\Dropbox (FRB SF)\Research\David\IGC_energy\replication_package\"
/*---------------------------------
		Ghana: Main survey
-----------------------------------*/
/********************************************************************************
					Round 1
question1: Do you own and operate your own business?
question2: Approximately what is your business's average monthly _profits_ (revenues minus costs)?
question3: How often have you experienced power outages / dumsor at your business in the last year?
question4: Suppose that power outages were permanently eliminated in Ghana. How might this affect the average profits of your business in the coming years?
question5: Suppose you could pay 2,000 GHC to permanently eliminate power outages at your business. Would be willing to do so?
*********************************************************************************/

use "./data/raw/ghana_main_survery_rawdata1.dta", clear

destring responsetime* weight, replace

gen surveyround = 1

* Any outages??? 
gen anyoutage = .
replace anyoutage = 0 if question3answer == "Never"
replace anyoutage = 1 if question3answer == "1 to 2 times"
replace anyoutage = 1 if question3answer == "3 to 5 times"
replace anyoutage = 1 if question3answer == "6 to 10 times"
replace anyoutage = 1 if question3answer == "More than 10 times"

sum anyoutage

* Frequency of outages; 
gen outagefreq = .
replace outagefreq = 1 if question3answer == "Never"
replace outagefreq = 2 if question3answer == "1 to 2 times"
replace outagefreq = 3 if question3answer == "3 to 5 times"
replace outagefreq = 4 if question3answer == "6 to 10 times"
replace outagefreq = 5 if question3answer == "More than 10 times"
replace outagefreq = 6 if question3answer == "I don't know"

tab outagefreq

* Profit of firms; 
gen profit = .
replace profit = 1000 if question2answer == "Less than 1,000 GHC"
replace profit = 1500 if question2answer == "1,000 to 2,000 GHC"
replace profit = 3500 if question2answer == "2,000 to 5,000 GHC"
replace profit = 10000 if question2answer == "5,000 to 15,000 GHC"
replace profit = 30000 if question2answer == "15,000 to 45,000 GHC"
replace profit = 67500 if question2answer == "45,000 to 90,000 GHC"
replace profit = 90000 if question2answer == "More than 90,000 GHC"

*Likely profits will increase if power cuts are eliminated
gen raiseownprofit = .
replace raiseownprofit = 1 if question4answer == "Likely to raise my profits"
replace raiseownprofit = 0 if missing(raiseownprofit) 

sum raiseownprofit

save "./data/temp/ghana_main_survey_round1.dta",replace

/*******************************************************************************
								Round 2
question1: Do you own and operate your own business?
question2: How often have you experienced power outages / dumsor at your business in the last year?
question3: Suppose that power outages are permanently eliminated at your business. How likely would this be to raise the average profits of your business in the coming years?
question4: Suppose that power outages are permanently eliminated in Ghana. How likely would this be to raise the number of other businesses operating in your industry?
question5: Suppose you could pay 8,000 GHC to permanently eliminate power outages at your business. Would be willing to do so?
**********************************************************************************/

use "./data/raw/ghana_main_survery_rawdata2.dta", clear

destring responsetime* weight, replace

gen surveyround = 2 

* Any outages??? 
gen anyoutage = .
replace anyoutage = 0 if question2answer == "Never"
replace anyoutage = 1 if question2answer == "1 to 2 times"
replace anyoutage = 1 if question2answer == "3 to 5 times"
replace anyoutage = 1 if question2answer == "6 to 10 times"
replace anyoutage = 1 if question2answer == "More than 10 times"

* Frequency of outages;
gen outagefreq = .
replace outagefreq = 1 if question2answer == "Never"
replace outagefreq = 2 if question2answer == "1 to 2 times"
replace outagefreq = 3 if question2answer == "3 to 5 times"
replace outagefreq = 4 if question2answer == "6 to 10 times"
replace outagefreq = 5 if question2answer == "More than 10 times"
replace outagefreq = 6 if question2answer == "I don't know"

* Very likely/Likely profits will increase if power cuts are eliminated???
gen raiseownprofit = .
replace raiseownprofit = 1 if question3answer == "Likely"
replace raiseownprofit = 1 if question3answer == "Very likely"
replace raiseownprofit = 0 if missing(raiseownprofit) 


* Other business will enter industry if power cuts are eliminated???
gen otherbusenter = .
replace otherbusenter = 1 if question4answer == "Likely"
replace otherbusenter = 1 if question4answer == "Very likely"
replace otherbusenter = 0 if missing(otherbusenter) 

save "./data/temp/ghana_main_survey_round2.dta",replace

/*******************************************************************************
							Round 3
question1: Do you own and operate your own business?
question2: How often have you experienced power outages / dumsor at your business in the last year?
question3: Suppose that power outages are permanently eliminated at your business. How likely would it be that you respond by making new investments in your business?
question4: Suppose that power outages are permanently eliminated at your business. How likely would it be that you respond by hiring more workers?
question5: Suppose you could pay 16,000 GHC to permanently eliminate power outages at your business. Would be willing to do so?
*********************************************************************************/

use "./data/raw/ghana_main_survery_rawdata3.dta", clear

destring responsetime* weight, replace

gen surveyround = 3

* Any outages??? 
gen anyoutage = .
replace anyoutage = 0 if question2answer == "Never"
replace anyoutage = 1 if question2answer == "1 to 2 times"
replace anyoutage = 1 if question2answer == "3 to 5 times"
replace anyoutage = 1 if question2answer == "6 to 10 times"
replace anyoutage = 1 if question2answer == "More than 10 times"

* Frequency of outages; 
gen outagefreq = .
replace outagefreq = 1 if question2answer == "Never"
replace outagefreq = 2 if question2answer == "1 to 2 times"
replace outagefreq = 3 if question2answer == "3 to 5 times"
replace outagefreq = 4 if question2answer == "6 to 10 times"
replace outagefreq = 5 if question2answer == "More than 10 times"
replace outagefreq = 6 if question2answer == "I don't know"

* Very likely/Likely firm will make new investments if power cuts are eliminated???
gen makenewinv = .
replace makenewinv = 1 if question3answer == "Likely"
replace makenewinv = 1 if question3answer == "Very likely"
replace makenewinv = 0 if missing(makenewinv) 

* Hire more workers if power cuts are eliminated???
gen hiremoreworkers = .
replace hiremoreworkers = 1 if question4answer == "Likely"
replace hiremoreworkers = 1 if question4answer == "Very likely"
replace hiremoreworkers = 0 if missing(hiremoreworkers) 

save "./data/temp/ghana_main_survey_round3.dta",replace

/*******************************************************************************
							Round 4
question1: Do you own and operate your own business?
question2: How long has your business been in operation?
question3: How often have you experienced power outages / dumsor at your business in the last year?
question4: Suppose that power outages are permanently eliminated at your business. How likely would it be that you respond by making new investments in your business?
question5: Suppose you could pay 16,000 GHC to permanently eliminate power outages at your business. Would be willing to do so?
********************************************************************************/

use "./data/raw/ghana_main_survery_rawdata4.dta", clear

destring responsetime* weight, replace

gen surveyround = 4
* Any outages??? 
gen anyoutage = .
replace anyoutage = 0 if question3answer == "Never"
replace anyoutage = 1 if question3answer == "1 to 2 times"
replace anyoutage = 1 if question3answer == "3 to 5 times"
replace anyoutage = 1 if question3answer == "6 to 10 times"
replace anyoutage = 1 if question3answer == "More than 10 times"

* Frequency of outages; 
gen outagefreq = .
replace outagefreq = 1 if question3answer == "Never"
replace outagefreq = 2 if question3answer == "1 to 2 times"
replace outagefreq = 3 if question3answer == "3 to 5 times"
replace outagefreq = 4 if question3answer == "6 to 10 times"
replace outagefreq = 5 if question3answer == "More than 10 times"
replace outagefreq = 6 if question3answer == "I don't know"

* Very likely/Likely firm will make new investments if power cuts are eliminated???
gen makenewinv = .
replace makenewinv = 1 if question4answer == "Likely"
replace makenewinv = 1 if question4answer == "Very likely"
replace makenewinv = 0 if missing(makenewinv) 

save "./data/temp/ghana_main_survey_round4.dta",replace

/*----------------------------------------
	Ghana: Placebo survey
------------------------------------------*/

* Round 1

use "./data/raw/ghana_placebo_rawdata1.dta", clear

destring responsetime* weight, replace

gen surveyround = 1

* Airport converting to solar will raise own profit
gen raiseownprofit = .
replace raiseownprofit = 1 if question3answer == "Likely" 
replace raiseownprofit = 1 if question3answer == "Very likely" 
replace raiseownprofit = 0 if missing(raiseownprofit)  

* Airport converting to solar will encourage new entrants
gen otherbusenter = .
replace otherbusenter = 1 if question4answer == "Likely" 
replace otherbusenter = 1 if question4answer == "Very likely" 
replace otherbusenter = 0 if missing(otherbusenter) 

save "./data/temp/ghana_placebo_round1.dta", replace

* Round 2
use "./data/raw/ghana_placebo_rawdata2", clear

destring responsetime* weight, replace

gen surveyround = 2

* Airport converting to solar will make firm increase investments
gen makenewinv = .
replace makenewinv = 1 if question3answer == "Likely"
replace makenewinv = 1 if question3answer == "Very likely"
replace makenewinv = 0 if missing(makenewinv) 

* Airport converting to solar will make firm employ more
gen hiremoreworkers = .
replace hiremoreworkers = 1 if question4answer == "Likely"
replace hiremoreworkers = 1 if question4answer == "Very likely"
replace hiremoreworkers = 0 if missing(hiremoreworkers) 

save "./data/temp/ghana_placebo_round2.dta", replace

* Round 3
use "./data/raw/ghana_placebo_rawdata3", clear

destring responsetime* weight, replace

gen  surveyround = 3

* Airport converting to solar will raise own profit
gen raiseownprofit = .
replace raiseownprofit = 1 if question3answer == "Likely" 
replace raiseownprofit = 1 if question3answer == "Very likely" 
replace raiseownprofit = 0 if missing(raiseownprofit)  

* Airport converting to solar will make firm employ more
gen hiremoreworkers = .
replace hiremoreworkers = 1 if question4answer == "Likely"
replace hiremoreworkers = 1 if question4answer == "Very likely"
replace hiremoreworkers = 0 if missing(hiremoreworkers) 

save "./data/temp/ghana_placebo_round3.dta", replace

/*-------------
	Round 4
---------------*/
use "./data/raw/ghana_placebo_rawdata4", clear

destring responsetime* weight, replace

gen  surveyround = 4

* Airport converting to solar will make firm increase investments
gen makenewinv = .
replace makenewinv = 1 if question3answer == "Likely"
replace makenewinv = 1 if question3answer == "Very likely"
replace makenewinv = 0 if missing(makenewinv) 

* Airport converting to solar will make firm employ more
gen hiremoreworkers = .
replace hiremoreworkers = 1 if question4answer == "Likely"
replace hiremoreworkers = 1 if question4answer == "Very likely"
replace hiremoreworkers = 0 if missing(hiremoreworkers) 

save "./data/temp/ghana_placebo_round4.dta", replace
//################################ Do file ENDs #################################
