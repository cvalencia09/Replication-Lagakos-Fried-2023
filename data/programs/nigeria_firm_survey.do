clear all

set more off


*Set to your own directory
cd "C:\Users\l1sxf17\Dropbox (FRB SF)\Research\David\IGC_energy\replication_package\"
/*---------------------------------
		Nigeria: main survey
-----------------------------------*/
/*******************************************************************************
									Round 1
question1: Do you own and operate your own business?
question2: How often have you experienced power outages at your business in the last year?
question3: Suppose that power outages were permanently eliminated in Nigeria. How likely would this be to raise the average profits of your business in the coming years?
question4: Suppose that power outages are permanently eliminated across Nigeria. How likely would this be to raise the number of other businesses operating in your industry?
question5: Suppose that you could pay N150,000 to permanently eliminate power outages at your business. Would be likely to do so?
********************************************************************************/

use "./data/raw/nigeria_main_survey_rawdata1.dta", clear

destring responsetime* weight, replace

gen  surveyround = 1

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

*Likely profits will increase if power cuts are eliminated
gen raiseownprofit = .
replace raiseownprofit = 1 if question3answer == "Likely"
replace raiseownprofit = 1 if question3answer == "Very likely"
replace raiseownprofit = 0 if missing(raiseownprofit) 

* Other business will enter industry if power cuts are eliminated???
gen otherbusenter = .
replace otherbusenter = 1 if question4answer == "Likely"
replace otherbusenter = 1 if question4answer == "Very likely"
replace otherbusenter = 0 if missing(otherbusenter) 

save "./data/temp/nigeria_main_survey_round1.dta",replace

/*******************************************************************************
									Round 2
question1: Do you own and operate your own business?
question2: How often have you experienced power outages at your business in the last year?
question3: Suppose that power outages were permanently eliminated at your business. How likely would this be to raise the average profits of your business in the coming years??
question4: Suppose that power outages are permanently eliminated across Nigeria. How likely would this be to raise the number of other businesses operating in your industry?
question5: Suppose that you could pay N600,000 to permanently eliminate power outages at your business. Would be likely to do so?
********************************************************************************/

use "./data/raw/nigeria_main_survey_rawdata2.dta", clear

destring responsetime* weight, replace

gen  surveyround = 2

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

save "./data/temp/nigeria_main_survey_round2.dta",replace

/*******************************************************************************
								Round 3
question1: Do you own and operate your own business?
question2: How often have you experienced power outages at your business in the last year?
question3: Suppose that power outages are permanently eliminated at your business. How likely would it be that you respond by making new investments in your business?
question4: Suppose that power outages are permanently eliminated at your business. How likely would it be that you respond by hiring more workers?
question5: Suppose you could pay N1,200,000 to permanently eliminate power outages at your business. Would be willing to do so?
********************************************************************************/

use "./data/raw/nigeria_main_survey_rawdata3.dta", clear

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

save "./data/temp/nigeria_main_survey_round3.dta",replace

/*******************************************************************************
							Round 4
question1: Do you own and operate your own business?
question2: How often have you experienced power outages at your business in the last year?
question3: Suppose that power outages are permanently eliminated at your business. How likely would it be that you respond by making new investments in your business?
question4: Suppose that power outages are permanently eliminated at your business. How likely would it be that you respond by hiring more workers?
question5: Suppose you could pay N2,400,000 to permanently eliminate power outages at your business. Would be willing to do so?
********************************************************************************/

use "./data/raw/nigeria_main_survey_rawdata4.dta", clear

destring responsetime* weight, replace

gen surveyround = 4

* Any outages??? 
gen anyoutage = .
replace anyoutage = 0 if question2answer == "Never"
replace anyoutage = 1 if question2answer == "1 to 2 times"
replace anyoutage = 1 if question2answer == "3 to 5 times"
replace anyoutage = 1 if question2answer == "6 to 10 times"
replace anyoutage = 1 if question2answer == "More than 10 times"

* Frequency of outages
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

save "./data/temp/nigeria_main_survey_round4.dta",replace

/*----------------------------------------
	Nigeria: Placebo survey
------------------------------------------*/

/******************************************************************************* 
								Round 1
question1: Do you own and operate your own business?
question2: Nigeria’s airports are considering converting to solar power. In your opinion, is this a good idea?
question3: Suppose that Nigeria's airports convert to solar power. How likely would it be that you respond by making new investments in your business?
question4: Suppose that Nigeria's airports convert to solar power. How likely would it be that you respond by hiring more workers?
question5: Suppose that you could pay an additional tax of N1,200,000 to convert Nigeria's airports to solar power? Would you be willing to do so?
********************************************************************************/

use "./data/raw/nigeria_placebo_rawdata1.dta", clear

destring responsetime* weight, replace

gen  surveyround = 1

* Very likely/likely airport converting to solar will make firm make new investments???
gen makenewinv = .
replace makenewinv = 1 if question3answer == "Likely" 
replace makenewinv = 1 if question3answer == "Very likely" 
replace makenewinv = 0 if missing(makenewinv)  

* Very likely/likely airport converting to solar will make firm hire more workers?
gen hiremoreworkers = .
replace hiremoreworkers = 1 if question4answer == "Likely" 
replace hiremoreworkers = 1 if question4answer == "Very likely" 
replace hiremoreworkers = 0 if missing(hiremoreworkers) 

save "./data/temp/nigeria_placebo_round1.dta", replace

/******************************************************************************* 
									Round 2
question1: Do you own and operate your own business?
question2: How long has your business been in operation?
question3: Nigeria’s airports are considering converting to solar power. In your opinion, is this a good idea?
question4: Suppose that Nigeria's airports convert to solar power. How likely would it be that you respond by making new investments in your own business?
question5: Suppose that you could pay an additional tax of N2,400,000 to convert Nigeria's airports to solar power? Would you be willing to do so?
********************************************************************************/

use "./data/raw/nigeria_placebo_rawdata2.dta", clear

destring responsetime* weight, replace

gen surveyround = 2

* Very likely/likely airport converting to solar will make firm increase investments???
gen makenewinv = .
replace makenewinv = 1 if question4answer == "Likely"
replace makenewinv = 1 if question4answer == "Very likely"
replace makenewinv = 0 if missing(makenewinv) 

save "./data/temp/nigeria_placebo_round2.dta", replace

/*******************************************************************************
									Round 3
question1: Do you own and operate your own business?
question2: Nigeria’s airports are considering converting to solar power. In your opinion, is this a good idea?
question3: Suppose that Nigeria's airports convert to solar power. How likely would it be that profits at your own business rise as a result?
question4: Suppose that Nigeria's airports convert to solar power. How likely would it be that the number of businesses operating in your industry rises as a result?
question5: Suppose that you could pay N150,000 to help convert Nigeria's airports to solar power. Would you be willing to do so?
********************************************************************************/					

use "./data/raw/nigeria_placebo_rawdata3.dta"

destring responsetime* weight, replace

gen surveyround = 3

* Very likely/likely airport converting to solar will raise own profit???
gen raiseownprofit = .
replace raiseownprofit = 1 if question3answer == "Likely" 
replace raiseownprofit = 1 if question3answer == "Very likely" 
replace raiseownprofit = 0 if missing(raiseownprofit)  

* Very likely/likely airport converting to solar will make other businesses enter industry???
gen otherbusenter = .
replace otherbusenter = 1 if question4answer == "Likely"
replace otherbusenter = 1 if question4answer == "Very likely"
replace otherbusenter = 0 if missing(otherbusenter) 

save "./data/temp/nigeria_placebo_round3.dta", replace

/*******************************************************************************
								Round 4
question1: Do you own and operate your own business?
question2: How long has your business been in operation?
question3: Nigeria’s airports are considering converting to solar power. In your opinion, is this a good idea?
question4: Suppose that Nigeria's airports convert to solar power. How likely would it be that you respond by hiring new workers at your business?
question5: Suppose that you could pay an additional tax of N600,000 to help convert Nigeria's airports to solar power. Would you be willing to do so?
********************************************************************************/

use "./data/raw/nigeria_placebo_rawdata4.dta", clear

destring responsetime* weight, replace

gen surveyround = 4

* Very likely/likely airport converting to solar will make firm employ more
gen hiremoreworkers = .
replace hiremoreworkers = 1 if question4answer == "Likely"
replace hiremoreworkers = 1 if question4answer == "Very likely"
replace hiremoreworkers = 0 if missing(hiremoreworkers) 

save "./data/temp/nigeria_placebo_round4.dta", replace
