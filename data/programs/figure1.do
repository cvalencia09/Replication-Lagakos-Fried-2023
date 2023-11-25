**********************************************************************
***** Program: figure1.do
***** By: Stephie Fried and David Lagakso
***** Date: Winter 2022
***** Purpose: Process data for figure 1
**********************************************************************

clear
set type double
cd "C:/Users/l1sxf17/Dropbox (FRB SF)/Research/David/IGC_energy/replication_package"

import excel "./data/raw/ExportedResults.xlsx", firstrow clear sheet("raw") 

rename Percentoffirmsexperiencingel outages
rename Economy country
keep outages country 
save ./data/temp/temp.dta, replace

import excel "./data/derived/country_codes.xlsx", firstrow clear
rename CountryName country
merge 1:1 country using ./data/temp/temp.dta
keep if _merge ==3
drop _merge

*Year for GDP per capita plot
gen year =2013

*Add 2013 GDP data from PWT
merge 1:1 year countrycode using ./data/raw/pwt91.dta
keep if _merge ==3
drop _merge
gen gdp_pop = rgdpe/pop
label variable gdp_pop "GDP per capita  (millions of dollars per million people)"


keep year gdp_pop outages country countrycode pop
drop if pop <= 5 //drop countries with less than 5 million people
gen gdp_pop_1000 = gdp_pop/1000
label variable gdp_pop_1000 "GDP per capita (billions of dollars per million people)"


destring(outages), replace force
drop if missing(outages)

export excel countrycode gdp_pop_1000 outages using "./data/output/fig1_data.xlsx", firstrow(variables) replace
 