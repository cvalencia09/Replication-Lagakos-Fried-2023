**********************************************************************
***** Program: platts.do
***** By: Stephie Fried, David Lagakos
***** Date: Winter 2019
***** Purpose: Calculate GDP per capita relative to Nigeria
***** Files created: 
***** Files used: allunits.dta
**********************************************************************
clear
set type double
cd "C:/Users/sdfried/Dropbox (ASU)/Research/David/IGC_energy/Data/calibration_submission/stata"


clear
use ./raw/pwt91.dta

keep if inlist(countrycode, "ETH", "GHA", "NGA", "TZA", "UGA") 
keep countrycode pop year cgdpe
gen gdp_pop = cgdpe/pop

sort year countrycode
gen gdp_pop_NGA = gdp_pop if countrycode =="NGA"
replace gdp_pop_NGA = gdp_pop_NGA[_n+1] if countrycode =="GHA"
replace gdp_pop_NGA = gdp_pop_NGA[_n+1] if countrycode =="ETH"
replace gdp_pop_NGA = gdp_pop_NGA[_n-1] if countrycode =="TZA"
replace gdp_pop_NGA = gdp_pop_NGA[_n-1] if countrycode =="UGA"

gen hello = gdp_pop/gdp_pop_NGA

sort country year
gen growth = gdp_pop[_n+1]/gdp_pop[_n] if year < 2014


gen gdp_pop_ref = gdp_pop if countrycode =="NGA" & year == 2014
sort gdp_pop_ref
replace gdp_pop_ref = gdp_pop_ref[_n-1] if missing(gdp_pop_ref)
gen rel_gdp_pop = gdp_pop/gdp_pop_ref

list rel_gdp_pop if countrycode =="ETH" & year ==2015
list rel_gdp_pop if countrycode =="GHA" & year ==2013
list rel_gdp_pop if countrycode =="TZA" & year ==2013
list rel_gdp_pop if countrycode =="UGA" & year ==2013


























