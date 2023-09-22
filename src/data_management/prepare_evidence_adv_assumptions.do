/* 

This file computes average price and and average advertising data to be used in our reduced form regressions to support our assumptions.

*/


* BASICS.
clear
capture log close
version 13
set more off
set mem 4g
eststo clear
* Header do-File with path definitions, those end up in global macros.
include project_paths
log using `"${PATH_OUT_DATA}/log/`1'.log"', replace

*********************************
* Load survey data 
use `"${PATH_OUT_DATA}/new_master_data.dta"', clear

* Create average retail prices and wholesale prices for a month weighted by market shares.
preserve
* Generate measure of price variance within a month
bysort month: egen price_sd_mon = sd(price)
* Generate measure ofprice variance within a year.
gen year = 1 if month<201301
replace year = 2 if month<201401 & year==.
replace year = 3 if month<201501 & year==.
replace year = 4 if month<201601 & year==.
replace year = 5 if month>=201601 & year==.
bysort year: egen price_sd_year = sd(price)
collapse price price_sd_mon price_sd_year wholesale_spot [iweight=mshare_contract], by(month) 
* Save average price data.
save `"${PATH_OUT_DATA}/avg_price_data.dta"', replace
restore

* Create average advertising measures for a month weighted by market shares.
preserve
*collapse ad_spending_mly ad_spending_norm [iweight=1/mshare_contract], by(month) 
* Should ad spending be weighted by market share? Not clear to me.
collapse ad_spending_mly ad_spending_norm , by(month) 

* Save average price data.
save `"${PATH_OUT_DATA}/avg_adv_data.dta"', replace
restore
