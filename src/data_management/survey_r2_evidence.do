/*
	This file explores some survey statistics for consumers awareness about the market environment. Moreover, it prints some statistics on aggregate advertising by firms. This information is not part of the paper, but only used
    in our reply to referee 2.

*/

* BASICS.
clear
capture log close
version 13
set more off
set mem 4g
* Header do-File with path definitions, those end up in global macros.
include project_paths
log using `"${PATH_OUT_DATA}/log/`1'.log"', replace

*********************************
* Load data for survey in 2012.
*********************************
* Load data for survey in 2012.
use `"${PATH_OUT_DATA}/surveydata2012_clean.dta"', clear
* Append data for all years that are cleaned and conformable.
append using `"${PATH_OUT_DATA}/surveydata2013_clean.dta"', force
append using `"${PATH_OUT_DATA}/surveydata2014_clean.dta"', force
append using `"${PATH_OUT_DATA}/surveydata2015_clean.dta"', force
append using `"${PATH_OUT_DATA}/surveydata2016_clean.dta"', force
* Sort data set.
sort year
order year
* Save data set.
save `"${PATH_OUT_DATA}/surveydata_allyears.dta"', replace
* Some summary demographic characteristics.
gen income_group = (income <=3)
replace income_group = 2 if income>3 & income <=6
replace income_group = 3 if income>6 & income <=9
gen senior = (age > 64)
tab sw_reliable
tab sw_cheaper
tab sw_service
tab info_liberalization
tab no_sw_effort
tab no_sw_how
tab no_sw_possibility 
tab conscious
tab aware_vtest
* These statistics are used in our reply to R2.2
tab no_sw_possibility
tab no_sw_how
tab info_liberalization
which they would switch
tab sw_saving
sum sw_saving_nb
* Look at potential savings compared to incumbent supplier.
use `"${PATH_OUT_DATA}/new_master_data.dta"', clear
bysort month: gen price_inc_aux = price if firm==1 & contract==1
bysort month: egen price_inc = max(price_inc_aux)
* Savings compared to incumbent.
gen savings_inc = 100*(price_inc - price)
sum savings_inc
* Load monthly advertising data.
use `"${PATH_OUT_DATA}/adv_nielsen_long_outside.dta"', clear
drop if firm==7
sum ad_spending_mly
list if ad_spending_mly==0
