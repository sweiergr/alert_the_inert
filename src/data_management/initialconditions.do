*******************************************************************************
* This file creates a dataset for the initial conditions based on (micromoments.do)
* These data are now based on January 2012 with sanity checks based on data from December 2011.
*******************************************************************************
* BASICS.
clear
capture log close
version 13
set more off
set mem 4g
* Header do-File with path definitions, those end up in global macros.
include project_paths
log using `"${PATH_OUT_DATA}/log/`1'.log"', replace

******************************************************************************
* Load survey data 
use `"${PATH_OUT_DATA}/survey_controlvars.dta"', clear
*construct labels for contract type
keep if year==2011 | year==2010
drop contract_type

gen contract_type=.
replace contract_type=1 if green_contract==0    
replace contract_type=2 if default==1 & year==2012
replace contract_type=3 if green_contract==1
label define contract_typel 1 "conventional" 2 "default" 3 "green" 4 "green/conv"
label values contract_type contract_typel
*some respondents got it wrong -> Lampiris and Essent can only be with a green contract
replace contract_type=3 if (firm==6 | firm==3) & !(contract_type==.)
*do not distinguish contract_types in outside option
replace contract_type=1 if firm==7 & !(contract_type==.)
******************************

******************************
*we only need the following variables in the dataset
keep year firm contract_type income senior vtest
*clean
order year firm contract_type income senior vtest
sort year firm contract_type income senior vtest
*for 2010 and 2011 data (not relevant anymore)
drop vtest
drop year
numlabel , add
******************************

******************************
*Drop the observations with missing values
drop if income==. 
drop if contract_type==.	//in micromoments we have kept them
******************************

******************************
*reshape to have the consumer types as columns
*new income variable
gen inc=1 if  income <=3
replace inc=2 if income >3 & income <=5
replace inc=3 if income >5 & income <=6
replace inc=4 if income >6
drop income
label def inclabel 1 "<	1500" 2 "1500-2500" 3 "2500-3750" 4 ">3750"
label values inc inclabel
numlabel , add
* reshape senior
gen senior2=0
replace senior2=1 if senior==0
rename senior senior1
collapse (sum) senior1 senior2, by(firm contract_type inc)
reshape long senior, i(firm contract_type inc) j(new)
rename senior observations
rename new age
label def agelabel 1 "senior" 2 "nonsenior"
label values age agelabel
numlabel , add
*create total by demographic combination
by inc age, sort: egen total = sum(observations)
*MS
gen share = observations/total
drop observations total
sort firm contract_type inc age
*matrix
gen combi = inc*10+age
label def combilabel 11 "senior inc_1" 21 "senior inc_2" 31 "senior inc_3" 41 "senior inc_4" 12 "nonsenior inc_1" 22 "nonsenior inc_2" 32 "nonsenior inc_3" 42 "nonsenior inc_4"
label values combi combilabel
numlabel , add
drop inc age
reshape wide share, i(firm contract_type) j(combi)
rename share11 senior_inc_1
rename share21 senior_inc_2
rename share31 senior_inc_3
rename share41 senior_inc_4
rename share12 nonsenior_inc_1
rename share22 nonsenior_inc_2
rename share32 nonsenior_inc_3
rename share42 nonsenior_inc_4
*******************************************************************************
******************************
* a) import aggregate data (only used for sanity checks and robustness checks)
preserve
use `"${PATH_OUT_DATA}/ms2011.dta"', clear
drop if mshare==.
drop mshare
replace mshare_contract = mshare_contract[_n]+mshare_contract[_n+1] if (firm==1|firm==2) & contract==1
rename mshare_contract agg_ms
drop if contract_type==2
tempfile agg_ms
save `agg_ms'
restore
******************************
******************************
* b) merge survey and aggregate
merge 1:1 firm contract using `agg_ms'
drop _merge
drop month
drop agg_ms
replace senior_inc_1 =0 if senior_inc_1==.
replace senior_inc_2 =0 if senior_inc_2==.
replace senior_inc_3 =0 if senior_inc_3==.
replace senior_inc_4 =0 if senior_inc_4==.
replace nonsenior_inc_1 =0 if nonsenior_inc_1==.
replace nonsenior_inc_2 =0 if nonsenior_inc_2==.
replace nonsenior_inc_3 =0 if nonsenior_inc_3==.
replace nonsenior_inc_4 =0 if nonsenior_inc_4==.
sort firm contract
* Order variables to make them compatible with MATLAB code.
order firm contract_type nonsenior* senior*
*******************************************************************************
save `"${PATH_OUT_DATA}/initialconditions.dta"', replace
export delimited using `"${PATH_OUT_DATA}/initialconditions.csv"', replace nolabel
******************************************************************************

