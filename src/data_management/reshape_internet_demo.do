/*
	Reshape internet usage data.
*/

clear
capture log close
version 13
set more off
set mem 4g
* Header do-File with path definitions, those end up in global macros.
include project_paths
log using `"${PATH_OUT_DATA}/log/`1'.log"', replace

**************************
* Import data set from Excel files in src/original_data folder.
use `"${PATH_OUT_DATA}/internet_demo"', clear
drop demo
*create average for non-seniors
sort year educ age
by year educ, sort: egen usage_mean=mean(usage) if !(age=="senior")
replace usage=usage_mean if !(age=="senior")
drop usage_mean
drop if age=="young"
gen senior=0
replace senior=1 if age=="senior"
drop age
preserve
drop if !(educ=="total")
drop educ
order year senior
* Save in Stata format.
save `"${PATH_OUT_DATA}/internet_usage_age.dta"', replace
* Export in csv-format to read into MATLAB.
export delimited using `"${PATH_OUT_DATA}/internet_usage_age.csv"', nolabel replace
restore
*******************************************************
drop if educ=="total"
order year senior educ
* Save in Stata format.
save `"${PATH_OUT_DATA}/internet_usage_educ.dta"', replace
* Export in csv-format to read into MATLAB.
export delimited using `"${PATH_OUT_DATA}/internet_usage_educ.csv"', nolabel replace
*******************************************************************************
