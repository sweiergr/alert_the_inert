/*
	This file reshapes the demographics data.
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


******************************************************************
* Age
use `"${PATH_OUT_DATA}/demo_age.dta"', clear

*create share by age group
sort region year
by region year: egen total=sum(number)
gen age_share=number/total
label var age_share "Share by age group, region and year"

*drop 
drop number total

*relabel variables
replace region="1" if region=="Flanders"
replace region="2" if region=="Wallonia"
replace region="3" if region=="Brussels"
destring region, replace
label define regionl 1 "Flanders" 2 "Wallonia" 3 "Brussels" 
label values region regionl

replace age_group="1" if age_group=="<18"
replace age_group="2" if age_group=="18_64"
replace age_group="3" if age_group==">65"
destring age_group, replace
label define age_groupl 1 "< 18" 2 "18-64" 3 ">=65" 
label values age_group age_groupl
numlabel, add


order year region
sort year region age_group
drop if age_group==1 | age_group==2
drop age_group
drop if region>1
drop region
rename age_share share_senior
label var share_senior "Share of Flemish population aged >65"
*save file
save `"${PATH_OUT_DATA}/demo_age_long.dta"', replace
******************************************************************

******************************************************************
* Education
use `"${PATH_OUT_DATA}/demo_edu.dta"', clear
drop if region=="Belgium"
*relabel variables
replace region="1" if region=="Flanders"
replace region="2" if region=="Wallonia"
replace region="3" if region=="Brussels"
destring region, replace
label define regionl 1 "Flanders" 2 "Wallonia" 3 "Brussels" 
label values region regionl
numlabel, add
label var share_higher_ed "Share of population aged 25-64 with higher education"
order year region
sort year region
drop if region>1
drop region
label var share_higher_ed "Share of Flemish population (25-64) with higher education"
save `"${PATH_OUT_DATA}/demo_edu_long.dta"', replace
******************************************************************

******************************************************************
* Income
use `"${PATH_OUT_DATA}/demo_inc.dta"', clear
drop region
rename cutoffEUR cutoff_rev
label var cutoff_rev "Revenue at the cutoff to next quartile"
rename share share_inc
label var share_inc "Share of BELGIAN population by income_cat"
destring year, replace
order year 
*create three income categories: 1=q1, 2=q2+a3, 3=q4
reshape wide cutoff_rev share_inc, i(year) j(quartile) string

rename share_incq1 share_inc1
gen share_inc2=share_incq2+share_incq3
rename share_incq4 share_inc3
drop share_incq*

rename cutoff_revq1 cutoff_rev1
gen cutoff_rev2=cutoff_revq3
rename cutoff_revq4 cutoff_rev3
drop cutoff_revq*

reshape long cutoff_rev share_inc, i(year) j(income_cat) string

destring income_cat, replace

label define income_catl 1 "low income" 2 "middle income" 3 "high income"
label values income_cat income_catl

numlabel, add

*save file - long version
save `"${PATH_OUT_DATA}/demo_inc_long.dta"', replace
******************************

******************************
*save file - wide version
reshape wide cutoff_rev share_inc, i(year) j(income_cat)
order year cutoff*

save `"${PATH_OUT_DATA}/demo_inc_wide.dta"', replace
******************************************************************
