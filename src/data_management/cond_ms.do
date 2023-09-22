/*

This file creates conditional market share distribution based on awareness

THIS FILE IS NOT USED ANYMORE IN THE MSOT RECENT VERSION OF THE PAPER.


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

******************************************************************************
* Load data 
use `"${PATH_OUT_DATA}/survey_controlvars.dta"', clear

*keep the relevant years
drop if year<=2007 //no question on green
drop if year <=2010 //issue with the default tariff
drop if year==2011 //we do not use it for now

*construct labels for contract type -> only 50% were asked
drop contract_type

gen contract_type=.
replace contract_type=1 if green_contract==0	
replace contract_type=2 if default==1 & year==2012 //only exists in year 2012	
replace contract_type=3 if green_contract==1
label define contract_typel 1 "conventional" 2 "default" 3 "green" 4 "green/conv"
label values contract_type contract_typel
//those with missing values are the observations, where we do not know 
//contract type for sure

numlabel , add

*some people got it wrong -> Lampiris and Essent can only be with a green contract
replace contract_type=3 if (firm==6 | firm==3) & !(contract_type==.)
*do not distinguish contract_types in outside option
replace contract_type=1 if firm==7 & !(contract_type==.)

*keep only those variables that we need
keep year firm contract_type vtest senior income
******************************************************************************

******************************************************************************
* Awareness

preserve

*do not build MS if we do not know the contract type for sure
drop if contract_type==.
//this drops half of th observations: 51% are left

tab year vtest

/*One issue witht the done_vtest question ist that not everyone was asked
 in every year
 -> however this is controlled for in the vtest variable
	See the reducedform_controlvars_survey.do file for more info*/ 
	
*when caclulating MS do not incl those persons that were not asked whether (s)he has done the vtest
drop if vtest==.
//drops another 567
	
*generate conditional market share
keep year firm contract_type vtest 

*TODO: for ENI no one (of the 8 obs) on a green contract in 2012 
*	   has done the vtest
*(in 2013, 7 out of 19 have done the vtest)
*browse if year==2012 & firm==4 & contract_type==3
*browse if year==2013 & firm==4 & contract_type==3
*LD: shall we give them a share of zero?

**** LD: Workaround for ENI having a 0 share for green-aware in 2012
gen id=_n
replace vtest=1 if id==797
drop id
**** NEEDS to be solved differently

sort firm contract_type year vtest
by firm contract_type year vtest: gen numerator=_N

sort year vtest
by year vtest: gen denominator=_N

gen ms_aware=numerator/denominator

*drop trash
drop numerator denominator

order year firm contract_type ms_aware vtest 

*drop duplicates
quietly by year firm contract_type vtest, sort:  gen dup = cond(_N==1,0,_n)
drop if dup>1
drop dup

*check whether sums up to 1
by year vtest, sort: egen sum=total(ms_aware)
sum sum, detail
*ok
drop sum

/*
*check whether all firms have all contract types by awareness type
drop if year==2011 //we do not use it for now
tab year contract_type if firm==1
	//ECS: no default at all as of 2014
tab year contract_type if firm==2
	//EDF: no default at all as of 2014
tab year contract_type if firm==3
	//ENECO: fine
tab year contract_type if firm==4
	//ENI: unbalanced in 2012 and 2011 -> TODO
tab year contract_type if firm==5
	//ESSENT: fine
tab year contract_type if firm==6
	//Lampiris: fine 
tab year contract_type if firm==3
	//Outside option: fine
*/

* SW: Sort data according to year-firm-contract
sort year firm contract_type vtest

**** SW: Workaround for dropping default contract of EDF.
replace ms_aware = ms_aware + ms_aware[_n+2] if firm==2 & contract_type==1 & year == 2012
*replace ms_aware = ms_aware + ms_aware[_n+2] if firm==2 & contract_type==1 & year < 2014
drop if firm==2 & contract_type==2


save `"${PATH_OUT_DATA}/cond_ms_awareness.dta"', replace
* Drop months that we currently do not use in the structural estimation.
*drop if year > 2013 | year < 2012
* Export required data to csv for reuse  in MATLAB.
export delimited using `"${PATH_OUT_DATA}/cond_ms_aware.csv"', replace nolabel
restore
******************************************************************************

******************************************************************************
* Senior

preserve

*do not build MS if we do not know the contract type for sure
drop if contract_type==.

*generate conditional market share
keep year firm contract_type senior 

sort firm contract_type year senior
by firm contract_type year senior: gen numerator=_N

sort year senior
by year senior: gen denominator=_N

gen ms_senior=numerator/denominator

*drop trash
drop numerator denominator

order year firm contract_type ms_senior senior 

*drop duplicates
quietly by year firm contract_type senior, sort:  gen dup = cond(_N==1,0,_n)
drop if dup>1
drop dup

*check whether sums up to 1
by year senior, sort: egen sum=total(ms_senior)
*ok
drop sum

* SW: Sort data according to year-firm-contract
sort year firm contract_type senior

**** SW: Workaround for dropping default contract of EDF.
*replace ms_senior = ms_senior + ms_senior[_n+2] if firm==2 & contract_type==1 & year < 2014
replace ms_senior = ms_senior + ms_senior[_n+2] if firm==2 & contract_type==1 & year == 2012
drop if firm==2 & contract_type==2


save `"${PATH_OUT_DATA}/cond_ms_senior.dta"', replace
* Drop months that we currently do not use in the structural estimation.
*drop if year > 2013 | year < 2012
* Export required data to csv for reuse  in MATLAB.
export delimited using `"${PATH_OUT_DATA}/cond_ms_senior.csv"', replace nolabel
restore
******************************************************************************

******************************************************************************
* income
preserve

*do not build MS if we do not know the contract type for sure
drop if contract_type==.

*many people do not answer the question
drop if income==. 

*construct new income variable
*hist income	
xtile quart = income, nq(4)
*tab quart
*tab income
gen income_new=.
replace income_new=1 if quart==1
replace income_new=2 if quart==2 | quart==3
replace income_new=3 if quart==4
tab income_new

/*
*TODO: strange - many people in income==1: do it by hand
gen income_new2=.
replace income_new2=1 if income<=3
replace income_new2=2 if income>=4 & income<=6
replace income_new2=3 if income>=7
tab income_new2
*/

drop income
rename income_new income_cat
label var income_cat "Income category of HH (1=in lowest, 3=in highest quartile)"
label define income_catl 1 "low income" 2 "middle income" 3 "high income"
label values income_cat income_catl
*rename income_new2 income_cat2

*generate conditional market share
keep year firm contract_type income_cat //income_cat2

*TODO: for ENECO, ESSENT, LAMPIRIS and OUTSIDE option, the ms are unbalanced
*LD: shall we give them a share of zero?

**** LD: Workaround for having a 0 share for green-aware in 2012
gen id=_n
	//ESSENT: unbalanced in 2012 -> no high income
	*browse if firm==5
	replace income_cat=3 if  id==724
	//Lampiris: unbalanced in 2012 -> no high income
	*browse if firm==6
	replace income_cat=3 if  id==734
	//ENECO: unbalanced in 2012 and 2013 -> no low in 2012 no high in 2013
	*browse if firm==3
	replace income_cat=1 if  id==694
	replace income_cat=3 if  id==905
	//Outside option: unbalanced in 2013 -> no high income
	*browse if firm==7
	replace income_cat=3 if  id==987
drop id
**** NEEDS to be solved differently

sort firm contract_type year income_cat
by firm contract_type year income_cat: gen numerator=_N

sort year income_cat
by year income_cat: gen denominator=_N

gen ms_income=numerator/denominator

*drop trash
drop numerator denominator

order year firm contract_type ms_income income_cat 

*drop duplicates
quietly by year firm contract_type income_cat, sort:  gen dup = cond(_N==1,0,_n)
drop if dup>1
drop dup

*check whether sums up to 1
by year income_cat, sort: egen sum=total(ms_income)
sum sum, detail
*ok
drop sum


/*check whether all firms have all contract types by awareness type
drop if year==2011 //we do not use it for now
tab year contract_type if firm==1
	//ECS: no default at all as of 2014
tab year contract_type if firm==2
	//EDF: no default at all as of 2014
tab year contract_type if firm==3
	//ENECO: unbalanced in 2012 and 2013 -> TODO
tab year contract_type if firm==4
	//ENI: fine
tab year contract_type if firm==5
	//ESSENT: unbalanced in 2012 		 -> TODO
tab year contract_type if firm==6
	//Lampiris: unbalanced in 2012 		 -> TODO 
tab year contract_type if firm==3
	//Outside option: unbalanced in 2012 -> TODO
*/


* SW: Sort data according to year-firm-contract
sort year firm contract_type income_cat

**** SW: Workaround for dropping default contract of EDF.
*replace ms_income = ms_income + ms_income[_n+3] if firm==2 & contract_type==1 & year < 2014
replace ms_income = ms_income + ms_income[_n+3] if firm==2 & contract_type==1 & year == 2012
drop if firm==2 & contract_type==2

save `"${PATH_OUT_DATA}/cond_ms_income.dta"', replace
* Drop months that we currently do not use in the structural estimation.
*drop if year > 2013 | year < 2012
* Export required data to csv for reuse  in MATLAB.
export delimited using `"${PATH_OUT_DATA}/cond_ms_income.csv"', replace nolabel
restore
******************************************************************************

