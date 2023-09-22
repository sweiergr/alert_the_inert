*******************************************************************************
* This file creates several control variables from the survey data to be used in reduced form regressions. 
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

*********************************
* Load data for survey in 2012.
use `"${PATH_OUT_DATA}/surveydata2011_clean.dta"', clear

* Append data for all years that are cleaned and conformable.
append using `"${PATH_OUT_DATA}/surveydata2004_clean.dta"', force
append using `"${PATH_OUT_DATA}/surveydata2005_clean.dta"', force
append using `"${PATH_OUT_DATA}/surveydata2006_clean.dta"', force
append using `"${PATH_OUT_DATA}/surveydata2007_clean.dta"', force
append using `"${PATH_OUT_DATA}/surveydata2008_clean.dta"', force
append using `"${PATH_OUT_DATA}/surveydata2009_clean.dta"', force
append using `"${PATH_OUT_DATA}/surveydata2010_clean.dta"', force
append using `"${PATH_OUT_DATA}/surveydata2012_clean.dta"', force
append using `"${PATH_OUT_DATA}/surveydata2013_clean.dta"', force
append using `"${PATH_OUT_DATA}/surveydata2014_clean.dta"', force
append using `"${PATH_OUT_DATA}/surveydata2015_clean.dta"', force
append using `"${PATH_OUT_DATA}/surveydata2016_clean.dta"', force

* Sort data set.
sort year
order year

*create firm variable
gen firm =.
replace firm=1 if supplier==12 //ECS
replace firm=2 if supplier==25 | supplier==32 //EDF
replace firm=3 if supplier==17 //Eneco
replace firm=4 if supplier==45 | supplier==27 //ENINuon
replace firm=5 if supplier==20 //Essent
replace firm=6 if supplier==24 //Lampiris
replace firm=7 if firm == .
label define firml 1 "ECS" 2 "EDF" 3 "Eneco" 4 "ENINuon" 5 "Essent" 6 "Lampiris" 7 "Other"
label values firm firml

*construct a contract_type variable (to assign a price)
*	- 50% of sample were not asked whether or not  green -> avg price green/conv
*	- 50% that were asked whether green or not green -> green or conv price
*	- the once on default get default price
gen contract_type_new=4
*add default contract only in 2012
replace contract_type_new=2 if default==1 & year<=2012 & !(year==2004 | year==2006 | year==2007)
*tab default
*tab year if default==.

*add green contract
*tab year green_contract
replace contract_type_new=1 if green_contract==2	
replace contract_type_new=3 if green_contract==1 & !(year==2004 | year==2006 | year==2007)
label define contract_type_newl2 1 "conventional" 2 "default" 3 "green" 4 "green/conv"
label values contract_type_new contract_type_newl2

numlabel , add

replace contract_type=4 if firm==7 & contract_type==2
replace contract_type=3 if firm==6 | firm==3
//NOTE: Lampiris and Eneco can only be on a green contract
//		these respondents got it wrong
replace contract_type=4 if firm==7

*******************************************************************************
* Create control variables
*****************
*income
replace income = . if income==10
gen mis_income=0
replace mis_income=1 if income==.
label var mis_income "Dummy =1 if person did not give income information (or info was not asked)"
*****************

*****************
*Importance of energycost in HH budget as a proxy for income?
tab energycost
gen affordability_weak=0
replace affordability_weak=1 if energycost==3 | energycost==4
label var affordability_weak "Dummy =1 if energycost are rather or very important in HH budget"

gen affordability_strong=0
replace affordability_strong=1 if energycost==4
label var affordability_strong "Dummy =1 if energycost are very important in HH budget"

replace affordability_weak=. if year==2004
replace affordability_strong=. if year==2004
*****************

*****************
* Education: just use binary variable for higher education.
gen higher_ed = 0
replace higher_ed = 1 if education==3 | education==4
label var higher_ed "Dummy =1 if person received uni or non-uni higher education"

gen lower_ed = 0
replace lower_ed = 1 if education==1
label var lower_ed "Dummy =1 if person received only primary education"
*****************

*****************
* Employment status as proxy for age, i.e. look at retirees
gen senior = 1 if employment==5
replace senior = 0 if senior==.
replace senior = . if employment==.
label var senior "Dummy =1 if person is retired ~ proxy for age"
*Workaround because some info does not exist in 2016
replace senior = 0 if age<45 & senior==1
//in 2016
replace senior = 0 if age<63 & year==2016
replace senior = 1 if age>=63 & year==2016
*****************

*****************
* working
gen working = 0
replace working = 1 if employment==1 
label var working "Dummy =1 if person is currently working"
*2016
replace working=1 if active16==1
*****************

*****************
* consumption //NOTE: not available in 2004 and 2005
gen mis_cons=0
replace mis_cons=1 if consumption==6
label var mis_cons "Dummy =1 if person did not give consumption information"
*Clean up consumption variable
replace consumption=. if consumption==6

* Recode consumption levels.
gen cons_level = 1 if consumption==1
replace cons_level = 2 if consumption==2 | consumption==3
replace cons_level = 3 if consumption==4 | consumption==5

//Consumption levels by person -> as a proxy?
gen cons_person=.
replace cons_person=450/size_hh if consumption==1
replace cons_person=1625/size_hh if consumption==2
replace cons_person=3200/size_hh if consumption==3
replace cons_person=8700/size_hh if consumption==4
replace cons_person=20000/size_hh if consumption==5
label var cons_person "Proxy for consumption by person"
*****************

*****************
gen gender_aux = 0 if gender==1
replace gender_aux = 1 if gender==2
drop gender
rename gender_aux gender
label var gender "Dummy =1 if person is a woman"
*****************

*****************
gen ownership_aux = 0
replace ownership_aux = 1 if ownership==1
drop ownership
rename ownership_aux ownership
label var ownership "Dummy =1 if person is owner (vs tenant)"
*****************

*****************
gen hometype_aux = .
replace hometype_aux = 0 if hometype==2
replace hometype_aux = 1 if hometype==1
drop hometype
rename hometype_aux hometype
label var hometype "Dummy =1 if person lives in house (vs apartment)"
*****************

*****************
* children //was not asked in 2004 and 2005
tab year size_hh if kids==. 
gen kids_aux = . 
replace kids_aux=0 if kids==. & year>2005
replace kids_aux=kids if !(kids==.)
tab kids_aux
drop kids
gen kids=kids_aux
drop kids_aux
label var kids "Number of kids living in HH"

*one is strange
replace kids=9 if kids==99
*****************

*****************
gen elecheating=. //was not asked <2012
replace elecheating=0 if year >=2012
replace elecheating=1 if heating1==1 | heating1==2 | heating2==2
label var elecheating "Dummy =1 if HH heats with elec"
tab elecheating
*****************

*****************
*green contract //Was not asked in 2004, was asked only to 50% of sample
gen green_contract_aux=0 if green_contract==2 | green_contract==3
replace green_contract_aux=1 if green_contract==1 
drop green_contract
gen green_contract=green_contract_aux
label var green_contract "Dummy =1 if HH has a green contract (only 50% of sample was asked)"
drop green_contract_aux
*****************

**************************
*incumbent dummy
gen incumbent = 1 if firm==1
replace incumbent = 0 if !(firm==1)
label var incumbent "Dummy =1 if HH has a contract with ECS"
**************************

*****************
*Has done the Vtest -> to be used in master_dataset to construct fully_aware_share
*Be careful different subset of people were asked the question whether they 
*are aware or have done the vtest
* This does not account for the basis --> Adjusted below
gen vtest=.
replace vtest=1 if done_vtest==1
replace vtest=0 if done_vtest==2
label var vtest "Dummy =1 if person has done the vtest"
tab vtest
gen new_vtest=0 if year>=2009
replace new_vtest=. if aware_vreg==. & (year==2012 | year==2013)
//these guys were NOT part of the sample to be asked about Vreg and Vtest
//all other with 0 were
replace new_vtest=1 if vtest==1
*drop if new_vtest==.
*tab year vtest
*tab year new_vtest
drop vtest
rename new_vtest vtest 
label var vtest "Dummy =1 if person has done the vtest (conditional on being asked)"
*****************

*****************
*Knows Vtest
* This does not account for the basis --> Adjusted below
gen knows_vtest=.
replace knows_vtest=1 if aware_vtest==1
replace knows_vtest=0 if aware_vtest==2
label var knows_vtest "Dummy =1 if person knows the vtest"
tab knows_vtest
tab year aware_vtest

gen new_knows_vtest=0 if year>=2009
replace new_knows_vtest=. if aware_vreg==. & (year==2012 | year==2013)
//these guys were NOT part of the sample to be asked about Vreg and Vtest
//all other with 0 were
replace new_knows_vtest=1 if knows_vtest==1
*drop if new_knows_vtest==.
*tab year knows_vtest
*tab year new_knows_vtest

drop knows_vtest
rename new_knows_vtest knows_vtest 
label var knows_vtest "Dummy =1 if person knows vtest (conditional on being asked)"
*****************

*****************
*Advertisement 

******
* (1) yearly spending by firm in 2012-2016
* Contstruct from adv_nielsen_long_outside -> no data after 201606!!
preserve

use `"${PATH_OUT_DATA}/adv_nielsen_long_outside"', clear

gen year=int(month/100)
drop month 

by year firm, sort: egen ad_spending_yrl=sum(ad_spending_mly)
label var ad_spending_yrl "Total real media spending by brand per year (in 1000 EUR) in Belgium"
*delete duplicates
qui by year firm, sort: gen dup=cond(_N==1,0,_n)
drop if dup>1
drop dup ad_spending_mly
order year firm
sort year firm
tempfile spending
save `spending'
restore
**
merge m:1 year firm using `spending'
drop _merge

******
* (2) yearly spending by sector
preserve
keep year firm ad_spending_yrl
drop if ad_spending_yrl==.
*drop duplicates
qui by year firm, sort: gen dup=cond(_N==1,0,_n)
drop if dup>1
drop dup
*create total by sector by year
by year, sort: egen ad_spending_sector=total(ad_spending_yrl)
label var ad_spending_sector "Total energy sector real media spending per year (in 1000 EUR) in Belgium"

*drop duplicates
qui by year, sort: gen dup=cond(_N==1,0,_n)
drop if dup>1
drop dup ad_spending_yrl

tempfile sector_spending
save `sector_spending'

restore
**
merge m:1 year using `sector_spending'
drop _merge
*****************
*****************
* prices
sort year
order year

*******
* (1) Annual bill for 3500 kWh HH
merge m:m year firm contract_type_new using `"${PATH_OUT_DATA}/new_prices_survey.dta"'
tab year _merge
*problem with Essent conventional
/* check manually: 
	- 2015: price_total 322.35
			min_avg_price 280.3225
	- 2016: price_total 305.133
			min_avg_price 276.5749
*/
replace price_total=322.35 if firm==5 & contract_type==1 & year==2015
replace min_avg_price=280.3225 if firm==5 & contract_type==1 & year==2015
replace price_total=305.133 if firm==5 & contract_type==1 & year==2016
replace min_avg_price=276.5749 if firm==5 & contract_type==1 & year==2016
drop if _merge==2 
drop _merge
*********************************
#delimit ;
order year firm incumbent green_contract vtest knows_vtest higher_ed lower_ed senior 
	gender size_hh kids ad_spending_yrl ad_spending_sector price_total min_avg_price 
	max_avg_price income mis_income affordability_weak affordability_strong 
	working cons_level mis_cons cons_person ownership hometype elecheating
	sw_relative sw_relative_excl
;
#delimit cr
*****************
*save dataset 
save `"${PATH_OUT_DATA}/survey_controlvars.dta"', replace
*******************************************************************************

