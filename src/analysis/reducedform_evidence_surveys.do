/*
	This file runs several reduced-form regs on the survey data to investigate 
	the importance of demographics for awareness and switching behavior.
	This file generates Table 1 in the main text of the paper.
*/

* BASICS.
clear
capture log close
version 13
set more off
set mem 4g
* Header do-File with path definitions, those end up in global macros.
include project_paths
log using `"${PATH_OUT_ANALYSIS}/log/`1'.log"', replace

*********************************
* Load survey data 
use `"${PATH_OUT_DATA}/survey_controlvars.dta"', clear

* Sort data set.
sort year
order year

*drop all years that we do not use in the master dataset
drop if year <2012
**************************

**************************
*rename variables for table
label var vtest "Fully informed"
label var size_hh "Household size"
label var gender "Woman"
label var senior "Senior"
label var higher_ed "Higher education"
label var lower_ed "Primary education"
label var income "Family net income"
label var affordability_strong "Energy costs important"
label var affordability_weak "Energy costs important"
label var incumbent "Incumbent supplier"
label var green_contract "Green contract"
label var year "Year"
**************************

est clear

******************************************************************
* AWARENESS on demographics
* (awareness is measured as being fully informed = "having done vtest")

local awareness_reg "size_hh gender senior higher_ed lower_ed income affordability_strong incumbent"
corr `awareness_reg'

probit vtest `awareness_reg' year
probit vtest `awareness_reg' year green_contract

*NOTE
est clear
eststo: probit vtest size_hh gender senior higher_ed lower_ed income affordability_weak
eststo: probit vtest size_hh gender senior higher_ed lower_ed income affordability_weak year 
eststo: probit vtest size_hh gender senior higher_ed lower_ed income affordability_weak year green_contract

esttab using `"${PATH_OUT_TABLES}/rf_vtest.tex"', ///
cells(b(star fmt(3)) se(par fmt(3))) replace ///
			legend label collabels(none) ///
			stats(N, fmt(0 2) label(Observations) ) ///
			eqlabels(none) /// 
star(* 0.10 ** 0.05 *** 0.01) nocons ///
addnotes("\textit{Data source}: VREG surveys 2012-2016." ///
 		 "\textit{Notes:} The table summarizes results from probit regressions with a dummy for" ///
		 "a consumer having used the PCW on  consumer characteristics. \emph{Household size}" ///
		 "is the number of people living in the household. \emph{Woman} and \emph{Senior} are dummies" ///
		 "for female consumers and respondents older than 65, respectively. \emph{Higher education} and" ///
		 "\emph{Primary education} are dummies for consumers with a higher education degree and primary" ///
		 "education degree only, respectively. \emph{Family net income} is the monthly net household income." ///
		 "\emph{Energy costs important} denotes consumers who state that energy costs are an important" ///
		 "part of their budget. \emph{Green contract} indicates consumers who currently receive energy" ///
		 "only from renewable sources. \emph{Year} captures a linear time trend." ///
		 "Standard errors in parentheses.") ///
mtitles("Fully informed" "Fully informed" "Fully informed")
est clear

******************************************************************
******************************************************************
* SWITCHING on awareness and demo
*	(both past switching behavior and future switching intentions)

*******************
*(1) history = "have you ever changed electricity supplier -> yes/no
gen switched=.
replace switched=1 if sw_history==1
replace switched=0 if sw_history==2
replace switched=0 if sw_history==. & contract==2
*tab year if sw_history==. & !(contract==2)
*these are the ones that do not know
label var switched "Dummy =1 if person has switched supplier in past"
label define switchedl 0 "no" 1 "yes"  
label values switched switchedl
*drop those that do not know
drop if sw_history==3

*****************
*rename for table
label var switched "Past switch"
*****************

local switched_reg "vtest size_hh gender senior higher_ed lower_ed income affordability_strong"
eststo: probit switched `switched_reg' year
probit switched `switched_reg' incumbent year
probit switched `switched_reg' incumbent year green_contract
eststo: probit switched `switched_reg' year
probit switched `switched_reg' year green_contract
probit switched `switched_reg' incumbent year green_contract ad_spending_sector
*******************
*******************
*(2) switching intention = "Do you consider switching in the next 6 months?"

gen switching_int=0
replace switching_int=1 if sw_intention==3 | sw_intention==4 | sw_intention_default==3 | sw_intention_default==4

*****************
*rename for table
label var switching_int "Int. switch"
*****************


local switched_reg "vtest size_hh gender senior higher_ed lower_ed income affordability_strong"
eststo: probit switching_int `switched_reg' year
probit switching_int `switched_reg' year green_contract

esttab, star(* 0.10 ** 0.05 *** 0.01) b(3) se(3) nocons compress ///
label stats(N , fmt(%9.0f) labels(`"Observations"'))

label var affordability_strong "Energy costs imp."

est clear
local switched_reg "vtest  senior  income "
eststo: probit switched `switched_reg'
eststo: probit switched `switched_reg' year 
eststo: probit switched `switched_reg' year size_hh higher_ed lower_ed affordability_strong 

eststo: probit switching_int `switched_reg'
eststo: probit switching_int `switched_reg' year 
eststo: probit switching_int `switched_reg' year size_hh higher_ed lower_ed affordability_strong 

esttab, star(* 0.10 ** 0.05 *** 0.01) b(3) se(3) nocons compress ///
label stats(N , fmt(%9.0f) labels(`"Observations"'))

esttab using `"${PATH_OUT_TABLES}/rf_switching.tex"', ///
star(* 0.10 ** 0.05 *** 0.01) nocons ///
cells(b(star fmt(3)) se(par fmt(3))) replace ///
			legend label collabels(none) ///
			stats(N, fmt(0 2) label(Observations) ) ///
			eqlabels(none) /// 
addnotes("\textit{Data source}: VREG surveys 2012-2016." ///
 		 "\textit{Notes:} The table summarizes results from probit regressions with a dummy for whether" ///
		 "a consumer has switched in the past (Columns 1-3) or intends to switch supplier" /// 
		 "(Columns 4-6) on consumer characteristics. The definition of the regressors is as in" ///
		 "Table B.1. Standard errors in parentheses.") ///
mtitles("Past sw." "Past sw." "Past sw." "Intention" "Intention" "Intention")
label var affordability_strong "Energy costs important" //LD: no need for reduced table size

******************************************************************
* BEING WITH INCUMBENT on demographics

local incumbent_reg "vtest size_hh gender senior higher_ed lower_ed income affordability_strong"
probit incumbent `incumbent_reg' year  
probit incumbent `incumbent_reg' year green_contract
******************************************************************
******************************************************************
* SUPPLIER CHOICE and demographics

local sup_choice_reg "vtest size_hh gender senior higher_ed lower_ed income affordability_strong"
mlogit firm `sup_choice_reg' year  
******************************************************************
******************************************************************
* GREEEN CONTRACT CHOICE and demographics
local green_choice_reg "vtest size_hh gender senior higher_ed lower_ed income affordability_strong incumbent"
probit green_contract `green_choice_reg' year  
******************************************************************
est clear

******************************************************************
* PRICE CHOICE and demographics
* (price measured as the average per supplier and contract_type in each year)
*************************
*rename for table
label var price_total "Average price"
*************************
local price_choice_reg "vtest size_hh gender senior higher_ed lower_ed income affordability_weak"
eststo: reg price_total `price_choice_reg' year 
eststo: reg price_total `price_choice_reg' year incumbent
eststo: reg price_total `price_choice_reg' year incumbent green_contract
esttab, star(* 0.10 ** 0.05 *** 0.01) b(3) se(3) nocons compress label stats(r2 N , fmt(%9.3f %9.0f) labels(`"R2"' `"Observations"'))
est clear

** Suggestions for new "better-choice regressions".
* Fill sw_frequency variable with zeros.
replace sw_frequency = 0 if sw_frequency ==.
label var sw_frequency "No of past switches"
local price_choice_reg "vtest size_hh gender senior higher_ed lower_ed income affordability_strong"
// corr `price_choice_reg'
local price_choice_reg_sml "vtest senior income"
local price_choice_reg_med "vtest senior income size_hh gender"
local price_choice_reg_lrg "vtest senior income size_hh  gender affordability_weak"
// reg min_avg_price `price_choice_reg' year green_contract
eststo: reg price_total `price_choice_reg_sml' 
eststo: reg price_total `price_choice_reg_med' 
eststo: reg price_total `price_choice_reg_lrg' 
eststo: reg price_total `price_choice_reg_lrg' sw_frequency
esttab, star(* 0.10 ** 0.05 *** 0.01) b(3) se(3) nocons compress label stats(r2 N , fmt(%9.3f %9.0f) labels(`"R2"' `"Observations"'))
* Save results to file.
esttab using `"${PATH_OUT_TABLES}/rf_betterchoice.tex"',  ///
cells(b(star fmt(3)) se(par fmt(3))) replace ///
			legend label collabels(none) ///
			stats(N, fmt(0 2) label(Observations) ) ///
			eqlabels(none) /// 
star(* 0.10 ** 0.05 *** 0.01) nocons ///
addnotes("\textit{Data source}: VREG surveys 2012-2016." ///
 		 "\textit{Notes:} The table summarizes results from OLS regressions of a consumer's average" ///
		 "monthly electricity expenditure on consumer characteristics. \emph{No of past switches}" ///
		 "captures how often a consumer has already switched suppliers. The remaning regressors" ///
		 "are defined as in Table B.1. Standard errors in parentheses.") ///
mtitles("Average price" "Average price" "Average price" "Average price")
*****************************************************************

est clear

* Redo some of the regressions for a more concise summary table for the main text.
* One regression illustrating who uses vtest website.
eststo: probit vtest size_hh gender senior higher_ed lower_ed income affordability_weak year green_contract
* Two regressions for illustrating which consumers consider switching and have siwtched frequently in the past.
eststo: probit switched `switched_reg' year size_hh higher_ed lower_ed affordability_strong 
eststo: probit switching_int `switched_reg' year size_hh higher_ed lower_ed affordability_strong 
* One regression illustrating who gets better deals.
eststo: reg price_total `price_choice_reg_lrg' sw_frequency
* Table for working paper.
esttab using `"${PATH_OUT_TABLES}/rf_maintext.tex"',  ///
cells(b(star fmt(3)) se(par fmt(3))) replace ///
			legend label collabels(none) ///
			stats(N, fmt(0 2) label(Observations) ) ///
			eqlabels(none) /// 
star(* 0.10 ** 0.05 *** 0.01) nocons ///
addnotes("\textit{Data source}: VREG surveys 2012-2016." ///
 		 "\textit{Notes:} The table summarizes reduced form regression results that analyze" ///
		 "heterogeneity in PCW usage and contract switching behavior. Standard errors in " ///
		 "parentheses. Detailed explanations and alternative specifications for each column " ///
		 "are presented in Appendix B.") ///
mtitles("Fully informed" "Past sw." "Intent sw." "Average price")

* Table for AEJ:Micro
esttab using `"${PATH_OUT_TABLES}/rf_maintext_aej.tex"',  ///
cells(b(fmt(3)) se(par fmt(3))) replace ///
			legend label collabels(none) ///
			stats(N, fmt(0 2) label(Observations) ) ///
			eqlabels(none) /// 
 nocons /// star(* 0.10 ** 0.05 *** 0.01)
mtitles("Fully informed" "Past switch" "Intention" "Average price")

********************************************************************************
* Save data set.
save `"${PATH_OUT_DATA}/rf_reg_surveys.dta"', replace

