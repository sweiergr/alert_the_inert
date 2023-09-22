/*
	Constructs a dataset of demographics 
	- fully informed -> share of survey respondents that did the vtest
	- age -> share of population aged > 65 in Flanders
	- education -> share of population (25-64) that has a atteint tertiary 
				   education in Flanders (not used in main model)
	- income -> income distribution across quartiles in Belgium
	

	Data on 
	- income and education did not exist for 2016 yet
		-> 2015 data is used in 2016 !! 
*/


clear
capture log close
version 13
set more off
set mem 4g
* Header do-File with path definitions, those end up in global macros.
include project_paths
log using `"${PATH_OUT_DATA}/log/`1'.log"', replace

********************************************************************************
* Information on number of Flemish households. This allows us to compute 
* various statistics to compute share of informed households.
import delimited `"${PATH_IN_DATA}/households_bel.csv"',  clear
* Rename and rescale variables as neccessary.
replace population_fl = population_fl * 1000000
gen build_fl = houses_fl + apartments_fl
gen hh_fl_1 = households_bel * share_pop_fl
gen hh_fl_2 = population_fl / hhsize_fl
* Save Stata file with information on number of Flemish households.
save `"${PATH_OUT_DATA}/hh_fl.dta"', replace

******************************************************************************
*AGE
********
use `"${PATH_OUT_DATA}/demo_age_long"', clear
gen month =year*100+1
*append 11 more month for age
forvalues i = 2(1)12{
	append using `"${PATH_OUT_DATA}/demo_age_long"', force
	replace month=year*100+`i' if month==.
}


******************************************************************************
*merge EDUCATION
********
merge m:1 year using `"${PATH_OUT_DATA}/demo_edu_long"'
replace share_higher_ed=37.2 if share_higher_ed==.
drop _merge

******************************************************************************
*merge INCOME
********
merge m:1 year using `"${PATH_OUT_DATA}/demo_inc_wide"'

drop cutoff*
rename share_inc1 share_inc_low
label var share_inc_low "Share of BELGIAN population in lower income quartile"
rename share_inc2 share_inc_middle
label var share_inc_middle "Share of BELGIAN population in both middle income quartiles"
rename share_inc3 share_inc_high
label var share_inc_high "Share of BELGIAN population in higher income quartile"
drop _merge
*we use the 2015 shares in 2016 (unlikely to change significantly form 2015 to 2016)
replace share_inc_low=12.3 if year==2016
replace share_inc_middle=46.1 if year==2016
replace share_inc_high=41.6 if year==2016
replace share_inc_low = share_inc_low / 100
replace share_inc_middle = share_inc_middle / 100
replace share_inc_high = share_inc_high / 100
tempfile master
save `master'


******************************************************************************
* merge FULLY INFORMED
********
*first need to construct a fully_informed share by year
preserve
use `"${PATH_OUT_DATA}/survey_controlvars.dta"', clear
keep year vtest
*do not use the ones with missing value for vtest in the calculations
drop if vtest==.
by year, sort: egen numerator=sum(vtest)
by year, sort: gen denominator=_N
*drop duplicates
quietly by year, sort:  gen dup = cond(_N==1,0,_n)
drop if dup>1
drop dup
*construct share
gen share_fully_inf=numerator/denominator
label var share_fully_inf "Share of survey repondents that did the vtest"
drop vtest num* denom*
tempfile fully_inf
save `fully_inf'
restore
********
merge m:m year using `fully_inf'
drop if _merge==2
drop _merge


********************************************************************************
* Calculation of vtest usage numers start here.
* Bring in variation into monthly vtest shares using the clicks on simulator
preserve
use `"${PATH_OUT_DATA}/macrodataraw_FL.dta"', clear
keep date year month simulator
* Convert date numbes to actual date type.
tostring month, gen(month_str) format(%20.0f)
gen emonth = date(month_str,"M")
tostring date, gen(date_str) format(%20.0f)
gen edate = date(date_str,"YM")
tostring year, gen(year_str) format(%20.0f)
gen eyear = date(year_str,"Y")

* Plot raw simulator usage numbers.
twoway (line simulator edate, sort title("Raw usage numbers of vtest"))
graph export `"${PATH_OUT_FIGURES}/simulator_raw.pdf"', replace

* According to VREG the simulator number in April 2012 must be an outlier.
* Replace it with average of surrounding numbers.
replace simulator = 0.5 * simulator[_n-1] + 0.5 * simulator[_n+1] if date==201204
* June 2012 is also likely an outlier. Do same imputation here.
replace simulator = 0.5 * simulator[_n-1] + 0.5 * simulator[_n+1] if date==201206
twoway (line simulator edate, sort title("Usage numbers of vtest - outliers cleaned"))
graph export `"${PATH_OUT_FIGURES}/simulator_clean.pdf"', replace
* Compute and plot hypothetical usage rate based on number of households in Flanders in 2012.
gen vtest_rate = simulator / 2646000
twoway (line vtest_rate edate, sort title("Share of vtest users (relative to no of HH)"))
graph export `"${PATH_OUT_FIGURES}/vtest_rate_hh.pdf"', replace
bysort year: sum vtest_rate
sum vtest_rate

********************************************************************************
*construct shares by month
* Generate survey group dummies.
gen yr_svy = 1 if date<=201208
replace yr_svy = 2 if date>201208 & date<=201308
replace yr_svy = 3 if date>201308 & date<=201408
replace yr_svy = 4 if date>201408 & date<=201508
replace yr_svy = 5 if date>201506
sort yr_svy month
bysort yr_svy: egen num_2_svy=sum(simulator)
gen share_sim_svy=simulator/num_2_svy
label var share_sim_svy "Monthly clicks on Vtest as a share of survey-yearly clicks"
* check whether sums up to 1
by yr_svy, sort: egen sum_mw=total(share_sim_svy)
sum sum_mw
drop sum_mw
twoway (line share_sim_svy edate, sort title("Monthly share of vtest users (relative to total clicks in survey period)"))
graph export `"${PATH_OUT_FIGURES}/vtest_monthsharesvperiod.pdf"', replace

* LD's old verison.
by year, sort: egen numerator2=sum(simulator)
gen share_sim=simulator/numerator2
label var share_sim "Monthly clicks on Vtest as a share of yearly clicks"
twoway (line share_sim edate, sort title("Monthly share of vtest users (relative to total clicks in calendar year)"))
graph export `"${PATH_OUT_FIGURES}/vtest_monthshareyear.pdf"', replace
* check whether sums up to 1
by year, sort: egen sum_mw=total(share_sim)
sum sum_mw
drop sum_mw

*clean
capture drop *year* *month* numerator
*drop simulator
rename date month
tempfile sim_shares
save `sim_shares'
restore
*******
merge m:1 month using `sim_shares'
drop if _merge==1 //was not merged for 2011 - although we have the data
drop if _merge==2 // nor for the very late periods: >=201610
drop _merge
sort month

gen share_inf_mtly=share_sim*share_fully_inf
bysort year: egen cumshare_inf_yr = sum(share_inf_mtly)
gen share_inf_mtly_svy=share_sim_svy*share_fully_inf
bysort yr_svy: egen cumshare_inf_svy = sum(share_inf_mtly_svy)
label var share_inf_mtly "Monthly share of fully informed (yrly vtest*mthly proportion of clics)"


* Compare evolution of vtest users implied by different survey/year classifications.
*check whether sum up to "share_fully_inf"
by year, sort: egen sum=total(share_inf_mtly) 
gen diff=round(share_fully_inf-sum)
sum diff
*ok

* We construct a few alternatives to the vtest numbers constructed above.
merge m:1 year using `"${PATH_OUT_DATA}/hh_fl.dta"'
drop if year==.
drop _merge

* Candidates for aggregate vtest usage numbers are the following:
gen share_sim_1 = simulator / 2646000
label var share_sim_1 "Fixed HH"
gen share_sim_2 = simulator / hh_fl_1
label var share_sim_2 "Actual HH"
label var share_inf_mtly "Survey - Year"
label var share_inf_mtly_svy "Survey - Period"
twoway (line share_inf_mtly edate, sort) (line share_inf_mtly_svy edate, sort) (line share_sim_2 edate, sort), xline(19237 19390) title(Comparison of vtest user numbers for different date classifications)
graph export `"${PATH_OUT_FIGURES}/vtest_agg_comp.pdf"', replace
*clean
capture drop share_sim sum diff
drop share_fully_inf
* We use this version based on raw clicks.
rename share_sim_2 share_fully_inf
drop if year<2012 | year>2016
bysort year: sum share_fully_inf
sum share_fully_inf
******************************************************************************
*clean dataset for saving
********
drop year
order month share_fully_inf
sort month 
* Save Stata file.
save `"${PATH_OUT_DATA}/demo_help.dta"', replace
* Drop months that we currently do not use in the structural estimation.
drop if month > 201606
* Export required data to csv for reuse  in MATLAB.
export delimited using `"${PATH_OUT_DATA}/demo_data.csv"', replace nolabel
**************************************************************************
