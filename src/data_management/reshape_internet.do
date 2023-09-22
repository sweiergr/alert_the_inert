/*
	Reshape internet data
	Adding the agg advertisement by month
	CASE A constructs agg advertisement based on raw data (no MA)
	CASE B constructs agg advertisement based on MA 
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
use `"${PATH_OUT_DATA}/internet"', clear

**************
* add monthly aggregated adverisement by month
preserve

**********
/* CASE a) we construct aggregate advertisement on raw data*/
use `"${PATH_OUT_DATA}/adv_nielsen"', clear

*cleaning the raw data follows straight "reshape_advertisingdata.do" 
rename Soussecteur subsector
rename Groupedannonceurs group
rename Annonceur advertiser
rename Marque brand
rename SubMedia mediatype
rename Annee year
rename Anneemois month
rename InvestissementtotalSommes ad_spending_brand
label var ad_spending_brand "Monthly media spending by brand in Belgium"
drop group advertiser subsector year
*destring month
replace month = subinstr(month, "-", "",.) 
destring month, replace
order month brand
sort month
*xx
gen firm=.
replace firm=1 if brand=="ELECTRABEL" | brand=="ENGIE"
replace firm=2 if brand=="EDF" | brand=="EDF LUMINUS" | brand=="LUMINUS" | brand=="LUMINUS.BE"
replace firm=3 if brand=="ENECO"
replace firm=4 if brand=="ENI" | brand=="NUON"
replace firm=5 if brand=="ESSENT"
replace firm=6 if brand=="LAMPIRIS"
*these 4 are spelled out in the "macro market shares"
replace firm=7 if brand=="BELPOWER" | brand=="EBEM" | brand=="ELEGANT" | brand=="OCTA +" 
*these is additionally spelled out in the "contract market shares" or price data
replace firm=7 if brand=="POWEO" | brand=="COMFORT ENERGY" | brand=="ELEXYS"
*missing in advertisement but spelled out in the price data:
*aspiravi, energy people, mega, klinkenberg
drop if firm==.
drop brand mediatype firm

*construct total mediaspending by month
collapse (sum) ad_spending_brand, by(month)

rename ad_spending_brand aggregate_adv_raw
label var aggregate_adv_raw "Aggregate ad spending per month (based on raw data)"

tempfile aggregate_adv_raw
save `aggregate_adv_raw'

restore

**********
* CASE b) we construct aggregate advertisement based on the constructed moving
*          averages of firm advertisement
preserve

use `"${PATH_OUT_DATA}/adv_nielsen_long"', clear
collapse (sum) ad_spending_mly, by(month)

rename ad_spending_mly aggregate_adv
label var aggregate_adv "Aggregate ad spending per month (based on MA of firm adv)"

**********

tempfile aggregate_adv
save `aggregate_adv'

restore

**************
* merge to internet file
merge 1:1 month using `aggregate_adv_raw'
drop if _merge==1 | _merge==2
drop _merge

merge 1:1 month using `aggregate_adv'
drop if _merge==1
drop _merge

* Safety check that months are ordered correctly.
sort month
**************************
* Add monthly percentage of people using vtest.
merge 1:m month using `"${PATH_OUT_DATA}/new_master_data.dta"', keepusing(vtest)
drop if _merge==1 | _merge==2
drop _merge

* Collapse observations to month-level.
collapse fixed_broadband-vtest, by(month)

* Save in Stata format.
save `"${PATH_OUT_DATA}/internet_adv.dta"', replace
* Export in csv-format to read into MATLAB.
export delimited using `"${PATH_OUT_DATA}/internet_adv.csv"', nolabel replace
*******************************************************************************
