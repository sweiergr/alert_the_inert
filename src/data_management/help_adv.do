/*
	THIS FILE IS NOT NEEDED ANYMORE!
	THIS DATASET IS NOT NEEDED ANYMORE -> we now have monthly data from Nielsen
	
	Constructs a temporary adv dataset with monthly data
	- 2014-2016. data is reported monthly by Nielsen
	- 2011-2013: data is constructed from yearly data reproted by Nielsen
				 (assumes a constant spending per month, i.e. /12)

	Constructs a second temporary adv dataset with monthly data
	- that aggregates the outside option differently (only those supppliers
	  that we know: belpower, elegant, octa+ etc.
	- still, monthly data in years 2011-13 are CONSTANT
*/


clear
capture log close
version 13
set more off
set mem 4g
* Header do-File with path definitions, those end up in global macros.
include project_paths
log using `"${PATH_OUT_DATA}/log/`1'.log"', replace

******************************************************************************
use `"${PATH_OUT_DATA}/adv_nielsen_long"', clear

*append 12 month per year 2011, 2012, 2013
forvalues i = 1(1)12{
	append using `"${PATH_OUT_DATA}/adv_nielsen_long_yrly"', force
	drop if year>=2014 & !(year==.)
	replace month=year*100+`i' if month==.
}

drop year ad_spending_yrl
sort month firm

*drop year 2011 for now
drop if month<=201112

*plot spending over time 
* ->looks ok (BUT need to replace const. spending in years 2011-2013)
gen yr=int(month/100)
gen m=mod(month,100)
gen date=ym(y,m)
label var date "Jan 2011 - June 2016"
*xtline ad_spending_mly, t(date) i(firm) overlay tlabel(, format(%tm))
drop yr m date

save `"${PATH_OUT_DATA}/adv_help"', replace
******************************************************************************


******************************************************************************
use `"${PATH_OUT_DATA}/adv_nielsen_long_outside"', clear

*prepare append

forvalues i = 1(1)12{
	append using `"${PATH_OUT_DATA}/adv_nielsen_long_yrly_outside"', force
	drop if year>=2014 & !(year==.)
	replace month=year*100+`i' if month==.
}

drop year ad_spending_yrl
sort month firm

*drop year 2011 for now
drop if month<=201112

*plot spending over time 
* ->looks ok (BUT need to replace const. spending in years 2011-2013)
gen yr=int(month/100)
gen m=mod(month,100)
gen date=ym(y,m)
label var date "Jan 2011 - June 2016"
*xtline ad_spending_mly, t(date) i(firm) overlay tlabel(, format(%tm))
drop yr m date

save `"${PATH_OUT_DATA}/adv_help_outside"', replace
* SW: add export to csv for importing data into MATLAB.
drop if month>201312
*export delimited using `"${PATH_OUT_DATA}/advertisingdata.csv"', replace nolabel
