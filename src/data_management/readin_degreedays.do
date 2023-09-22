/*
	Read in temperature time series and save as Stata-dta-file.
	Note that this is not used in the final model specification.

*/

clear
capture log close
version 13
set more off
set mem 4g
* Header do-File with path definitions, those end up in global macros.
include project_paths
log using `"${PATH_OUT_DATA}/log/`1'.log"', replace

* Import data set from Excel files in src/original_data folder.
import excel using `"${PATH_IN_DATA}/BEL_HDDCDD.xls"', sheet("Data") firstrow clear
* Adjust variable names.
	foreach v of varlist B-P {
	   local x : variable label `v'
	   rename `v' year`x'
	   rename year`x' degreedays`x'
	}
* Rename month
gen month=.
replace month=1 if Mois=="Janvier"
replace month=2 if Mois=="Fevrier"
replace month=3 if Mois=="Mars"
replace month=4 if Mois=="Avril"
replace month=5 if Mois=="Mai"
replace month=6 if Mois=="Juin"
replace month=7 if Mois=="Juillet"
replace month=8 if Mois=="Aout"
replace month=9 if Mois=="Septembre"
replace month=10 if Mois=="Octobre"
replace month=11 if Mois=="Novembre"
replace month=12 if Mois=="Decembre"
drop if Mois=="Total" 
order month
drop Mois
reshape long degreedays, i(month) j(year)
gen yrmonth = year*100+month
drop if year<2004
drop month
rename yrmonth month
order month degreedays
sort year month
by year: egen degreedays_year = total(degreedays)
order month year 
* Save as Stata file.
save `"${PATH_OUT_DATA}/degreedays"', replace
