/*
	Read in data on whoelsale natural gas prices.

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
* Create file name.
local fname_in = `"${PATH_IN_DATA}/wholesaleGasPrices.xlsx"'
local fname_out = `"${PATH_OUT_DATA}/wholesaleGasPrices"'
import excel "`fname_in'", sheet("Wholesale gas prices ") cellrange(A12:E171) clear

* Label variables.
rename A date
rename B wgp_NL
rename C wgp_US
rename D wgp_Asia
rename E wgp_GER
capture drop date_aux

* Reformat month variable.
gen month_num = month(date)
gen year_num = year(date)
tostring month_num, gen(month_str)
tostring year_num, gen(year_str) 
replace month_str = "0" + month_str if strlen(month_str)==1
gen date_str = year_str + month_str
destring date_str, gen(month)

* Drop unnecessary date variables.
drop date month_num month_str year_num year_str date_str
order month
sort month
* Save as dta-file in bld/out/data folder.
save "`fname_out'", replace




