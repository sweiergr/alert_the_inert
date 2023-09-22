/*
	Read in price data table and and save as Stata-dta-file.
	This file is not used anymore in the most recent specification.

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
import excel using `"${PATH_IN_DATA}/prices_mod.xlsx"', firstrow clear

foreach v of varlist I-AJ {
   local x : variable label `v'
   rename `v' month`x'
}
drop month
* Save as Stata file.
save `"${PATH_OUT_DATA}/pricedata_FL"', replace

* Import price data for gas in Flanders.
* Create file name.
	local fname_in = `"${PATH_IN_DATA}/prices_"'+"gas"+"_"+ "fl"+".xlsx"
	local fname_out = `"${PATH_OUT_DATA}/pricedata_"'+"gas"+"_"+ "fl"+".dta"
	import excel "`fname_in'", sheet("Contracts") firstrow clear

	* Adjust variable names.
	foreach v of varlist I-AJ {
	   local x : variable label `v'
	   rename `v' month`x'
	   rename month`x' gas_fl_month`x'
	}
	drop gas_fl_month
drop if Supplier==""
capture drop AK
	* Save as dta-file in bld/out/data folder.
	save "`fname_out'", replace

* Import price data for electricity and gas in other regions.
local type_list gas elec
local region_list bxl wl
foreach type of local type_list {
	foreach region of local region_list{
	* Create file name.
	local fname_in = `"${PATH_IN_DATA}/prices_"'+"`type'"+"_"+ "`region'"+".xlsx"
	local fname_out = `"${PATH_OUT_DATA}/pricedata_"'+"`type'"+"_"+ "`region'"+".dta"
	import excel "`fname_in'", sheet("Contracts") firstrow clear

	* Adjust variable names.
	foreach v of varlist I-AJ {
	   local x : variable label `v'
	   rename `v' month`x'
	   rename month`x' `type'_`region'_month`x'
	}
	drop `type'_`region'_month
drop if Supplier==""
capture drop AK
	* Save as dta-file in bld/out/data folder.
	save "`fname_out'", replace
}
}

