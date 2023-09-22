/*
	Read in advertisement data sets (1) from UBA, (2) from Nielsen. 
	Not all of these data sets are used in the final estimation specification.

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
* Import UBA dataset from Excel files in src/original_data folder.
foreach t in bel north south uba{
	* Create file name.
	local fname_in = `"${PATH_IN_DATA}/advertisement_"'+"`t'"+".xlsx"
	local fname_out = `"${PATH_OUT_DATA}/advertisement_"'+"`t'"
	import excel "`fname_in'", sheet("data") firstrow clear
	reshape long yr, i(group) j(year)
	rename yr ad_spending
	drop if group=="regional"
	* Save as dta-file in bld/out/data folder.
	save "`fname_out'", replace
}

******************************************************************************
*Import Nieslen dataset (monthly 2014-2016)
import excel `"${PATH_IN_DATA}/adv_nielsen_mly.xls"', sheet("data") firstrow clear
* Save as dta-file in bld/out/data folder.
saveold `"${PATH_OUT_DATA}/adv_nielsen"', replace
******************************************************************************
*Import Nieslen dataset (yearly 2011-201s)
import excel `"${PATH_IN_DATA}/adv_nielsen_yrly.xls"', sheet("Adex-result") firstrow clear
drop H
* Save as dta-file in bld/out/data folder.
saveold `"${PATH_OUT_DATA}/adv_nielsen_yrly"', replace

