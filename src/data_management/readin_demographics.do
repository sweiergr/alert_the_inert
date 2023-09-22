/*
	Read in demographics save as Stata-dta-file.

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
* Import different demographic datasets from Excel files in src/original_data folder.
foreach t in age edu inc{
	* Create file name.
	local fname_in = `"${PATH_IN_DATA}/demo_"'+"`t'"+".xlsx"
	local fname_out = `"${PATH_OUT_DATA}/demo_"'+"`t'"
	import excel "`fname_in'", sheet("data") firstrow clear
	* Save as dta-file in bld/out/data folder.
	save "`fname_out'", replace
}

