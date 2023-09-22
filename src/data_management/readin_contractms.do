/*
	Read in contract market share data set and and save as Stata-dta-file.

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
local fname_in = `"${PATH_IN_DATA}/contract_ms_long.xlsx"'
local fname_out = `"${PATH_OUT_DATA}/contract_acspts"'
import excel "`fname_in'", sheet("long_totals_corrected") firstrow clear
*rename var
rename contract_acspoints contract_acspts
label var contract_acspts "Number of access points severd by contract in quarter"
label var new_id "name of supplier and contract"
order quarter firm contract_type contract_acspts
* Save as dta-file in bld/out/data folder.
save "`fname_out'", replace

