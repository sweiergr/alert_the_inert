/*
	Read in macro/market share data sets and and save as Stata-dta-file.
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
local fname_in = `"${PATH_IN_DATA}/data_macro_FL.xlsx"'
local fname_out = `"${PATH_OUT_DATA}/macrodataraw_FL"'
import excel "`fname_in'", sheet("data_FL") cellrange(A2:U59) firstrow clear
* Save as dta-file in bld/out/data folder.
save "`fname_out'", replace



