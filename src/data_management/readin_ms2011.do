/*
	This file readsin the MS in Dec 2011, which is needed to forward all MS
	to t+1 without loosing first period. This is needed link the data to our model framework.

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
import excel using `"${PATH_IN_DATA}/ms2011.xlsx"', sheet("data") firstrow clear
* Save as Stata file.
save `"${PATH_OUT_DATA}/ms2011"', replace
