/*
	Read in Internet usage rates by demograpic type from OECD data.	
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
import excel using `"${PATH_IN_DATA}/internet_demo.xlsx"', sheet(data) firstrow clear
label var demo "demographic split: age and education"
label var usage "individuals using the internet in last 12 month (5) (source OECD)"
* Save as Stata file.
save `"${PATH_OUT_DATA}/internet_demo"', replace
