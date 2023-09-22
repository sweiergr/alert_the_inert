/*
	Read in CPI.
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
import excel using `"${PATH_IN_DATA}/CPI_bel.xlsx"', firstrow clear
rename A month
rename priceindex CPI
rename yearlypriceindex CPI_yrl
keep month CPI CPI_yrl
sort month 
* Save as Stata file.
save `"${PATH_OUT_DATA}/CPI_bel"', replace
