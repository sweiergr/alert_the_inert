/*
	Read in wholesale price table and and save as Stata-dta-file.
*/

clear
capture log close
version 13
set more off
set mem 4g
* Header do-File with path definitions, those end up in global macros.
include project_paths
log using `"${PATH_OUT_DATA}/log/`1'.log"', replace
* CONTRACTS - daily, delivery in next trimester
* Import data set from csv file in src/original_data folder.
import delimited using `"${PATH_IN_DATA}/wholesale_contracts.csv"', delimiters(";") varnames(1)
* Save as Stata file.
save `"${PATH_OUT_DATA}/wholesale_contracts_raw"', replace
clear
*SPOT MARKET - quarterly hourly
* Import data set from csv file in src/original_data folder.
import delimited using `"${PATH_IN_DATA}/wholesale_spot.csv"', delimiters(";") varnames(1)
* Save as Stata file.
save `"${PATH_OUT_DATA}/wholesale_spot_raw"', replace
