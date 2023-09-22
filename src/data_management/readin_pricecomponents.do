/*
	Read in price data split up into components and and save as  Stata-dta-file.

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
import excel using `"${PATH_IN_DATA}/prices_component.xlsx"', firstrow clear

drop if Supplier=="" | Supplier=="Supplier"
drop if Type=="GAZ"
drop if Place=="Wallonia" | Place=="Brussels"
drop Type Place 
drop EP EQ

* Save as Stata file.
save `"${PATH_OUT_DATA}/pricedata_component_FL"', replace
