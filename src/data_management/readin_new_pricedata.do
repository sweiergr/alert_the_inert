/*
	Read in price data table and and save as Stata-dta-file.
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
import excel using `"${PATH_IN_DATA}/newprices.xlsx"', firstrow clear

*******************************************************************************
* Electricity in Flanders
preserve

drop if Supplier==""
drop if Type=="GAZ"
drop if Place=="Wallonia" | Place=="Brussels"
drop Type Place 

* Save as Stata file.
save `"${PATH_OUT_DATA}/new_pricedata_FL"', replace

restore

*******************************************************************************
* Elec in Wallonia
preserve

drop if Supplier==""
drop if Type=="GAZ"
drop if Place=="Flanders" | Place=="Brussels"
drop Type Place 

* Save as Stata file.
save `"${PATH_OUT_DATA}/new_pricedata_elec_wl.dta"', replace

restore

*******************************************************************************
* Elec in Brussels
preserve

drop if Supplier==""
drop if Type=="GAZ"
drop if Place=="Flanders" | Place=="Wallonia"
drop Type Place 

* Save as Stata file.
save `"${PATH_OUT_DATA}/new_pricedata_elec_bxl.dta"', replace

restore

*******************************************************************************
* Gas in Flanders
preserve

drop if Supplier==""
drop if Type=="ELEC"
drop if Place=="Wallonia" | Place=="Brussels"
drop Type Place 

* Save as Stata file.
save `"${PATH_OUT_DATA}/new_pricedata_gas_fl.dta"', replace

restore

*******************************************************************************
* Gas in Wallonia
preserve

drop if Supplier==""
drop if Type=="ELEC"
drop if Place=="Flanders" | Place=="Brussels"
drop Type Place 

* Save as Stata file.
save `"${PATH_OUT_DATA}/new_pricedata_gas_wl.dta"', replace

restore

*******************************************************************************
* Gas in Brussels
preserve

drop if Supplier==""
drop if Type=="ELEC"
drop if Place=="Flanders" | Place=="Wallonia"
drop Type Place 

* Save as Stata file.
save `"${PATH_OUT_DATA}/new_pricedata_gas_bxl.dta"', replace

restore
*******************************************************************************
