/*
	Reshape macro/market share data to be compatible with our model framework.
	It also creates the separate switching rate dataset.
*/

clear
capture log close
version 13
set more off
set mem 4g
* Header do-File with path definitions, those end up in global macros.
include project_paths
log using `"${PATH_OUT_DATA}/log/`1'.log"', replace

********************************
* Load macro data in wide format.
* Create file name.
local fname_in = `"${PATH_OUT_DATA}/macrodataraw_FL"'
use "`fname_in'", clear
* Aggregate data for small firms into outside good.
* We put all firms with market sahre constantly below 1% into the outside good.
replace Other = Other + Belpower + EBEM + Elegant + OCTA + WATZ
drop Belpower EBEM Elegant OCTA WATZ
* Drop information on passive consumers, network and simulator visits for now.
drop passive* network simulator
* Drop region as long as we focus on Flanders.
drop region 
* Drop duplicate information on date.
drop year month
rename date month
label var month "YearMonth"

********************************
* Rescale switching variable into decimals and save in separate file.
replace switching = switching / 100
preserve
keep month switching
label var switching "Monthly switching rates (supplier switch) as reported by VREG"
save `"${PATH_OUT_DATA}/switching_rates_FL"', replace
restore
drop switching
********************************

* Reshape market share data into long format.
* Initialize count variable
local i = 0
foreach var of varlist ECS-Other{
	local i = `i' + 1
	* Rescale shares into deicmals.
	replace `var' = `var' / 100
	* Rename variables to be compatible with reshape command.
	rename `var'  mshare`i'
	}
* Reshape into market shares into long format.
reshape long mshare, i(month) j(firm)

* Might be worthwhile to check whether it's possible to keep firm label, but not crucial.
label define firml 1 "ECS" 2 "EDF" 3 "Eneco" 4 "ENINuon" 5 "Essent" 6 "Lampiris" 7 "Other"
label values firm firml
*add numberlabel
numlabel , add
* Save market share data set.
order month firm // contract_type
sort month firm //contract_type
save `"${PATH_OUT_DATA}/marketshares_long"', replace
********************************
