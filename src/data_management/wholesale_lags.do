/*
	Create a dataset with lagged wholesale price IVs (i.e. interaction with dependence on wholesale price), wholesale prices are transformed into 2012-EUR.
*/

clear
capture log close
version 13
set more off
set mem 4g

* Header do-File with path definitions, those end up in global macros.
include project_paths
log using `"${PATH_OUT_DATA}/log/`1'.log"', replace
***********THE BEGINNING OF THIS DOFILE COPIES "RESHAPE_WHOLESALEPRICES"
***********ONLY EXCEPTION starts at 2012-6 months
* Load wholesale contract data.
use `"${PATH_OUT_DATA}/wholesale_contracts_raw"', clear
//THESE TWO LINES ARE DIFFERENT TO "RESHAPE_WHOLESALEPRICES"
* drop data before 201107
drop if year <2011
drop if year <2012 & month <7 

* Drop duplicate information on date:
drop day month year yearmonthdaycsv 
drop day_type week_number weekday_char weekday_num
rename yearmonth month

* Average data by month (no weights)
rename apxfbepeq1 contractprice
destring contractprice, replace dpcomma
drop if contractprice==.
sort month
by month: egen wholesale_contract=mean(contractprice)
label var wholesale_contract "Monthly average of daily price of elctricity contracts (delivery next quarter)"

*drop daily variables
drop yearmonthday contractprice

*drop duplicates
sort month wholesale_contract
quietly by month wholesale_contract:  gen dup = cond(_N==1,0,_n)
*tab dup
drop if dup>1
drop dup

tempfile wholesale_contracts
save `wholesale_contracts'


* Load wholesale spot data.
use `"${PATH_OUT_DATA}/wholesale_spot_raw"', clear

****************
//THESE TWO LINES ARE DIFFERENT TO "RESHAPE_WHOLESALEPRICES"
* drop data before 201107
drop if year <2011
drop if year <2012 & month <7 
****************

* Drop duplicate information on date:
drop day month year yearmonthdaycsv day_type week_number weekday_char
drop weekday_num quarter_num hour quarter_in_hour pop_elia	pop_belpex

*rename yearmonth identifyier as in other datasets
rename yearmonth month

* Average data by month (no weights)
rename bpxbbpxprc spotprice
destring spotprice, replace dpcomma
drop if spotprice==.
sort month
by month: egen wholesale_spot=mean(spotprice)
label var wholesale_spot "Monthly average of quarterhourly price on electricity spot market Belpex"

*drop daily variables
drop yearmonthday quarter spotprice

*drop duplicates
sort month wholesale_spot
quietly by month wholesale_spot:  gen dup = cond(_N==1,0,_n)
*tab dup
drop if dup>1
drop dup

* Merge both wholesaleprices
merge 1:1 month using `wholesale_contracts'
drop _merge

*************************
*deflate prices
merge m:1 month using `"${PATH_OUT_DATA}/CPI_bel"'
drop if _merge==2
drop _merge

rename wholesale_spot wholesale_spot_nominal
rename wholesale_contract wholesale_contract_nominal

gen wholesale_spot=wholesale_spot_nominal/CPI
label var wholesale_spot "Monthly average of real quarter-hourly price on elec spot market Belpex"

gen wholesale_contract=wholesale_contract_nominal/CPI
label var wholesale_contract "Monthly average of real daily price of elec contracts (delivery next quarter)"

drop CPI CPI_yrl wholesale_spot_nominal wholesale_contract_nominal
*************************

*************************
//THE FOLLOWING SECTION IS DIFFERENT FROM "RESHAPE_WHOLESALEPRICES"

*create wholesale price lags - 6 lags 
forval i = 1/6 {
	gen ws_spot_lag_`i' = wholesale_spot[_n-`i']
	label var ws_spot_lag_`i' "`i' month lag of real mthly avg price on elec spot market Belpex"
}
forval i = 1/6 {
	gen ws_future_lag_`i' = wholesale_contract[_n-`i']
	label var ws_future_lag_`i' "`i' month lag of real mthly avg  price on elec contract (delivery next quarter)"
}
******************************************************************************
* create 2 different datasets (one for spot lags, one for future lags) with 594 observations as in new_masterdata
*drop data before 2012
drop if month <201201
*drop months in 2016 for which we do not have data in the master dataset
drop if month>201606
merge 1:m month using `"${PATH_OUT_DATA}/new_master_data"', keepusing(month firm contract_type)
drop _merge
order month firm contract_type
sort month firm contract_type

******************************
*create the interaction with wholesale dependence
*borrowed from merge_new_masterdata.do

* Merge supplier dependence on wholesale market
merge m:1 firm using `"${PATH_OUT_DATA}/dep_wholesale"'
drop _merge
order month firm contract_type
sort month firm contract_type
* Rescale dependence
replace share_purchase=share_purchase+1
* Generate instruments based on interaction of wholesale contract prices and production share
gen ws_future_int = wholesale_contract * share_purchase
forval i = 1/6 {
	gen ws_future_lag_int_`i' = ws_future_lag_`i' * share_purchase
	}
* only keep the interactions
keep month firm contract_type ws_future_lag_int*

sort month firm contract_type
save `"${PATH_OUT_DATA}/ws_future_int_lags"', replace
* Export data set to csv file for import in MATLAB.
export delimited using `"${PATH_OUT_DATA}/ws_future_int_lags.csv"', novarnames nolabel replace

