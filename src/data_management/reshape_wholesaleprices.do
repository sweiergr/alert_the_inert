	/*
		Reshape wholesale price data (spot and forward) to be compatible with 
		other data. In addition, transform nominal into real (Base: jan 2012)
	*/

	clear
	capture log close
	version 13
	set more off
	set mem 4g

	* Header do-File with path definitions, those end up in global macros.
	include project_paths
	log using `"${PATH_OUT_DATA}/log/`1'.log"', replace

	* Load wholesale contract data.
	use `"${PATH_OUT_DATA}/wholesale_contracts_raw"', clear

	* drop data before 2012 for now
	drop if year <2012

	* Drop duplicate information on date:
	drop day month year yearmonthdaycsv 
	drop day_type week_number weekday_char weekday_num

	*rename yearmonth identifyier as in other datasets
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

* drop data before 2012 for now
drop if year <2012

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

sort month
save `"${PATH_OUT_DATA}/wholesale_prices"', replace

