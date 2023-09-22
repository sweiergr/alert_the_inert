/*
	Constructs number of elec and gas consumers by supplier (2011-2016)
	- get share (mly) by supplier in elec and gas market
	- get total size (yrl) of elec and gas market
	To be used to normalize advertisement data. In the most recent version, this is 
	not used anymore, except for a robustness check.

*/

clear
capture log close
version 13
set more off
set mem 4g
* Header do-File with path definitions, those end up in global macros.
include project_paths
log using `"${PATH_OUT_DATA}/log/`1'.log"', replace
use `"${PATH_OUT_DATA}/marketshares_long"', clear
* Merge gas market shares (monthly) and size of market (from BREG)
merge m:1 month firm using `"${PATH_OUT_DATA}/data_macro_gas_FL"'
//the 2011 elec shares are from the using dataset
replace mshare=mshare_2011 if mshare==.
drop _merge mshare_2011
sort month firm
*calculate number of customers 
gen nb_elec=mshare*size_elec
gen nb_gas=gas_mshare*size_gas
gen nb_both=nb_elec+nb_gas
*calculate number of non-customers (those with other suppliers)
gen nb_nonelec=(1-mshare)*size_elec
gen nb_nongas=(1-gas_mshare)*size_gas
gen nb_nonboth=nb_nonelec+nb_nongas
*label variables
label var nb_elec "Number of customer in elec market"
label var nb_gas "Number of customer in gas market"
label var nb_both "Number of customer in both markets (elec and gas)"
label var nb_nonelec "Number of non-customer in elec market"
label var nb_nongas "Number of non-customer in gas market"
label var nb_nonboth "Number of non-customer in both markets (elec and gas)"
*plot number of customers over time
gen yr=int(month/100)
gen m=mod(month,100)
gen date=ym(y,m)
label var date "Jan 2012 - June 2016"
drop yr m date
*Clean
drop mshare gas_mshare size_elec size_gas
* Save data.
order month firm 
sort month firm
save `"${PATH_OUT_DATA}/nb_customer"', replace
********************************
