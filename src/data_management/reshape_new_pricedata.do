/*
	Reshape price data to be compatible with market share data and	transform nominal into real prices (Base: Jan 2012)
*/

clear
capture log close
version 13
set more off
set mem 4g

* Header do-File with path definitions, those end up in global macros.
include project_paths
log using `"${PATH_OUT_DATA}/log/`1'.log"', replace
**************************
* Import data set from Excel files in src/original_data folder.
use `"${PATH_OUT_DATA}/new_pricedata_FL"', clear
**************************
*Prepare firms
rename Supplier firm
drop if firm==""

replace firm="ENINuon" if firm=="ENI" | firm=="Nuon" | firm=="NUON"
replace firm="ECS" if firm=="ELECTRABEL"
replace firm="EDF" if firm=="LUMINUS"
replace firm="Other" if !(firm=="ENECO" | firm=="ESSENT" | firm=="LAMPIRIS" ///
	|firm=="ENINuon" | firm=="ECS" | firm=="EDF")
**************************

**************************
*For now, clean contract components that we do not need
drop Subscription*
drop Rate*
drop Extra*
**************************

**************************
*Clean contract attributes
*Propose to only use prices for single  rates (vs double-hourly rates)
sort firm VarFi
drop if MonoBi=="Bi" | MonoBi=="BI"
drop MonoBi
*For now, only keep duration=1 (although different durations exist)
drop Duration
*Drop quantity that are not compatibly with 3500kWh/year
destring Quantitymin, replace

replace Quantitymin=0 if Quantitymin==.
drop if Quantitymin>3500
replace Quantitymax=5001 if Quantitymax==.
drop if Quantitymax<3500
drop Quantitymin Quantitymax
**************************

**************************
*TODO: for now drop the VARIO contracts from EDF, as they have no contract MS
drop if contract_type=="delete_noms"
*TODO: two other have no contract MS either drop them for now
drop if inms=="no"
**************************
**************************
*drop variables that are not needed
drop Name VarFix ID inms
rename new_id contract_id
**************************
**************************
*for now, we do not need the contract id
drop contract_id
**************************
**************************
*construct a simple average over contracts by contract_type and by firm
collapse Total*, by(contract_type firm)
replace contract_type="conventional" if firm=="Other"
* EDF default
replace Total201201=Total201205 if firm=="EDF" & contract_type=="default"
replace Total201202=Total201205 if firm=="EDF" & contract_type=="default"
replace Total201203=Total201205 if firm=="EDF" & contract_type=="default"
replace Total201204=Total201205 if firm=="EDF" & contract_type=="default"
* ECS green
replace Total201604=Total201603 if firm=="ECS" & contract_type=="green"
replace Total201605=Total201603 if firm=="ECS" & contract_type=="green"
replace Total201606=Total201603 if firm=="ECS" & contract_type=="green"
replace Total201607=Total201603 if firm=="ECS" & contract_type=="green"
replace Total201608=Total201603 if firm=="ECS" & contract_type=="green"
replace Total201609=Total201603 if firm=="ECS" & contract_type=="green"
* ENI green
replace Total201606=Total201605 if firm=="ENINuon" & contract_type=="green"
replace Total201607=Total201605 if firm=="ENINuon" & contract_type=="green"
replace Total201608=Total201605 if firm=="ENINuon" & contract_type=="green"
replace Total201609=Total201605 if firm=="ENINuon" & contract_type=="green"
*************************

*************************
* Reshape price data into long format.
reshape long Total, i(firm contract_type) j(month)
rename Total price_total
label var price_total "Simple average of contract price by supplier and contract_type"
label var month "YearMonth"
*************************

*************************
*Encode firm variable
sort firm month
encode firm, gen(firm_num)
drop firm
rename firm_num firm
label define firml 1 "ECS" 2 "EDF" 3 "Eneco" 4 "ENINuon" 5 "Essent" 6 "Lampiris" 7 "Other"
label values firm firml

*Encode contract_type variable
sort contract_type month
encode contract_type, gen(contract_type_num)
drop contract_type
rename contract_type_num contract_type
label define contract_typel 1 "conventional" 2 "default" 3 "green" 
label values contract_type contract_typel

*add numberlabel
numlabel , add
*************************

*************************
*deflate prices
merge m:1 month using `"${PATH_OUT_DATA}/CPI_bel"'
drop if _merge==2
drop _merge

rename price_total price_nominal

gen price_total=price_nominal/CPI
label var price_total "Simple average of real contract price by supplier and contract_type"

drop CPI CPI_yrl price_nominal
*************************

order month firm contract_type
sort month firm contract_type

*clean
drop if price_total==. | price_total==0

* check whether we have each combi (firm / contract_type) in each period.
tab month contract_type

* Save as Stata file.
save `"${PATH_OUT_DATA}/new_pricedata_FL_long"', replace



*******************************************************************************
* Additional: prepare yearly avg prices to be merged with survey 

* construct a simple yearly average over contract types by supplier
*	- 50% of sample were not asked whether or not  green -> avg price green/conv
*	- 50% that were asked whether green or not green -> green or conv price
*	- the once on default get default price

* create the mix type: avg price green/conv
sort firm month
reshape wide price_total, i(firm month) j(contract_type)

gen price_total4=(price_total1 + price_total3)/2
replace price_total4=price_total1 if price_total3==.
replace price_total4=price_total3 if price_total1==.
reshape long price_total, i(firm month) j(contract_type_new)

* create prices for green vs conventional vs default
gen year=int(month/100)

collapse price_total, by(contract_type_new firm year)
order year firm
sort year firm contract_type_new
drop if price_total==.
label var price_total "Yearly electrictiy bill of 3500kWh HH (avg by firm and contract_type)"
label define contract_type_newl2 1 "conventional" 2 "default" 3 "green" 4 "green/conv"
label values contract_type_new contract_type_newl2
numlabel, add

*construct default price by year for each of the two suppliers
*attach avg default price to suppliers that are not default
gen default=price_total if contract_type==2
sort year
by year: egen default_price=mean(default)
label var default_price "Yearly electrictiy bill of 3500 HH who pays default price"
drop default
replace default_price=price_total if contract_type==2

*construct min average price by year
sort firm year
egen min_avg_price = min(price_total), by(year)
label var min_avg_price "Yearly electrictiy bill of 3500 kWh HH who pays smallest of avg contract prices"

*construct max average price by year
sort firm year
egen max_avg_price = max(price_total), by(year)
label var max_avg_price "Yearly electrictiy bill of 3500 kWh HH who pays highest of avg contract prices"

*************************
*deflate prices
gen month=year*100+1
merge m:1 month using `"${PATH_OUT_DATA}/CPI_bel"'
drop if _merge==2
drop _merge CPI month

rename price_total price_total_nominal
gen price_total=price_total_nominal/CPI_yrl
label var price_total "Yearly real electrictiy bill of 3500kWh HH (avg by firm and contract_type)"

rename default_price default_price_nominal
gen default_price=default_price_nominal/CPI_yrl
label var default_price "Yearly real electrictiy bill of 3500 HH who pays default price"

rename min_avg_price min_avg_price_nominal
gen min_avg_price=min_avg_price_nominal/CPI_yrl
label var min_avg_price "Yearly real electrictiy bill of 3500 kWh HH who pays smallest of avg contract prices"

rename max_avg_price max_avg_price_nominal
gen max_avg_price=max_avg_price_nominal/CPI_yrl
label var max_avg_price "Yearly real electrictiy bill of 3500 kWh HH who pays smallest of avg contract prices"
drop CPI_yrl price_total_nominal default_price_nominal min_avg_price_nominal max_avg_price_nominal
*************************
sort contract_type firm year
save `"${PATH_OUT_DATA}/new_prices_survey.dta"', replace
*******************************************************************************

