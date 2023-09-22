/*
	Merge the different data files into one master file for read-in into MATLAB
	for the structural estimation.
*/

clear
capture log close
version 13
set more off
set mem 4g
* Header do-File with path definitions, those end up in global macros.
include project_paths
log using `"${PATH_OUT_DATA}/log/`1'.log"', replace
* Merge REAL price data with market shares
use `"${PATH_OUT_DATA}/new_pricedata_FL_long"', clear
* Merge aggregate supplier market shares (from VREG, monthly)
merge m:1 month firm using `"${PATH_OUT_DATA}/marketshares_long"'
drop _merge
rename price_total price
* Merge contract market shares (based on VREG quarterly data)
merge m:1 month firm contract_type using `"${PATH_OUT_DATA}/contract_ms"'
drop _merge
* We assume Essent_price_con = Essent_price_green if missing
sort month firm contract_type
* Shift all market shares in data by one period to align with our model framework.
replace price = price[_n+1] if firm==5 & contract_type==1 & month>=201404
replace mshare = mshare[_n+1] if firm==5 & contract_type==1 & month>=201404
sort firm contract_type month
gen mshare_new = mshare[_n+1] if month<201609 & !(contract_type==2)
replace mshare_new = mshare[_n+1] if month<201212 & contract_type==2
drop mshare
rename mshare_new mshare
gen mshare_contract_new = mshare_contract[_n+1] if month<201609 & !(contract_type==2)
replace mshare_contract_new = mshare_contract[_n+1] if month<201212 & contract_type==2
drop mshare_contract
rename mshare_contract_new mshare_contract
drop if month==201609
drop if month==201212 & contract_type==2

*add REAL wholesale price instruments
merge m:1 month using `"${PATH_OUT_DATA}/wholesale_prices"'
drop _merge
//NOTE: wholesale_spot ends in 201607

* Merge REAL gas prices in FL and WL (we do not use bxl).
local region_list fl wl
	foreach region of local region_list{
	* Create file name.
	local fname_in = `"${PATH_OUT_DATA}/new_pricedata_long_gas_"'+"`region'"+".dta"
	sort month firm
	merge m:1 month firm using `fname_in'
	drop _merge
}

* Merge REAL elec prices in WL (we do not use bxl).
local fname_in = `"${PATH_OUT_DATA}/new_pricedata_long_"'+"elec"+"_"+ "wl"+".dta"
sort month firm contract_type
merge m:1 month firm contract_type using `fname_in'
drop if _merge==2 //these are the defaults after 201302
drop _merge
sort month firm contract_type
* We assume Essent_price_con = Essent_price_green if missing
sort month firm contract_type
replace elec_wl_price = elec_wl_price[_n+1] if firm==5 & contract_type==1 & month>=201404
*add switching rates
merge m:1 month using `"${PATH_OUT_DATA}/switching_rates_FL"'
drop _merge

*Hausman instruments
/* NOTE: 
	- for now we drop the bxl prices. 
	- the WL ELEC instrument can be used by contract_type.
	- the GAS instruments can be used by firm (though FL and WL often the same)
*/
gen gp_inst = (gas_fl_price + gas_wl_price)/2
gen ep_inst = elec_wl_price

*Generate dummy for green contract.
gen green = 1 if contract_type==3
replace green = 0  if green==.
*Generate additional regressor for incumbent.
gen incumbent = 1 if firm==1
replace incumbent = 0 if incumbent==.
* add vtest shares
merge m:1 month using `"${PATH_OUT_DATA}/demo_help.dta"', keepusing(share_fully_inf)
drop _merge
rename share_fully_inf vtest
* add REAL advertisement spending
merge m:1 month firm using `"${PATH_OUT_DATA}/adv_nielsen_long_outside.dta"'
drop _merge
* Drop months for which we do not have full data
drop if month>201606
sort month firm contract_type
order month firm contract_type mshare mshare_contract incumbent price green wholesale_spot wholesale_contract
********************************************************************************
* merge default to conventional
* (1) weighted average price
sort month firm contract_type
gen weight_conv=mshare_contract/(mshare_contract+mshare_contract[_n+1]) if contract_type==1 & month<201212 & (firm==1|firm==2)
gen weight_def=mshare_contract[_n+1]/(mshare_contract+mshare_contract[_n+1]) if contract_type==1 & month<201212 & (firm==1|firm==2)
replace price = price*weight_conv + price[_n+1]*weight_def if firm==1 & contract_type==1 & month<201212
replace price = price*weight_conv + price[_n+1]*weight_def if firm==2 & contract_type==1 & month<201212
drop weight_conv weight_def
* (2) sum of mshare
sort month firm contract_type
replace mshare_contract = mshare_contract + mshare_contract[_n+1] if firm==1 & contract_type==1 & month<201212 
replace mshare_contract = mshare_contract + mshare_contract[_n+1] if firm==2 & contract_type==1 & month<201212
* (3) drop default
drop if contract_type==2
********************************************************************************
* Merge supplier dependence on wholesale market (rough proxy, but OK since just used as IV)
merge m:1 firm using `"${PATH_OUT_DATA}/dep_wholesale"'
drop _merge
* Rescale dependence measure.
replace share_purchase=share_purchase+1
* Rescale price data on monthly level.
replace price = price / 12
sort month firm contract_type
* Construct lagged wholesale prices.
gen ws_spot_lag = wholesale_spot[_n-11]
gen ws_future_lag = wholesale_contract[_n-11]
* Construct additional price instruments based on interaction of purchase share and deviation from long-term wholesale price
sum wholesale_spot
gen ws_centered = wholesale_spot - r(mean)
gen ws_spot_lag_centered = ws_spot_lag - r(mean)
* Generate instruments based on interaction of wholesale prices and production share.
gen ws_spot_lag_int = ws_spot_lag_centered * share_purchase
gen ws_future_lag_int = ws_future_lag * share_purchase
gen ws_spot_int = ws_centered * share_purchase
gen ws_future_int = wholesale_contract * share_purchase
* Drop centered wholesale prices.
drop ws_centered ws_spot_lag_centered
* Merge green price component by supplier/contract_type
merge m:1 month firm contract_type using `"${PATH_OUT_DATA}/pricedata_component"'
drop if _merge==2
drop _merge
drop var_component fix_component
********************************************************************************
* Run price regressions to get means of price beliefs.
/*
	We are regressing monthly prices on explanatory variables that 
	consumers use to form price belief expectations.
	In particular, we assume that consumers know firm FE, the markup for green 
	electricity, and seasonal effects and time trends.
	
	We run a few specifications to see how sensitive our estimated price belief distribution is to what exactly we incorporate in the auxiliary model.
	1. Very parsimonious: no firm fixed effects.
	2. Currently our preferred specification: Firm FE, green markup term, seasonal and time trend effects
	3. Extended specification: incorporating 2. and additional statistics from wholesale market, with which most consumers are probably not familiar.
*/

* Run regressions of price on various regressors and fixed effects.
gen month_sq = month^2
gen month_of_year = real(substr(string(month),-2,.))
* Generate period dummies.
gen period = 1 if month <=201209
replace period = 2 if month > 201209 & month<=201306
replace period = 3 if month > 201306 & month<=201403
replace period = 4 if month > 201403 & month<=201412
replace period = 5 if month > 201412 & month<=201509
replace period = 6 if month > 201509 & month<=201606
* Generate dummy for VAT reform.
gen vat_dec = (month>=201404)
replace vat_dec = 2 if month>201508
* Check for evolution of price dispersion over time.
* Standard deviation per period.
* 1 : 4.125814
* 2 : 3.31435
* 3 : 3.304667
* 4 : 2.335758
* 5 : 2.335235
* 6 : 2.726935
bysort period: sum price
* Standard deviations by vat period (not used in main specification).
* 1 : 3.667797
* 2 : 2.262177
* 3 : 2.703723
bysort vat_dec: sum price
* Convert date numbes to actual date type.
tostring month, gen(month_str) format(%20.0f)
gen emonth = date(month_str,"YM")
* Plot evolution of prices for conventional contracts.
* vertical lines for period splits: 19237 - 19510 - 19783 - 20058 - 20332 
* vertical lines for vat splits: 19783 - 20301
twoway (line price emonth, sort) if contract_type==1, ytitle(Monthly price)  xline(19237 19510 19783 20058 20332) by(, title(Evolution of prices per period for conventional contracts)) by(firm)
graph save `"${PATH_OUT_FIGURES}/price_evo_conv_period"', replace
* Plot evolution of prices for green contracts.
twoway (line price emonth, sort) if contract_type==3, ytitle(Monthly price)  xline(19237 19510 19783 20058 20332) by(, title(Evolution of prices per period for green contracts)) by(firm)
graph save `"${PATH_OUT_FIGURES}/price_evo_green_period"', replace
* Model 1: 
reg price green i.month_of_year month month_sq 
predict price_pred_1, xb
predict price_eps_1, resid
gen pred_error_1 = price - price_pred_1
hist pred_error_1
graph save `"${PATH_OUT_FIGURES}/pred_error_1"', replace
* Model 2:
reg price i.firm green i.month_of_year month month_sq 
predict price_pred_2, xb
predict price_eps_2, resid
gen pred_error_2 = price - price_pred_2
hist pred_error_2
graph save `"${PATH_OUT_FIGURES}/pred_error_2"', replace
* Model 3:
reg price i.firm green gas_fl_price wholesale_spot wholesale_contract i.month_of_year month month_sq 
predict price_pred_3, xb
predict price_eps_3, resid
gen pred_error_3 = price - price_pred_3
hist pred_error_3
graph save `"${PATH_OUT_FIGURES}/pred_error_3"', replace
* Model 4: month fixed effects
reg price i.firm green gas_fl_price wholesale_spot wholesale_contract i.month
predict price_pred_4, xb
predict price_eps_4, resid
gen pred_error_4 = price - price_pred_4
hist pred_error_4
graph save `"${PATH_OUT_FIGURES}/pred_error_4"', replace
* Model 5: period fixed effects
reg price i.firm green gas_fl_price wholesale_spot wholesale_contract month month_sq i.period i.vat_dec 
predict price_pred_5, xb
predict price_eps_5, resid
gen pred_error_5 = price - price_pred_5
hist pred_error_5
graph save `"${PATH_OUT_FIGURES}/pred_error_5"', replace
* Model 6: period fixed effects but reduced set of regressors
reg price i.firm green month month_sq i.period i.vat_dec 
predict price_pred_6, xb
predict price_eps_6, resid
gen pred_error_6 = price - price_pred_6
hist pred_error_6, 
graph save `"${PATH_OUT_FIGURES}/pred_error_6"', replace
* Compare implied prediction errors.
sum pred_error*
* Based on these statistics, we go for the prediction from model 6. 
* Compare different predictions for price belief distribution.
sum price price_pred* price_eps* pred_error*
* Correlation between various predicted and actual prices.
corr price price_pred*
* Save the different price belief specifications to a separate file.
preserve
keep month firm contract_type mshare mshare_contract price price_pred* price_eps* price_pred*
save `"${PATH_OUT_DATA}/price_belief_est"', replace
restore
* Clean up for saving main data set.
* Only keep predicted price from prefered model, currently model 2.
drop period vat_dec emonth month_str month_sq month_of_year pred_error* price_pred_1-price_pred_5 price_eps*
* END OF PRICE BELIEF REGRESSIONS.
********************************************************************************

********************************************************************************
* Save as Stata file.
save `"${PATH_OUT_DATA}/new_master_data"', replace
* Export to CSV format for importing in MATLAB.
*drop if firm==7
* Main data set for importing into MATLAB.
export delimited using `"${PATH_OUT_DATA}/new_master_data.csv"', nolabel replace

