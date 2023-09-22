/*

This file runs reduced-form regressions on the macro data for finding 
evidence for state-dependence in our data.
Moreover, we run regressionsto provide evidence for validity of Hausman IV assumption
in our application.

*/

* BASICS.
clear
capture log close
version 13
set more off
set mem 4g
eststo clear
* Header do-File with path definitions, those end up in global macros.
include project_paths
log using `"${PATH_OUT_ANALYSIS}/log/`1'.log"', replace

*********************************
* Load survey data 
use `"${PATH_OUT_DATA}/new_master_data.dta"', clear
sort month firm contract_type
* Creating lagged market shares on firm and contract level.
gen mshare_firm_lag = mshare[_n-11]
gen mshare_contract_lag = mshare_contract[_n-11]
gen mshare_contract_forward = mshare_contract[_n+11] 
* Create indicator for informedness of population.
gen vtest_inc = vtest*incumbent
* Construct lagged cost shifters.
gen wspot_lag = wholesale_spot[_n-11]
gen gasFL_lag = gas_fl_price[_n-11]
gen gasWL_lag = gas_wl_price[_n-11]
* Scale prices in 100-EUR.
replace price = price / 100
* Generate relative ad spending per subscriber.
gen as_rel = ad_spending_norm

* Run regressions of market shares on explanatory variables.
reg mshare_contract i.firm price green ad_spending_mly
reg mshare_contract i.firm price green ad_spending_mly mshare_contract_lag
reg mshare_contract i.firm price green ad_spending_mly mshare_contract_lag vtest vtest_inc
reg mshare_contract i.firm price green ad_spending_mly mshare_firm_lag
reg mshare_contract i.firm price green ad_spending_mly mshare_firm_lag vtest vtest_inc
* Adjusted this to forwarded market shares in new data set.
reg mshare_contract_forward i.firm price green 
reg mshare_contract_forward i.firm price green mshare_contract
reg mshare_contract_forward i.firm price green mshare_contract vtest vtest_inc
reg mshare_contract_forward  incumbent price green mshare_contract
reg mshare_contract_forward  incumbent price green mshare_contract
ivreg mshare_contract_forward  incumbent (price=wholesale_spot gas_fl_price gas_wl_price) green mshare_contract 
ivreg mshare_contract_forward  incumbent (price=wholesale_spot gas_fl_price gas_wl_price) green mshare_contract 
ivreg mshare_contract_forward  incumbent (price=wholesale_spot gas_fl_price gas_wl_price) green mshare_contract ad_spending_mly 

lab var mshare_contract_forward "Market Share"
lab var mshare_contract "Lagged share"
lab var price "Price"
lab var as_rel "Advertising"
lab var incumbent "Incumbent"
lab var green "Green contract"
* Rescale advertising variable to avoid 0.000-display of coefficients.
replace as_rel = as_rel / 1000
lab var mshare_contract "Market share"
lab var mshare_contract_lag "Lagged share"
********************************************************************************
* These are the regressions currently reported in the paper.
ivregress gmm mshare_contract  incumbent (price=wholesale_spot gas_fl_price gas_wl_price) green as_rel
eststo no_sd_inc
estadd local lagged_ms "No" 
estadd local firm_fe "No"
estadd local lag_inst "No"
ivregress gmm mshare_contract  i.firm (price=wholesale_spot gas_fl_price gas_wl_price) green as_rel
eststo no_sd_fed
estadd local lagged_ms "No" 
estadd local firm_fe "Yes"
estadd local lag_inst "No"
ivregress gmm mshare_contract incumbent (price=wholesale_spot gas_fl_price gas_wl_price) green mshare_contract_lag as_rel
eststo sd_inc
estadd local lagged_ms "Yes" 
estadd local firm_fe "No"
estadd local lag_inst "No"
ivregress gmm mshare_contract  i.firm (price = wspot_lag gasFL_lag gasWL_lag wholesale_spot gas_fl_price gas_wl_price) mshare_contract_lag green  as_rel
eststo sd_fe
estadd local lagged_ms "Yes" 
estadd local firm_fe "Yes"
estadd local lag_inst "No"
ivregress gmm mshare_contract  incumbent green as_rel ///
	(mshare_contract_lag price = wspot_lag  wholesale_spot ///
	  gasWL_lag gas_wl_price  ///
	  gasFL_lag gas_fl_price )
eststo sd_inc_inst
estadd local lagged_ms "Yes" 
estadd local firm_fe "No"
estadd local lag_inst "Yes"

* Write table to tex-file.
esttab using `"${PATH_OUT_TABLES}/rfmacro.tex"', ///
keep(price incumbent green as_rel mshare_contract_lag) ///
cells(b(star fmt(3)) se(par fmt(3))) replace ///
			legend label collabels(none) ///
			eqlabels(none) /// 
star(* 0.10 ** 0.05 *** 0.01) nocons ///
notes stats(r2 N lagged_ms firm_fe lag_inst, fmt(%9.3f %9.0f) labels(`"R2"' `"Observations"' `"Lagged Share"'  `"Firm Fixed Effects"' `"Lagged Share Inst."')) ///
addnotes("\textit{Data source}: Panel (2012-2016) of contract-level market shares provided by VREG." ///
 		 "\textit{Notes:} The table summarizes results from regressing contract-level market shares" ///
 		 " on contract characteristics and lagged market shares." ///
		 "Standard errors in parenthesis.") ///
mtitles("Market share" "Market share" "Market share" "Market share" "Market share")

* Another reduced form test for state dependence. Not reported in final paper.
ivregress gmm mshare_contract  incumbent green as_rel ///
		wholesale_spot wspot_lag   ///
		(price =  gasFL_lag gasWL_lag  gas_wl_price gas_fl_price )
* Evidence on relevance of instruments.
ivreg2 mshare incumbent green (price  =  wholesale_spot gas_wl_price ws_spot_lag ws_spot_int), ffirst robust

********************************************************************************
* This is used to justify validity of Hausman IV in our application.
gen month_sq = month^2
gen year = 2012 if month<201300
replace year = 2013 if month>201300 & month <201400
replace year = 2014 if month>201400 & month<201500
replace year = 2015 if month>201500 & month<201600
replace year = 2016 if month>201600
* Merge new wholesale gas prices.
merge m:1 month using `"${PATH_OUT_DATA}/wholesaleGasPrices"'
drop if _merge==2
drop _merge
* Scale price variables for better comparability of coefficients.
replace gas_wl_price = gas_wl_price / 12
replace gas_fl_price = gas_fl_price / 12
replace elec_wl_price = elec_wl_price / 12
replace price = price * 100
sum *price*
gen reg_campaign_3 = (month>=201209 & month < 201212)
gen reg_campaign_6 = (month>=201209 & month < 201303)
gen reg_campaign_12 = (month>=201209 & month < 201306)

* Regress retail prices in Flanders on campaign dummy.
reg price i.firm wgp_NL wgp_GER i.contract_type wholesale_spot wholesale_contract month i.year reg_campaign_3, robust
reg price i.firm wgp_NL wgp_GER i.contract_type wholesale_spot wholesale_contract  month i.year reg_campaign_6, robust
reg price i.firm  wgp_NL wgp_GER i.contract_type wholesale_spot wholesale_contract  month i.year reg_campaign_12, robust
* Regress retail prices in Wallonia on campaign dummy.
reg elec_wl_price wgp_NL wgp_GER i.firm  i.contract_type wholesale_spot wholesale_contract i.year month reg_campaign_3, robust cluster(firm)
reg elec_wl_price wgp_NL wgp_GER i.firm  i.contract_type wholesale_spot wholesale_contract i.year month reg_campaign_6, robust cluster(firm)
reg elec_wl_price  wgp_NL wgp_GER i.firm  i.contract_type wholesale_spot wholesale_contract i.year month reg_campaign_12, robust cluster(firm)
* Update variable labels.
label var price "FL Elec. Ret. Price"
label var gas_fl_price "FL Gas Price"
label var wgp_NL "Natural Gas Wholesale NL"
label var wgp_GER "Natural Gas Wholesale GER"
label var wholesale_spot "Elec. Wholesale Spot"
label var wholesale_contract "Elec. Wholesale Fut."
label var reg_campaign_3 "Post-Regulator Campaign 3"
label var reg_campaign_6 "Post-Regulator Campaign 6"
label var reg_campaign_12 "Post-Regulator Campaign 12"
* Regress retail gas prices in Wallonia on campaign dummy.
eststo clear
eststo: reg gas_fl_price i.firm wgp_NL wgp_GER i.year month wholesale_spot wholesale_contract reg_campaign_3  if contract_type==3, robust
eststo: reg gas_fl_price i.firm wgp_NL wgp_GER i.year month wholesale_spot wholesale_contract reg_campaign_6 if contract_type==1, robust 
eststo: reg gas_fl_price i.firm wgp_NL wgp_GER i.year month wholesale_spot wholesale_contract reg_campaign_12  if contract_type==1, robust
eststo: reg gas_fl_price i.firm wgp_NL price wgp_GER i.year month wholesale_spot wholesale_contract reg_campaign_6 if contract_type==1, robust 

eststo: reg gas_wl_price i.firm wgp_NL wgp_GER i.year month wholesale_spot wholesale_contract reg_campaign_3  if contract_type==3, robust
eststo: reg gas_wl_price i.firm wgp_NL wgp_GER i.year month wholesale_spot wholesale_contract reg_campaign_6 if contract_type==1, robust 
eststo: reg gas_wl_price i.firm wgp_NL wgp_GER i.year month wholesale_spot wholesale_contract reg_campaign_12  if contract_type==1, robust
eststo: reg gas_wl_price i.firm wgp_NL price wgp_GER i.year month wholesale_spot wholesale_contract reg_campaign_6 if contract_type==1, robust 
esttab, star(* 0.10 ** 0.05 *** 0.01) b(3) se(3) compress ///
label stats(N , fmt(%9.0f) labels(`"Observations"'))
* LaTeX-table to show evidence for validity of Hausman IV in our application.
esttab using `"${PATH_OUT_TABLES}/hausman_IV_evidence_1.tex"', ///
star(* 0.10 ** 0.05 *** 0.01) ///
keep(wgp_NL wgp_GER wholesale_spot wholesale_contract reg_campaign_3 reg_campaign_6 reg_campaign_12 price) ///
cells(b(star fmt(3)) se(par fmt(3))) replace ///
			legend label collabels(none) ///
			stats(N, fmt(0 2) label(Observations) ) ///
			eqlabels(none) /// 
addnotes("\textit{Notes:} The table summarizes regression results to provide supporting evidence for " ///
		 " the validity of the Hausman IV assumption. ") ///
mtitles("Gas FL" "Gas FL" "Gas FL" "Gas FL" "Gas WL" "Gas WL"  "Gas WL"  "Gas WL")
eststo clear
* Regress retail prices in Flanders on campaign dummy.
eststo: reg price i.firm wgp_NL wgp_GER i.contract_type wholesale_spot  month i.year reg_campaign_3, robust
* Regress retail prices in Wallonia on campaign dummy.
eststo: reg elec_wl_price wgp_NL wgp_GER i.firm  i.contract_type wholesale_spot  i.year month reg_campaign_3, robust cluster(firm)
eststo: reg gas_wl_price i.firm wgp_NL wgp_GER i.year month wholesale_spot reg_campaign_3  if contract_type==1, robust
eststo: reg gas_wl_price i.firm wgp_NL wgp_GER i.year month wholesale_spot  reg_campaign_6 if contract_type==1, robust 
eststo: reg gas_wl_price i.firm wgp_NL wgp_GER i.year month wholesale_spot  reg_campaign_12  if contract_type==1, robust
esttab, star(* 0.10 ** 0.05 *** 0.01) b(3) se(3) stats() compress /// label stats(N , fmt(%9.0f) labels(`"Observations"'))

esttab using `"${PATH_OUT_TABLES}/hausman_IV_evidence_2.tex"', ///
star(* 0.10 ** 0.05 *** 0.01) ///
keep(wgp_NL wgp_GER wholesale_spot reg_campaign_3 reg_campaign_6 reg_campaign_12) ///
cells(b(star fmt(3)) se(par fmt(3))) replace ///
			legend label collabels(none) /// stats(N, fmt(0 2) label(Observations) ) ///
			stats() ///
			eqlabels(none) /// 
addnotes("\textit{Notes:} The table summarizes regression results to provide supporting evidence for " ///
		 " the validity of our Hausman IV. ") ///
mtitles("Elec FL" "Elec WL" "Gas WL" "Gas WL"  "Gas WL")
