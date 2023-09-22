/*

This file runs reduced-form regressions to support our identification assumptions.
Specifically, we analyze the relationship between firm advertising, PCW usage, and the regulator campaign.

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
* Load monthly internet-advertising data. 
use `"${PATH_OUT_DATA}/internet_adv.dta"', clear

* Merge with data on average advertising and prices.
* Merge average prices.
merge 1:1 month using `"${PATH_OUT_DATA}/avg_price_data.dta"'
drop _merge
* Merge average advertising.
merge 1:1 month using `"${PATH_OUT_DATA}/avg_adv_data.dta"'
drop _merge
* Scale advertising in millions.
replace aggregate_adv_raw = aggregate_adv_raw / 1000000
replace aggregate_adv = aggregate_adv / 1000000
replace ad_spending_mly = ad_spending_mly /1000000

* Generate January dummies.
gen month_of_year = real(substr(string(month),-2,.))
gen dummy_jan = 0
replace dummy_jan = 1 if month_of_year==1
* Generate quarter of year dummies.
gen qu_of_year = 1 if month_of_year<=3
replace qu_of_year = 2 if month_of_year>3 & month_of_year<=6
replace qu_of_year = 3 if month_of_year>6 & month_of_year<=9
replace qu_of_year = 4 if month_of_year>9

* Generate measure ofprice variance within a year.
gen year = 1 if month<201301
replace year = 2 if month<201401 & year==.
replace year = 3 if month<201501 & year==.
replace year = 4 if month<201601 & year==.
replace year = 5 if month>=201601 & year==.


** here we want to show that PCW usage is negatively affected by advertising
** because it decreases benefits from search. 
** Moreover, we want to see that PCW usage also is positively correlated with our 
** search cost shifters, i.e. internet penetration and regulator campaign.
* Regress PCW usage on search cost shifters and aggregate advertising.
* What do we expect?
* - broandband_internet coefficient is positive.
* - advertising coefficient is negative
* - regulator campaign positive
* - price coefficient (mean): positive?
* - price coefficient (variance): positive?
reg vtest month i.month_of_year price_sd_year price fixed_broadband regulator ad_spending_mly, robust
estadd local firm_fe "No"
estadd local month_year_fe "Yes"
eststo ada_vtest_shifters

* These graphs are not use din the paper anymore since we only rely on regressions to get a clearer picture.
* Graphs capture too manby different things.
* Plot relationship between PCW usage and ad spending: ideally see a negative relationship here.
lowess vtest ad_spending_mly, ytitle(PCW Usage) xtitle(Aggregate Monthly Advertising Expenditure) title(PCW usage vs. Ad Spending)
graph export `"${PATH_OUT_FIGURES}/ada_evidence_vtestadspending.pdf"', replace
* Plot relationship between internet penetration and PCW usage, ideally see positive relationship here.
drop if regulator==1
drop if month<=201302
lowess vtest fixed_broadband, ytitle(PCW Usage) xtitle(Internet Penetration Rate) title(PCW usage vs. Internet Penetration)
graph export `"${PATH_OUT_FIGURES}/ada_evidence_vtestinternet.pdf"', replace

* Save LOWESS values to format graph nicely in Python.
lowess vtest fixed_broadband, nograph gen(vt_int_fit)
lowess vtest ad_spending_mly, nograph gen(vt_ads_fit)
* Export data for graphs in csv format.
export delimited vtest ad_spending_mly fixed_broadband vt_int_fit vt_ads_fit using `"${PATH_OUT_ANALYSIS}/ada_evidence_data.csv"', nolabel replace
lowess ad_spending_mly fixed_broadband, ytitle(Ad Spending) xtitle(Internet Penetration Rate) title(Ad Spending usage vs. Internet Penetration)
graph export `"${PATH_OUT_FIGURES}/ada_evidence_adspendinginternet.pdf"', replace


* Load monthly firm-level master data to analyze firm-specific advertising. 
use `"${PATH_OUT_DATA}/new_master_data.dta"', clear
* Generate month of year dummies.
gen month_of_year = real(substr(string(month),-2,.))
gen dummy_jan = 0
replace dummy_jan = 1 if month_of_year==1
* Collapse data to firm-level.
collapse mshare price wholesale_spot vtest ad_spending_mly ad_spending_norm share_purchase month_of_year dummy_jan, by(month firm)
* Rescale advertising expenditures.
replace ad_spending_mly = ad_spending_mly /1000000
* Rename price variable before merging average retail price in a given month.
rename price  price_firm
* Merge with data on average advertising and prices.
* Merge average prices.
merge m:1 month using `"${PATH_OUT_DATA}/avg_price_data.dta"'
drop _merge
* Merge average advertising.
merge m:1 month using `"${PATH_OUT_DATA}/avg_adv_data.dta"'
drop _merge
* Merge with data on internet and regulator campaign.
merge m:1 month using `"${PATH_OUT_DATA}/internet_adv.dta"'
drop _merge
* Label variables for tables.
label var ad_spending_mly "Ad spending"
label var wholesale_spot "Wholesale price"
label var regulator "Regulator campaign"
label var fixed_broadband "Internet penetration"
label var price "Retail price (mean)"
label var price_sd_year "Retail price (SD)"
label var vtest "PCW Usage"
* Regress aggregate advertising on our assumed search cost shifters.
reg ad_spending_mly i.firm regulator fixed_broadband i.month_of_year wholesale_spot price, robust cluster(firm)
estadd local firm_fe "Yes"
estadd local month_year_fe "Yes"
// estadd local lag_inst "No"
eststo ada_adspending_internet
esttab using `"${PATH_OUT_TABLES}/rf_ada_evidence.tex"', replace ///
keep(regulator fixed_broadband wholesale_spot price price_sd_year ad_spending_mly) ///
star(* 0.10 ** 0.05 *** 0.01) ///
notes stats(r2 N firm_fe month_year_fe , fmt(%9.3f %9.0f) labels(`"R2"' `"Observations"' `"Firm FE"' `"Month-of-Year FE"')) ///
cells(b(star fmt(3)) se(par fmt(3))) ///
			legend collabels(none) ///
			eqlabels(none) label /// 
addnotes("\textit{Data source}: VREG surveys 2012-2016 and Nielsen MDB." ///
 		 "\textit{Notes:} Column (1) summarizes the results from an OLS regression of the" ///
 		 "monthly aggregate share of PCW users on various potential shifters of the" ///
 		 "expected benefits and costs of search. Column (2) summarizes the results" ///
 		 "from an OLS regression of monthly firm-specific advertising on a similar set" ///
 		 "of shifters. All regressions include month-of-year fixed effects. Column (2)" ///
 		 "also incorporates firm fixed effects. Robust standard errors in parentheses.") ///
mtitles("PCW usage" "Ad spending")