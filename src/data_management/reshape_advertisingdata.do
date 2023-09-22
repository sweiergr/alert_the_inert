/*
	Reshape advertising data to be compatible with market share data.
	(1) using Nielsen data (monthly, 2011-16)
		- version1 with outside good regrouping only those suppliers we know
		- version2 with outside good regrouping all other spending in energy
	(2) using Nielsen data (yearly, 2011-15)
	(3) using UBA data (yearly, 2013-2015)
	
	Transform nominal into real prices (base: Jan 2012)
	Construct moving average in (1) version1 and version2
	Normalize wrt number of consumers (elec+gas)
	
*/

clear
capture log close
version 13
set more off
set mem 4g
* Header do-File with path definitions, those end up in global macros.
include project_paths
log using `"${PATH_OUT_DATA}/log/`1'.log"', replace

******************************************************************************
* NIELSEN DATA - monthly (2011-2016)
use `"${PATH_OUT_DATA}/adv_nielsen"', clear
* clean
rename Soussecteur subsector
rename Groupedannonceurs group
rename Annonceur advertiser
rename Marque brand
rename SubMedia mediatype
rename Annee year
rename Anneemois month
rename InvestissementtotalSommes ad_spending_brand
label var ad_spending_brand "Monthly media spending by brand in Belgium"
drop group advertiser
*we are not interested in other combustibles: for now just keep elec and gas
tab subsector
drop if subsector=="COMBUSTIBLES LIQUIDES - MAZOUT" | subsector=="COMBUSTIBLES SOLIDES - AUTRES"
*destring month
replace month = subinstr(month, "-", "",.) 
destring month, replace
order month brand
sort month
preserve
***********************************************
*Alternative definition of outside option -> only those suppliers we know
drop subsector year
*assign ad_spending to firm
gen firm=.
replace firm=1 if brand=="ELECTRABEL" | brand=="ENGIE"
replace firm=2 if brand=="EDF" | brand=="EDF LUMINUS" | brand=="LUMINUS" | brand=="LUMINUS.BE"
replace firm=3 if brand=="ENECO"
replace firm=4 if brand=="ENI" | brand=="NUON"
replace firm=5 if brand=="ESSENT"
replace firm=6 if brand=="LAMPIRIS"
*these 4 are spelled out in the "macro market shares"
replace firm=7 if brand=="BELPOWER" | brand=="EBEM" | brand=="ELEGANT" | brand=="OCTA +" 
*these is additionally spelled out in the "contract market shares" or price data
replace firm=7 if brand=="POWEO" | brand=="COMFORT ENERGY" | brand=="ELEXYS"
*missing in advertisement but spelled out in the price data:
*aspiravi, energy people, mega, klinkenberg
drop if firm==.
label define firml 1 "ECS" 2 "EDF" 3 "Eneco" 4 "ENINuon" 5 "Essent" 6 "Lampiris" 7 "Other"
label values firm firml
numlabel, add
*construct ad_spending by firm
collapse (sum) ad_spending_brand, by(month firm)
rename ad_spending_brand ad_spending_mly
label var ad_spending_mly "Monthly (nominal) media spending (in EUR) in Belgium"
*deflate prices
merge m:1 month using `"${PATH_OUT_DATA}/CPI_bel"'
drop if _merge==2
drop _merge
rename ad_spending_mly ad_spending_mly_nominal
gen ad_spending_mly=ad_spending_mly_nominal/CPI
label var ad_spending_mly "Monthly real media spending (in EUR) in Belgium"
drop CPI CPI_yrl ad_spending_mly_nominal
*add zeros if no ad_spending was observed
reshape wide ad_spending_mly, i(month) j(firm)
forvalues i = 1(1)7{
		replace ad_spending_mly`i'=0 if ad_spending_mly`i'==.
}
reshape long ad_spending_mly, i(month) j(firm)
*plot spending over time
gen yr=int(month/100)
gen m=mod(month,100)
gen date=ym(y,m)
label var date "Jan 2011 - June 2016"
drop yr m date

*save the data to construct advertising spending table
tempfile adv_table_total
save `adv_table_total'

*normalize by number of costumers
merge m:1 month firm using `"${PATH_OUT_DATA}/nb_customer"'
drop if _merge==2
drop _merge

gen ad_spending_norm=ad_spending_mly/nb_both
label var ad_spending_norm "Euro spent by customer, monthly"
drop nb_elec nb_gas nb_both
*normalize by number of non-costumers
gen ad_spending_noncust_norm=ad_spending_mly/nb_nonboth
label var ad_spending_noncust_norm "Euro spent by non-customer, monthly"
drop nb_nonelec nb_nongas nb_nonboth
*Create moving average
*(1) no leads, 6 lags
reshape wide ad_spending_mly ad_spending_norm ad_spending_noncust_norm, i(month) j(firm)
forvalues i = 1(1)7{
	tsset month
	tssmooth ma ma_adv_exp`i'=ad_spending_mly`i', window(6 1 0)
	tssmooth ma ma_adv_norm`i'=ad_spending_norm`i', window(6 1 0)
	tssmooth ma ma_adv_noncust_norm`i'=ad_spending_noncust_norm`i', window(6 1 0)
	replace ad_spending_mly`i' = ma_adv_exp`i'
	replace ad_spending_norm`i' = ma_adv_norm`i'
	replace ad_spending_noncust_norm`i' = ma_adv_noncust_norm`i'
	drop ma_adv_exp`i' ma_adv_norm`i' ma_adv_noncust_norm`i'
}
reshape long ad_spending_mly ad_spending_norm ad_spending_noncust_norm, i(month) j(firm)
sort firm month
label var ad_spending_mly "Monthly real media spending (in EUR) in Belgium, moving average"
label var ad_spending_norm "Euro spent by customer in Belgium, monthly (moving average)"
label var ad_spending_noncust_norm "Euro spent by non-customer in Belgium, monthly (moving average)"

*drop months that we do not need in the analysis
drop if month<201201

*plot moving average time
gen yr=int(month/100)
gen m=mod(month,100)
gen date=ym(y,m)
label var date "Jan 2012 - June 2016"
*xtline ad_spending_mly, t(date) i(firm) overlay tlabel(, format(%tm))
drop yr m date

*plot nb spent by customer time
gen yr=int(month/100)
gen m=mod(month,100)
gen date=ym(y,m)
label var date "Jan 2012 - June 2016"
*xtline ad_spending_norm, t(date) i(firm) overlay tlabel(, format(%tm))
drop yr m date

*plot nb spent by non-customer time
gen yr=int(month/100)
gen m=mod(month,100)
gen date=ym(y,m)
label var date "Jan 2012 - June 2016"
*xtline ad_spending_noncust_norm, t(date) i(firm) overlay tlabel(, format(%tm))
drop yr m date

* Save as Stata file.
sort month firm
save `"${PATH_OUT_DATA}/adv_nielsen_long_outside"', replace
* Add export to csv for importing data into MATLAB.
export delimited using `"${PATH_OUT_DATA}/advertisingdata.csv"', replace nolabel
***********************************************

****************
*Table on total adv spending
use `adv_table_total', clear

*calculate yearly spending per firm for table in datasection
gen year=int(month/100)
sort firm year
collapse (rawsum) ad_spending_mly, by(year firm)
gen fi = round(ad_spending_mly/1000000, 0.01)
drop ad_spending_mly
drop if year==2011

reshape wide fi, i(year) j(firm)

rename fi1 ECS
rename fi2 EDF
rename fi3 Eneco
rename fi4 Eni
rename fi5 Essent
rename fi6 Lampiris
rename fi7 Other

eststo clear
by year: eststo: estpost summarize ///
    ECS EDF Eneco Eni Essent Lampiris Other
esttab using `"${PATH_OUT_TABLES}/adspending.tex"', replace ///
	mtitle("2012" "2013" "2014" "2015" "2016") ///
	cells("sum") booktabs noobs nonum nolabel nodepvar collabels(none)
****************

***********************************************
restore

*assign ad_spending to firm (put all others into "OTHER" but drop gas and other combustibles
drop if !(subsector=="ENERGIE")
drop subsector year

gen firm=7 //Ohter
replace firm=1 if brand=="ELECTRABEL" | brand=="ENGIE"
replace firm=2 if brand=="EDF" | brand=="EDF LUMINUS" | brand=="LUMINUS" | brand=="LUMINUS.BE"
replace firm=3 if brand=="ENECO"
replace firm=4 if brand=="ENI"  | brand=="NUON"
replace firm=5 if brand=="ESSENT"
replace firm=6 if brand=="LAMPIRIS"

label define firml 1 "ECS" 2 "EDF" 3 "Eneco" 4 "ENINuon" 5 "Essent" 6 "Lampiris" 7 "Other"
label values firm firml

numlabel, add

*construct ad_spending by firm
collapse (sum) ad_spending_brand, by(month firm)
rename ad_spending_brand ad_spending_mly
label var ad_spending_mly "Monthly (nominal) media spending (in EUR) in Belgium"

*deflate prices
merge m:1 month using `"${PATH_OUT_DATA}/CPI_bel"'
drop if _merge==2
drop _merge
rename ad_spending_mly ad_spending_mly_nominal
gen ad_spending_mly=ad_spending_mly_nominal/CPI
label var ad_spending_mly "Monthly real media spending (in EUR) in Belgium"
drop CPI CPI_yrl ad_spending_mly_nominal

*add zeros if no ad_spending was observed
reshape wide ad_spending_mly, i(month) j(firm)
forvalues i = 1(1)7{
		replace ad_spending_mly`i'=0 if ad_spending_mly`i'==.
}
reshape long ad_spending_mly, i(month) j(firm)

*plot spending over time
gen yr=int(month/100)
gen m=mod(month,100)
gen date=ym(y,m)
drop if yr<2012
label var date "Jan 2012 - June 2016"
*xtline ad_spending_mly, t(date) i(firm) overlay tlabel(, format(%tm))
drop yr m date

*normalize by number of costumers
merge m:1 month firm using `"${PATH_OUT_DATA}/nb_customer"'
drop if _merge==2
drop _merge

gen ad_spending_norm=ad_spending_mly/nb_both
label var ad_spending_norm "Euro spent by customer, monthly"
drop nb_elec nb_gas nb_both

*normalize by number of non-costumers
gen ad_spending_noncust_norm=ad_spending_mly/nb_nonboth
label var ad_spending_noncust_norm "Euro spent by non-customer, monthly"
drop nb_nonelec nb_nongas nb_nonboth

*Create moving average
*(1) no leads, 6 lags
reshape wide ad_spending_mly ad_spending_norm ad_spending_noncust_norm, i(month) j(firm)
forvalues i = 1(1)7{
	tsset month
	tssmooth ma ma_adv_exp`i'=ad_spending_mly`i', window(6 1 0)
	tssmooth ma ma_adv_norm`i'=ad_spending_norm`i', window(6 1 0)
	tssmooth ma ma_adv_noncust_norm`i'=ad_spending_noncust_norm`i', window(6 1 0)
	replace ad_spending_mly`i' = ma_adv_exp`i'
	replace ad_spending_norm`i' = ma_adv_norm`i'
	replace ad_spending_noncust_norm`i' = ma_adv_noncust_norm`i'
	drop ma_adv_exp`i' ma_adv_norm`i' ma_adv_noncust_norm`i'
}
reshape long ad_spending_mly ad_spending_norm ad_spending_noncust_norm, i(month) j(firm)
sort firm month
label var ad_spending_mly "Monthly real media spending (in EUR) in Belgium, moving average"
label var ad_spending_norm "Euro spent by customer in Belgium, monthly (moving average)"
label var ad_spending_noncust_norm "Euro spent by non-customer in Belgium, monthly (moving average)"

*drop months that we do not need in the analysis
drop if month<201201

*plot moving average time
gen yr=int(month/100)
gen m=mod(month,100)
gen date=ym(y,m)
label var date "Jan 2012 - June 2016"
*xtline ad_spending_mly, t(date) i(firm) overlay tlabel(, format(%tm))
drop yr m date

*plot nb spent by customer time
gen yr=int(month/100)
gen m=mod(month,100)
gen date=ym(y,m)
label var date "Jan 2012 - June 2016"
*xtline ad_spending_norm, t(date) i(firm) overlay tlabel(, format(%tm))
drop yr m date

*plot nb spent by non_customer time
gen yr=int(month/100)
gen m=mod(month,100)
gen date=ym(y,m)
label var date "Jan 2012 - June 2016"
*xtline ad_spending_noncust_norm, t(date) i(firm) overlay tlabel(, format(%tm))
drop yr m date
*
sort month firm

* Save as Stata file.
save `"${PATH_OUT_DATA}/adv_nielsen_long"', replace
******************************************************************************


*ALTERNATIVE ADVERTISING MEASURES
******************************************************************************
* NIELSEN DATA - yearly (2011-2015)
use `"${PATH_OUT_DATA}/adv_nielsen_yrly"', clear
* drop 
rename Soussecteur subsector
rename Groupedannonceurs group
rename Annonceur advertiser
rename Marque brand
rename SubMedia mediatype
rename Annee year
rename InvestissementtotalSommes ad_spending_brand
label var ad_spending_brand "Yearly media spending by brand in Belgium"
*do we need all these variables. Do they indicate the same?
gen test=0
replace test=1 if group==advertiser
*browse if test==0
*only GDF SUEZ GROUP seems to have different announcers
* Electrabel, GDF Suez Belgium, Suez environment, tractebel engineering
drop test group

gen test=0
replace test=1 if advertiser==brand
*sort advertiser
*browse if test==0
*Only EDF Luminus has different brands Luminus - Luminus.be
drop test advertiser

*we are not interested in gas etc
drop if !(subsector=="ENERGIE")
drop subsector

*destring year
destring year, replace
order year brand
sort year

*for now drop year==2011
drop if year==2011

preserve
***********************************************
*alternative definition of outside option -> only those suppliers we know
*assign ad_spending to firm
gen firm=.
replace firm=1 if brand=="ELECTRABEL" | brand=="ENGIE"
replace firm=2 if brand=="EDF" | brand=="EDF LUMINUS" | brand=="LUMINUS" | brand=="LUMINUS.BE"
replace firm=3 if brand=="ENECO"
replace firm=4 if brand=="ENI"  | brand=="NUON"
replace firm=5 if brand=="ESSENT"
replace firm=6 if brand=="LAMPIRIS"
replace firm=7 if brand=="BELPOWER" | brand=="EBEM" | brand=="ELEGANT" | brand=="OCTA +" 
*these are additionally spelled out in the "contract market shares"
replace firm=7 if brand=="POWEO" | brand=="POWECO" 
*missing: Aspiravi

drop if firm==.

label define firml 1 "ECS" 2 "EDF" 3 "Eneco" 4 "ENINuon" 5 "Essent" 6 "Lampiris" 7 "Other"
label values firm firml

numlabel, add

*construct ad_spending by firm
collapse (sum) ad_spending_brand, by(year firm)
rename ad_spending_brand ad_spending_yrl
label var ad_spending_yrl "Yearly (nominal) media spending (in EUR) in Belgium"

*deflate prices
gen month=year*100+1
merge m:1 month using `"${PATH_OUT_DATA}/CPI_bel"'
drop if _merge==2
drop _merge CPI month
rename ad_spending_yrl ad_spending_yrl_nominal
gen ad_spending_yrl=ad_spending_yrl_nominal/CPI_yrl
label var ad_spending_yrl "Yearly real media spending (in EUR) in Belgium"
drop CPI CPI_yrl ad_spending_yrl_nominal

*add zeros if no ad_spending was observed
reshape wide ad_spending_yrl, i(year) j(firm)
forvalues i = 1(1)7{
		replace ad_spending_yrl`i'=0 if ad_spending_yrl`i'==.
}
reshape long ad_spending_yrl, i(year) j(firm)

*TODO: construct monthly spendings: FOR NOW assume a constant share per month
gen ad_spending_mly=ad_spending_yrl/12
label var ad_spending_mly "Monthly media spending - ASSUMES CONSTANT SHARE"

* Save as Stata file.
sort year firm
save `"${PATH_OUT_DATA}/adv_nielsen_long_yrly_outside"', replace
***********************************************

***********************************************
restore

*assign ad_spending to firm
gen firm=7 //Ohter
replace firm=1 if brand=="ELECTRABEL" | brand=="ENGIE"
replace firm=2 if brand=="EDF" | brand=="EDF LUMINUS" | brand=="LUMINUS" | brand=="LUMINUS.BE"
replace firm=3 if brand=="ENECO"
replace firm=4 if brand=="ENI" | brand=="NUON"
replace firm=5 if brand=="ESSENT"
replace firm=6 if brand=="LAMPIRIS"

label define firml 1 "ECS" 2 "EDF" 3 "Eneco" 4 "ENINuon" 5 "Essent" 6 "Lampiris" 7 "Other"
label values firm firml

numlabel, add

*construct ad_spending by firm
collapse (sum) ad_spending_brand, by(year firm)
rename ad_spending_brand ad_spending_yrl
label var ad_spending_yrl "Yearly (nominal) media spending (in EUR) in Belgium"

*deflate prices
gen month=year*100+1
merge m:1 month using `"${PATH_OUT_DATA}/CPI_bel"'
drop if _merge==2
drop _merge CPI month
rename ad_spending_yrl ad_spending_yrl_nominal
gen ad_spending_yrl=ad_spending_yrl_nominal/CPI_yrl
label var ad_spending_yrl "Yearly real media spending (in EUR) in Belgium"
drop CPI CPI_yrl ad_spending_yrl_nominal

*add zeros if no ad_spending was observed
reshape wide ad_spending_yrl, i(year) j(firm)
forvalues i = 1(1)7{
		replace ad_spending_yrl`i'=0 if ad_spending_yrl`i'==.
}
reshape long ad_spending_yrl, i(year) j(firm)

*TODO: construct monthly spendings: FOR NOW assume a constant share per month
gen ad_spending_mly=ad_spending_yrl/12
label var ad_spending_mly "Monthly real media spending - ASSUMES CONSTANT SHARE"

* Save as Stata file.
sort year firm
save `"${PATH_OUT_DATA}/adv_nielsen_long_yrly"', replace
******************************************************************************


******************************************************************************
* UBA DATA

* Import data set from Excel files in src/original_data folder.
use `"${PATH_OUT_DATA}/advertisement_bel"', clear

* Rename some variables.
rename group firm

* Drop trash.
drop if firm==""

* Summarize small firms in outside good.
replace firm="Other" if firm=="BELPOWER" | firm=="EBEM" | firm=="ELEGANT" | firm=="OCTAPLUS" | firm=="WATZ" | firm=="all_other"
* ENINuon is one firm
replace firm="ENINuon" if firm=="ENI" | firm=="Nuon" | firm=="NUON"
* Electrabel is ECS
replace firm="ECS" if firm=="ELECTRABEL"
* EDF is Luminus.
replace firm="EDF" if firm=="LUMINUS"

*construct a yearly sum of spending of outside good
collapse (sum) ad_spending, by(firm year)
label var ad_spending "Yearly media spending by group in Belgium (2012 imputed)"

*deflate prices
gen month=year*100+1
merge m:1 month using `"${PATH_OUT_DATA}/CPI_bel"'
drop if _merge==2
drop _merge CPI month
rename ad_spending ad_spending_nominal
gen ad_spending=ad_spending_nominal/CPI_yrl
label var ad_spending "Yearly real media spending by group in Belgium (2012 imputed)"
drop CPI ad_spending_nominal

* Encode firm variable and construct outside good analgously to market shares.
encode firm, gen(firm_num)
drop firm
rename firm_num firm
label define firml 1 "ECS" 2 "EDF" 3 "Eneco" 4 "ENINuon" 5 "Essent" 6 "Lampiris" 7 "Other"
label values firm firml

*add numberlabel
numlabel , add

sort firm year
order firm year

* Save as Stata file.
save `"${PATH_OUT_DATA}/advertisement_bel_long"', replace
