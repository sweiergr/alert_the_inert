/*
	Reshape price components data to be compatible with market share data and survey.

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
*some have missing ID
egen ID_new = concat(Supplier Name VarFix Duration), p(" - ")
order ID ID_new
**************************
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
* We drop the VARIO contracts from EDF, as they have no contract MS
drop if contract_type=="delete_noms"
*We also drop two other have no contract MS.
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
*some are strings because they contain the word "MISSING"
       foreach var of varlist _all {
                capture confirm string variable `var'
                if !_rc {
                        replace `var'="" if `var'=="MISSING"
						destring `var', replace
                }
        }
**************************

**************************
order ID_new firm contract_type Rate*
foreach var of varlist Rate1201201-Rate2201609 {
replace `var'=9999 if `var'==.
}
**************************

**************************
*reshape data 
drop Total*

reshape long Subscription Rate1 Rate2 Extra1 Extra2, i(ID_new) j(month)
drop Rate2
replace Rate1=. if Rate1==9999

rename Subscription fix_component
rename Rate1 var_component
*careful if Extra2 is empty!!
gen extra1_help=0
gen extra2_help=0
replace extra1_help=Extra1 if !(Extra1==.)
replace extra2_help=Extra2 if !(Extra2==.)
gen extra_component=extra1_help+extra2_help
replace extra_component=. if extra_component==0
drop extra1_help extra2_help
drop Extra1 Extra2 
order month firm contract_type ID_new fix_component
drop ID_new
**************************

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
*Give outside_option the conventional
replace contract_type=1 if firm==7
order month firm contract_type
*add numberlabel
numlabel, add
*************************
**************************
* Manual replacements for miscoded information, see reshape_new_pricedata
sort month firm contract_type

replace fix_component = 1 if firm==2 & contract_type==2 & month<=201204
replace var_component = 1 if firm==2 & contract_type==2 & month<=201204
replace extra_component = 1 if firm==2 & contract_type==2 & month<=201204

replace fix_component = 1 if firm==5 & contract_type==1 & month>=201404
replace var_component = 1 if firm==5 & contract_type==1 & month>=201404
replace extra_component = 1 if firm==5 & contract_type==1 & month>=201404

replace fix_component = 1 if firm==4 & contract_type==3 & month>=201606
replace var_component = 1 if firm==4 & contract_type==3 & month>=201606
replace extra_component = 1 if firm==4 & contract_type==3 & month>=201606

replace fix_component = 1 if firm==1 & contract_type==3 & month>=201604
replace var_component = 1 if firm==1 & contract_type==3 & month>=201604
replace extra_component = 1 if firm==1 & contract_type==3 & month>=201604

drop if fix_component==. & var_component==. & extra_component==.
count if fix_component==.
count if var_component==.
count if extra_component==.
**************************

*************************
*create monthly averages by firm & contract_type
sort firm contract_type month
by firm contract_type month, sort: egen fix_component_avg = mean(fix_component)
by firm contract_type month, sort: egen var_component_avg = mean(var_component)
by firm contract_type month, sort: egen extra_component_avg = mean(extra_component)

drop fix_component var_component extra_component

*drop duplicates
sort firm contract_type month
quietly by firm contract_type month:  gen dup = cond(_N==1,0,_n)
*tab dup
drop if dup>1
drop dup

*drop missing values
drop if fix_component==. & var_component==. & extra_component==.
*************************
*************************
sort firm contract_type month
replace fix_component = fix_component[_n+4]  if firm==2 & contract_type==2 & month<=201204
replace var_component = var_component[_n+4] if firm==2 & contract_type==2 & month<=201204
replace extra_component = extra_component[_n+4] if firm==2 & contract_type==2 & month<=201204
sort month firm contract_type 
replace fix_component = fix_component[_n+1] if firm==5 & contract_type==1 & month>=201404
replace var_component = var_component[_n+1] if firm==5 & contract_type==1 & month>=201404
replace extra_component = extra_component[_n+1] if firm==5 & contract_type==1 & month>=201404
replace fix_component = fix_component[_n-1] if firm==4 & contract_type==3 & month>=201606
replace var_component = var_component[_n-1] if firm==4 & contract_type==3 & month>=201606
replace extra_component = extra_component[_n-1] if firm==4 & contract_type==3 & month>=201606
replace fix_component = fix_component[_n-1]  if firm==1 & contract_type==3 & month>=201604
replace var_component = var_component[_n-1] if firm==1 & contract_type==3 & month>=201604
replace extra_component = . if firm==1 & contract_type==3 & month>=201604
*************************

*Now account for the fact that ECS has green component on all contracts
sort firm month contract_type
replace extra_component_avg=extra_component_avg[_n-1] if extra_component_avg==. & firm==1
*correct for this in the variable component
replace var_component_avg = var_component_avg - extra_component_avg if firm==1 & contract_type==3

*************************
*deflate prices
merge m:1 month using `"${PATH_OUT_DATA}/CPI_bel"'
drop if _merge==2
drop _merge
rename fix_component_avg fix_nominal
rename var_component_avg var_nominal
rename extra_component_avg extra_nominal
gen fix_component_avg=fix_nominal/CPI
gen var_component_avg=var_nominal/CPI
gen extra_component_avg=extra_nominal/CPI
drop CPI CPI_yrl fix_nominal var_nominal extra_nominal
*************************

*************************
*plot it to see differences - var and extra
preserve
collapse (mean) fix_component_avg var_component_avg extra_component_avg, by(month firm)

*prep data
sort month
gen newdate=date(string(month),"YM")
gen date=newdate
format %td date	

*graph
twoway (line var_component_avg date if firm==1, lwidth(thin)) ///
	(line extra_component_avg date if firm==1, lwidth(thin) yaxis(2)), ///
	xtitle("") ///
	ytitle("") ytitle("",axis(2)) title("ECS")  ///
	ylabel(0(2)10) ///
	ylabel(0(2)10,axis(2)) ///
	legend(label(1 variable component (real)) label(2 green component (real)) ) ///
	name(ECS,replace) 
capture window manage close graph
twoway (line var_component_avg date if firm==2 , lwidth(thin)) ///
	(line extra_component_avg date if firm==2 , lwidth(thin) yaxis(2)), ///
	xtitle("") ///
	ytitle("") ytitle("",axis(2))  title("EDF")  ///
	ylabel(0(2)10) ///
	ylabel(0(2)10,axis(2)) ///
	legend(off) name(EDF,replace)
capture window manage close graph 
twoway (line var_component_avg date if firm==3 , lwidth(thin)) ///
	(line extra_component_avg date if firm==3 , lwidth(thin) yaxis(2)), ///
	xtitle("") ///
	ytitle("") ytitle("",axis(2)) title("Eneco")  ///
	ylabel(0(2)10) ///
	ylabel(0(2)10,axis(2)) ///
	legend(off) name(Eneco,replace) 
capture window manage close graph
twoway (line var_component_avg date if firm==4 , lwidth(thin)) ///
	(line extra_component_avg date if firm==4 , lwidth(thin) yaxis(2)), ///
	xtitle("") ///
	ytitle("") ytitle("",axis(2)) title("ENI")  ///
	ylabel(0(2)10) ///
	ylabel(0(2)10,axis(2)) ///
	legend(off) name(ENI,replace) 
capture window manage close graph
twoway (line var_component_avg date if firm==5 , lwidth(thin)) ///
	(line extra_component_avg date if firm==5 , lwidth(thin) yaxis(2)), ///
	xtitle("") ///
	ytitle("") ytitle("",axis(2)) title("Essent")  ///
	ylabel(0(2)10) ///
	ylabel(0(2)10,axis(2)) ///
	legend(off) name(Essent,replace)
capture window manage close graph 
twoway (line var_component_avg date if firm==6, lwidth(thin)) ///
	(line extra_component_avg date if firm==6, lwidth(thin) yaxis(2)), ///
	xtitle("") ///
	ytitle("") ytitle("",axis(2)) title("Lampiris")  ///
	ylabel(0(2)10) ///
	ylabel(0(2)10,axis(2)) ///
	legend(off)	name(Lampiris,replace) 
capture window manage close graph
capture grc1leg ECS EDF Eneco ENI Essent Lampiris, name(g3,replace) ///
	title("Price components in EUR cent / kWh") ///
	subtitle("averages per supplier (excl. fix component)")
capture window manage close graph

restore
*************************

*drop default
drop if contract_type==2

*labelling
label var fix_component_avg "Avg fix price component, EUR/month (supplier /contract type)"
label var var_component_avg "Avg of variable price component, cent/kWh (supplier /contract type)"
label var extra_component_avg "Avg of extra price component, cent/kWh (supplier /contract type)"

*************************
* Save as Stata file.
save `"${PATH_OUT_DATA}/pricedata_component"', replace

