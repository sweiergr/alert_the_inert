/*
	Reshape contract market share data and perpare it to merge with price data.

*/

clear
capture log close
version 13
set more off
set mem 4g
* Header do-File with path definitions, those end up in global macros.
include project_paths
log using `"${PATH_OUT_DATA}/log/`1'.log"', replace
* Import data set from Excel files in src/original_data folder.
use `"${PATH_OUT_DATA}/contract_acspts"', clear

**************************
*destring year variable
gen year = substr(quarter,1,4)
destring year, replace
gen quarter2 = substr(quarter,-1,.)
destring quarter2, replace
drop quarter
gen quarter = year*10+quarter2
label var quarter "Year-quarter (YYYYq)"
drop year quarter2
order quarter firm contract_type contract_acspts
**************************

**************************
*Prepare firms
replace firm="Other" if contract_type=="outsideoption" | contract_type=="outside_social"
replace firm="DNB" if contract_type=="dnb" | contract_type=="dnb_social"
replace firm="ENINuon" if firm=="EniNuon"
*Prepare contract type of outsideoption
replace contract_type="conventional" if firm=="Other"  
**************************

**************************
*Clean dataset from social and dnb
*For now, do not keep the social contracts and the network operator contracts
drop if contract_type=="social" | contract_type=="outside_social" | contract_type=="dnb_social"
drop if contract_type=="dnb"
*others
drop if contract_type=="delete" //this is the EDF vario group
**************************

**************************
*Merge some contracts
order quarter firm new_id
sort quarter firm contract_type new_id

*ECS: basis_var
replace contract_acspts = contract_acspts+contract_acspts[_n+1]+contract_acspts[_n+2] if new_id=="ecs_basis_var"
drop if new_id=="merge_ecs_basis_var"

*EDF: standard
replace contract_acspts = contract_acspts+contract_acspts[_n+1] if new_id=="edf_standard"
drop if new_id=="merge_edf_standard"

*ENI: weekendplus
replace contract_acspts = contract_acspts+contract_acspts[_n+1] if new_id=="eni_weekendplus_1yr" & quarter<=20123
drop if new_id=="merge_eni_weekendplus"

*Essent: green_fix_3yr
replace new_id="essent_green_fix_3yr_merge" if new_id=="merge_essent_green_fix_3yr"
sort quarter firm contract_type new_id
replace contract_acspts = contract_acspts+contract_acspts[_n+1] if new_id=="essent_green_fix_3yr" & quarter>=20122
drop if new_id=="essent_green_fix_3yr_merge"

*Essent: green_fix_1yr
replace new_id="essent_green_fix_1yr_merge" if new_id=="merge_essent_green_fix_1yr"
sort quarter firm contract_type new_id
replace contract_acspts = contract_acspts+contract_acspts[_n+1] if new_id=="essent_green_fix_1yr" & quarter>=20122
drop if new_id=="essent_green_fix_1yr_merge"

*Lampiris: 
replace new_id="lampiris_fix_merge" if new_id=="merge_lampiris_fix"
sort quarter firm contract_type new_id
replace contract_acspts = contract_acspts+contract_acspts[_n+1]+contract_acspts[_n+2] if new_id=="lampiris_fix"
drop if new_id=="lampiris_fix_merge"
**************************

**************************
*for now drop the variables that we do not need
keep quarter firm contract_type contract_acspts
**************************

**************************
* encode firm
gen firm_help=.
replace firm_help=1 if firm=="ECS"
replace firm_help=2 if firm=="EDF"
replace firm_help=3 if firm=="Eneco"
replace firm_help=4 if firm=="ENINuon"
replace firm_help=5 if firm=="Essent"
replace firm_help=6 if firm=="Lampiris"
replace firm_help=7 if firm=="Other"
count if firm_help==.
drop firm
rename firm_help firm
label define firml 1 "ECS" 2 "EDF" 3 "Eneco" 4 "ENINuon" 5 "Essent" 6 "Lampiris" 7 "Other"
label values firm firml
numlabel, add
order quarter firm contract_type
**************************

**************************
* encode contract_type
gen contract_type_help=.
replace contract_type_help=1 if contract_type=="conventional"
replace contract_type_help=2 if contract_type=="default"
replace contract_type_help=3 if contract_type=="green"
count if contract_type_help==.
drop contract_type
rename contract_type_help contract_type
label define contract_typel 1 "conventional" 2 "default" 3 "green"
label values contract_type contract_typel
numlabel, add
order quarter firm contract_type
********************************************************************************
********************************************************************************
*create market share for each contract type (excl social tariffs and dnp)
* Create a within_firm market shares by contract_type -> to be merged to masterdata 
*denominator
by quarter firm, sort: egen firm_acspts=sum(contract_acspts)
label var firm_acspts "firm acces points by quarter (excl. social and dnb)"

*numerator
by quarter firm contract_type, sort: egen contract_type_acspts=sum(contract_acspts)
label var contract_type_acspts "number of acces points by contract_type, firm and quarter"

*drop duplicates
quietly by quarter firm contract_type, sort:  gen dup = cond(_N==1,0,_n)
drop if dup>1
drop dup 

*generate share
gen contract_ms = contract_type_acspts/firm_acspts
label var contract_ms "contract market share by firm and quarter (excl. social and dnb)"
drop if contract_type==2 & quarter>=20131

********************************************************************************

*clean
keep quarter firm contract_type contract_ms

*prepare time variable for merge
reshape wide contract_ms, i(firm contract_type) j(quarter)

*2012-2015
foreach i of numlist 2/5{
	gen contract_ms201`i'01=contract_ms201`i'1
	gen contract_ms201`i'02=contract_ms201`i'1
	gen contract_ms201`i'03=contract_ms201`i'1
	drop contract_ms201`i'1
	gen contract_ms201`i'04=contract_ms201`i'2
	gen contract_ms201`i'05=contract_ms201`i'2
	gen contract_ms201`i'06=contract_ms201`i'2
	drop contract_ms201`i'2
	gen contract_ms201`i'07=contract_ms201`i'3
	gen contract_ms201`i'08=contract_ms201`i'3
	gen contract_ms201`i'09=contract_ms201`i'3
	drop contract_ms201`i'3
	gen contract_ms201`i'10=contract_ms201`i'4
	gen contract_ms201`i'11=contract_ms201`i'4
	gen contract_ms201`i'12=contract_ms201`i'4
	drop contract_ms201`i'4
}

*2016 has only 3 quarters
foreach i of numlist 6{
	gen contract_ms201`i'01=contract_ms201`i'1
	gen contract_ms201`i'02=contract_ms201`i'1
	gen contract_ms201`i'03=contract_ms201`i'1
	drop contract_ms201`i'1
	gen contract_ms201`i'04=contract_ms201`i'2
	gen contract_ms201`i'05=contract_ms201`i'2
	gen contract_ms201`i'06=contract_ms201`i'2
	drop contract_ms201`i'2
	gen contract_ms201`i'07=contract_ms201`i'3
	gen contract_ms201`i'08=contract_ms201`i'3
	gen contract_ms201`i'09=contract_ms201`i'3
	drop contract_ms201`i'3
}

reshape long contract_ms, i(firm contract_type) j(month)
********************************************************************************

*clean
order month firm
sort month firm contract_type
drop if contract_type==2 & month>=201301

*******************************************************************************
*creat monthly contract ms - using the macro data
*merge
merge m:1 month firm using `"${PATH_OUT_DATA}/marketshares_long"'
drop if _merge==2
drop _merge
order month firm contract_type

gen mshare_contract=contract_ms*mshare
label var mshare_contract "Supplier market share by contract type (as reported in contract acspt)"
drop  contract_ms mshare
* Save as Stata file.
save `"${PATH_OUT_DATA}/contract_ms"', replace
********************************************************************************
