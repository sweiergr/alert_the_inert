*******************************************************************************
* This file creates dataset for the micro-moments
*******************************************************************************
* BASICS.
clear
capture log close
version 13
set more off
set mem 4g
* Header do-File with path definitions, those end up in global macros.
include project_paths
log using `"${PATH_OUT_DATA}/log/`1'.log"', replace

******************************************************************************
* Load survey data 
use `"${PATH_OUT_DATA}/survey_controlvars.dta"', clear
*construct labels for contract type -> only 50% were asked
drop if year<=2007 //no question on green (<2007)
drop if year<=2010 //issue in default (<2010)
drop if year==2011 // not needed
drop contract_type
gen contract_type=.
replace contract_type=1 if green_contract==0    
replace contract_type=2 if default==1 & year==2012
replace contract_type=3 if green_contract==1
label define contract_typel 1 "conventional" 2 "default" 3 "green" 4 "green/conv"
label values contract_type contract_typel
*some people got it wrong -> Lampiris and Essent can only be with a green contract
replace contract_type=3 if (firm==6 | firm==3) & !(contract_type==.)
*do not distinguish contract_types in outside option
replace contract_type=1 if firm==7 & !(contract_type==.)
********************************************************************************
* Print some statistics on new relative switching frequency variable.
sum sw_relative*
bysort senior: sum sw_relative* sw_frequency
bysort income: sum sw_relative* sw_frequency
*******************************************************************************
*we only need the following variables in the dataset
keep year firm contract_type income senior vtest sw_relative sw_relative_excl
*clean
order year firm contract_type income senior vtest
sort year firm contract_type income senior vtest
*******************************************************************************

*Drop some observations with missing values
drop if vtest==.
drop if income==. 

*Workaround for dropping "default contract" of EDF.
replace contract_type = 1 if firm==2 & contract_type==2 
numlabel , add
replace contract_type = 1 if contract_type == 2
replace contract_type = 3 if firm==3 | firm==6
replace contract_type = 1 if firm==7
replace contract_type = 0 if contract_type==.
* Compute correlation between using PCW and choosing green.
gen green_pcw = (contract_type==3) * vtest
label var green_pcw "Correlation between using PCW and choosing green contract." 
* Put the sw_relative variable at the end of the dataset 
order year firm contract_type income senior vtest green_pcw sw_relative sw_relative_excl
*comment out the sw_relative_excl for now
drop sw_relative_excl

*******************************************************************************
save `"${PATH_OUT_DATA}/micromoments.dta"', replace
* Export required data to csv for reuse  in MATLAB.
export delimited using `"${PATH_OUT_DATA}/micromoments.csv"', replace nolabel
******************************************************************************
