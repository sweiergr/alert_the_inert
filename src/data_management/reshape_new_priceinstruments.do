/*
	Reshape price data from other regions for use as Hausman IVs. Also, transform nominal into real (Base: jan 2012)
*/

clear
capture log close
version 13
set more off
set mem 4g
* Header do-File with path definitions, those end up in global macros.
include project_paths
log using `"${PATH_OUT_DATA}/log/`1'.log"', replace
********************************************************************************
* Loop over price data for ELEC in other regions.
local region_list wl 
foreach region of local region_list{
	* Create file name.
	local fname_in = `"${PATH_OUT_DATA}/new_pricedata_elec"'+"_"+"`region'"+".dta"
	local fname_out = `"${PATH_OUT_DATA}/new_pricedata_long_elec"'+"_"+"`region'"+".dta"
	* Import data set from /bld/out/data folder.
	use `fname_in', clear
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

	*Drop quantity that are not compatibly with 3500 / 23260 kWh/year
	destring Quantitymin, replace

	*ELEC: *3500 average consumption
	replace Quantitymin=0 if Quantitymin==.
	drop if Quantitymin>3500
	replace Quantitymax=5001 if Quantitymax==.
	drop if Quantitymax<3500
	drop Quantitymin Quantitymax
	**************************
	
	**************************
	*We drop the VARIO contracts from EDF, as they have no contract MS
	drop if contract_type=="delete_noms"
	* Two other have no contract MS, therefore drop them.
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
	* Manual corrections on miscodings.
	* EDF
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

	* Reshape price data into long format.
	reshape long Total, i(firm contract_type) j(month)
	rename Total elec_`region'_price
	label var elec_`region'_price "Simple average of contract price by supplier and contract_type"
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

	rename elec_`region'_price price_`region'_nominal

	gen elec_`region'_price=price_`region'_nominal/CPI
	label var elec_`region'_price "Simple average of real contract price by supplier and contract_type"

	drop CPI CPI_yrl price_`region'_nominal
	*************************

	**************************
	sort month firm 
	order month firm
	* Save as Stata file.
	save `fname_out', replace
}
********************************************************************************

********************************************************************************
* Loop over price data for GAS in all regions.

* We only use data from Flandes and Wallonia, but not from Brussels.
local region_list wl fl
foreach region of local region_list{

	* Create file name.
	local fname_in = `"${PATH_OUT_DATA}/new_pricedata_gas"'+"_"+"`region'"+".dta"
	local fname_out = `"${PATH_OUT_DATA}/new_pricedata_long_gas"'+"_"+"`region'"+".dta"
	* Import data set from /bld/out/data folder.
	use `fname_in', clear
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

	*Drop quantity that are not compatibly with 3500 / 23260 kWh/year
	destring Quantitymin, replace

	*GAS: *23260 average consumption
	replace Quantitymin=0 if Quantitymin==.
	drop if Quantitymin>23260
	replace Quantitymax=1000000 if Quantitymax==.
	drop if Quantitymax<23260
	drop Quantitymin Quantitymax
	**************************

	**************************
	*drop variables that are not needed
	drop Name VarFix ID inms
	drop new_id contract_type //not needed for gas
	**************************

	**************************
	*construct a simple average over contracts by contract_type and by firm
	collapse Total*, by(firm)
	//Do all have entries in each month? -> Yes they do

	* Reshape price data into long format.
	reshape long Total, i(firm) j(month)
	*************************

	**************************
	* Rearrange variables.
	rename Total gas_`region'_price
	label var gas_`region'_price `"Price of gas in `region' (simple average by supplier)"'
	label var month "YearMonth"

	*Encode firm variable
	sort firm month
	encode firm, gen(firm_num)
	drop firm
	rename firm_num firm
	label define firml 1 "ECS" 2 "EDF" 3 "Eneco" 4 "ENINuon" 5 "Essent" 6 "Lampiris" 7 "Other"
	label values firm firml

	sort month firm

	*add numberlabel
	numlabel , add
	**************************

	*************************
	*deflate prices
	merge m:1 month using `"${PATH_OUT_DATA}/CPI_bel"'
	drop if _merge==2
	drop _merge

	rename gas_`region'_price gas_`region'_nominal

	gen gas_`region'_price=gas_`region'_nominal/CPI
	label var gas_`region'_price "Real price of gas in `region' (simple average by supplier)"

	drop CPI CPI_yrl gas_`region'_nominal
	*************************

	sort month firm 
	order month firm
	* Save as Stata file.
	save `fname_out', replace
}
********************************************************************************